#include "utils.h"
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <cstdarg>
#include <unordered_map>
#include <vector>
#include "commprof.h"
#include <algorithm>
#include "symbols.h"
#include <inttypes.h>
#include "create_db.h"
#include <iostream>

int prof_enabled = 1;
int local_cid= 0;
int my_coms = 1;
int ac;
char *av[MAX_ARGS];
std::unordered_map<MPI_Request, MPI_Comm> requests_map;
std::vector<prof_attrs*> local_communicators;

/* Tool date */
int mpisee_major_version = 0;
int mpisee_minor_version = 1;
char mpisee_build_date[sizeof(__DATE__)] = __DATE__;
char mpisee_build_time[sizeof(__TIME__)] = __TIME__;
double total_time = 0.0;

extern "C" {
int
namedel(MPI_Comm comm, int keyval, void *attr, void *s)
{
  prof_attrs *com = (prof_attrs*)attr;
  free(com);
  return MPI_SUCCESS;
}
}

extern "C" {
int
namekey(void)
{
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    PMPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}
}

/*
 * Create a new communicator structure and add the parent
 * communicator's name prefix to it
 */
prof_attrs*
get_comm_name(MPI_Comm comm)
{
    int flag;
    prof_attrs *communicator = NULL, *com_info;
    communicator = (prof_attrs*) malloc(sizeof(prof_attrs));
    if (communicator == NULL){
        mcpt_abort("malloc get_comm_name failed\nAborting...\n");
    }
    memset(communicator, 0, sizeof(prof_attrs));
    if ( comm != MPI_COMM_WORLD ){
        PMPI_Comm_get_attr(comm, namekey(), &com_info, &flag);
        if ( flag ){
            strcpy(communicator->name, com_info->name);
        }
        else{
            mcpt_abort("Flag in file:%s line:%d invalid\nAborting\n",
                       __FILE__,__LINE__);
        }
    }
    else{
        strcpy(communicator->name, "W");
    }
    return communicator;
}


void
init_comm(char *buf, prof_attrs** communicator, MPI_Comm comm, MPI_Comm* newcomm){
    size_t length;
    int comm_size,i,j;
    if ( buf == NULL || communicator == NULL ||
         comm == MPI_COMM_NULL || newcomm == NULL){
        mcpt_abort("Newcomm called with NULL\n");
    }
    PMPI_Comm_size(*newcomm, &comm_size);
    length = strlen((*communicator)->name);
    strcpy(&(*communicator)->name[length], buf);
    (*communicator)->size = comm_size;
    for (i = 0; i < NUM_OF_PRIMS; i++) {
        for (j = 0; j < NUM_BUCKETS; j++) {
            (*communicator)->buckets_time[i][j] = 0.0;
            (*communicator)->buckets_msgs[i][j] = 0;
        }

    }
    local_communicators.push_back(*communicator);
    my_coms++;
    local_cid++;
    return;
}

//Change these values modify the buckets
int choose_bucket(int64_t bytes) {
    // These numbers must correspond to the buckets variable
    if (bytes < (1LL << 7)) {       // 1 << 7 is 128
        return 0;
    } else if (bytes < (1LL << 10)) { // 1 << 10 is 1024
        return 1;
    } else if (bytes < (1LL << 13)) { // 1 << 13 is 8192
        return 2;
    } else if (bytes < (1LL << 16)) { // 1 << 16 is 65536
        return 3;
    } else if (bytes < (1LL << 20)) { // 1 << 20 is 1048576
        return 4;
    } else if (bytes < (1LL << 25)) { // 1 << 25 is 33554432
        return 5;
    } else {
        return 6; // For bytes >= 33554432
    }
  // The above is more optimized
  // int index;
  // int64_t tmp;
  // for (index = 0; index < NUM_BUCKETS-1; index++) {
  //   tmp = buckets[index];
  //   tmp = 1 << tmp;
  //   if (tmp > bytes) {
  //       break;
  //   }
  //   if (index == (NUM_BUCKETS - 1) && (bytes > tmp)) {
  //     index += 1;
  //     break;
  //   }
  // }
  // return index;
}



prof_attrs*
profile_this(MPI_Comm comm, int64_t count,MPI_Datatype datatype,int prim,
             double t_elapsed,int root){
    int size,flag,bucket_index;
    prof_attrs *communicator = NULL;
    int64_t sum = 0;
    /* int rank; */
    if ( comm == MPI_COMM_NULL  )
        return communicator;
    flag = 0;
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    if ( datatype != MPI_DATATYPE_NULL ){
        PMPI_Type_size(datatype, &size);
        sum = count * size;
    }
    else{
        sum = count;
    }
    if (flag) {
        bucket_index = choose_bucket(sum);
        communicator->buckets_msgs[prim][bucket_index] += 1;
        communicator->buckets_time[prim][bucket_index] += t_elapsed;
    }
    else{
        fprintf(stderr, "mpisee: empty flag when profiling %s - this might be a bug\n",prim_names[prim]);
    }
    return communicator;
}

int
MPI_Pcontrol(const int level, ...)
{

    int mpi_errno = MPI_SUCCESS;
    va_list list;

    /* ... body of routine ...  */

    va_start(list, level);
    va_end(list);

    if ( level == 1)
        prof_enabled = 1;
    else if ( level == 0 )
        prof_enabled = 0;
    else
        printf("mpisee: MPI_Pcontrol called with invalid value: %d\nProfiling enabled = %d\n",level,prof_enabled);
    /* ... end of body of routine ... */
    return mpi_errno;
}

int
_MPI_Init(int *argc, char ***argv){
    int ret,rank,size;
    int i,j,rc;
    prof_attrs *communicator;
    ret = PMPI_Init(argc, argv);
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);

    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*1024);
        appname = get_appname();
        printf("MPI_Init: mpisee Profiling Tool\nProfiling application\
 %s\n",appname);
 #ifdef MPICH_NAME
        printf("MPICH library used\n");
 #endif
 #ifdef OMPI_MAJOR_VERSION
        printf("OpenMPI library used\n");
 #endif
        fflush(stdout);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    if ( communicator == NULL ){
        mcpt_abort("malloc failed at line %s\n",__LINE__);
    }

    strcpy(communicator->name, "W");
    communicator->size = size;
    for ( i = 0; i<NUM_OF_PRIMS; i++ ){
        for (j = 0; j < NUM_BUCKETS; j++) {
            communicator->buckets_time[i][j] = 0.0;
            communicator->buckets_msgs[i][j] = 0;
        }
    }
    rc = PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);

    local_communicators.push_back(communicator);
    if ( rc != MPI_SUCCESS ){
        mcpt_abort("Comm_set_attr failed at line %s\n",__LINE__);
    }
    if ( argc != NULL )
        ac = *argc;
    total_time = MPI_Wtime();
    return ret;
}


static int
_MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
    int ret,rank,size;
    int i,j,rc;
    prof_attrs *communicator;
    ret = PMPI_Init_thread(argc, argv, required, provided);
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);

    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*1024);
        appname = get_appname();
        printf("MPI_Init_thread: mpisee Profiling Tool\nProfiling\
application %s\n",appname);
        fflush(stdout);
 #ifdef MPICH_NAME
        printf("MPICH library used\n");
 #endif
 #ifdef OMPI_MAJOR_VERSION
        printf("OpenMPI library used\n");
 #endif
        fflush(stdout);
    }
    /* Debuggin file please remove when running */
    /* dname = (char*) malloc (32); */
    /* sprintf(dname, "MCPT_%d", rank); */
    /* dbg_file = fopen(dname,"w"); */
    /********************************************/

    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    if ( communicator == NULL ){
        mcpt_abort("malloc failed at line %s\n",__LINE__);
    }

    strcpy(communicator->name, "W");
    communicator->size = size;
    for ( i = 0; i<NUM_OF_PRIMS; i++ ){
        for (j = 0; j < NUM_BUCKETS; j++) {
            communicator->buckets_time[i][j] = 0.0;
            communicator->buckets_msgs[i][j] = 0;
        }
    }
    local_cid++;
    local_communicators.push_back(communicator);
    rc = PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    if ( rc != MPI_SUCCESS ){
        mcpt_abort("Comm_set_attr failed at line %s\n",__LINE__);
    }
    if ( argc != NULL )
        ac = *argc;
    return ret;
}


int
MPI_Init_thread(int *argc, char ***argv, int required, int *provided)
{
    if ( argc != NULL )
        getProcCmdLine (&ac, av);
    return _MPI_Init_thread(argc, argv, required, provided);
}


extern "C" {
void
F77_MPI_INIT_THREAD (int *required, int *provided, int *ierr)
{
    char **tmp;
    int ret;
    getProcCmdLine (&ac, av);
    tmp = av;
    ret = _MPI_Init_thread(&ac, (char***)&tmp , *required, provided);
    *ierr = ret;
    return;
}
}

/* void
mpi_init_ (int *ierr){
  int ret = 0;
  char **tmp;
  getProcCmdLine (&ac, av);
  tmp = av;
  ret = _MPI_Init (&ac, (char ***) &tmp);
  *ierr = ret;
  return;
} */

extern "C" {
void
F77_MPI_INIT (int *ierr)
{
  int ret = 0;
  char **tmp;
  getProcCmdLine (&ac, av);
  tmp = av;
  ret = _MPI_Init (&ac, (char ***) &tmp);
  *ierr = ret;
  return;
}
}

int
MPI_Init(int *argc, char ***argv)
{
    if ( argc != NULL  ){
        getProcCmdLine (&ac, av);
        //getRunCmd(ac, av);

    }
    return _MPI_Init(argc, argv);
}

/*
 * The naming scheme for the MPI_Comm_create is done by collecting the
 * maximum number of communicators from all processes and using it as id.
 * Use parent communicator name as prefix
 */
int
MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    int ret,r;
    int comms,rank,min_rank;
    char *buf = NULL;
    prof_attrs *communicator;
    ret = PMPI_Comm_create(comm, group, newcomm);
    /*
     * Synchronize with other processess and get the maximum number
     * of communicators to use it an id. This must be done in the
     * parent communicator. Even if Comm_split fails for a process
     * that process must call this Allreduce
     */
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        return ret;
    }
    if ( group == MPI_GROUP_NULL || group == MPI_GROUP_EMPTY )
        return ret;
    PMPI_Comm_rank(comm, &rank);
    PMPI_Allreduce(&rank, &min_rank, 1, MPI_INT, MPI_MIN, *newcomm);
    /* Use the parent's name as a prefix for the newly created communicator */
    communicator = get_comm_name(comm);
    buf = (char *)malloc(sizeof(char) * 16);
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Append prefix+suffix and initialize the data for the new communicator */
    // sprintf(buf, "_c%d.%d", my_coms, min_rank);
    r = snprintf(buf, 16, "_c%d.%d", my_coms, min_rank);
    //truncation error handling
    if (r < 0 || r >= 16) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    init_comm(buf, &communicator, comm, newcomm);
    free(buf);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}


extern "C" {
void
F77_MPI_COMM_CREATE(MPI_Fint  * comm, MPI_Fint  * group, MPI_Fint  *comm_out ,
                    MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm, c_comm_out;
    MPI_Group c_group;
    c_comm = MPI_Comm_f2c(*comm);
    c_group = MPI_Group_f2c(*group);

    ret= MPI_Comm_create(c_comm, c_group, &c_comm_out);
    *ierr = ret;
    if( ret == MPI_SUCCESS )
        *comm_out = MPI_Comm_c2f(c_comm_out);
    return;

}
}

/*
 * The naming scheme for the MPI_Comm_split requires two numbers
 * 1. current number of communicators
 * 2. a unique id for each new communicator via split
 * The prefix of this name is always the parent's name
 */
int
MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    int ret,comms,r;
    prof_attrs *communicator;
    char *buf = NULL;
    int rank, min_rank;
    /* Call the actual split */
    ret = PMPI_Comm_split(comm, color, key, newcomm);
    /*
     * Synchronize with other processess and get the maximum number
     * of communicators to use it an id. This must be done in the
     * parent communicator. Even if Comm_split fails for a process
     * that process must call this Allreduce
     */
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    if ( *newcomm == MPI_COMM_NULL  ){
        return ret;
    }
    PMPI_Comm_rank(comm, &rank);
    /*
     * Get the minimum rank. Note that the rank is the rank from the parent
     * communicator. Every new communicator will have an id based on the
     * minimum rank of a process in the parent communicator
     */
    PMPI_Allreduce(&rank, &min_rank, 1, MPI_INT, MPI_MIN, *newcomm);
    /* Use the parent's name as a prefix for the newly created communicator */
    communicator = get_comm_name(comm);
    buf = (char *)malloc(sizeof(char) * 16);
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 16, "_s%d.%d", my_coms, min_rank);
    //truncation error handling
    if (r < 0 || r >= 16) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    // sprintf(buf, "_s%d.%d", my_coms, min_rank);
    /* Append prefix+suffix and initialize the data for the new communicator */
    init_comm(buf, &communicator, comm, newcomm);
    free(buf);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}

extern "C" {
void
F77_MPI_COMM_SPLIT(MPI_Fint  * comm, int  * color, int  * key,
                   MPI_Fint  *comm_out , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm,c_comm_out;
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_split(c_comm, *color, *key, &c_comm_out);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *comm_out = MPI_Comm_c2f(c_comm_out);
    return;
}
}


int
MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
    int ret,comms,r;
    prof_attrs *communicator;
    char *buf = NULL;
    ret = PMPI_Comm_dup(comm, newcomm);
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;
    communicator = get_comm_name(comm);
    buf = (char *)malloc(sizeof(char) * 8);
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 8, "_d%d", my_coms);
    //truncation error handling
    if (r < 0 || r >= 8) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    /* Append prefix+suffix and initialize the data for the new communicator */
    // sprintf(buf,"_d%d",my_coms);
    init_comm(buf, &communicator, comm, newcomm);
    free(buf);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}


extern "C" {
void
F77_MPI_COMM_DUP(MPI_Fint  * comm, MPI_Fint  *comm_out , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm newcomm;
    MPI_Comm c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_dup(c_comm, &newcomm);
    *ierr = ret;
    if ( ret == MPI_SUCCESS  )
        *comm_out =  MPI_Comm_c2f(newcomm);
    return;
}
}

/* idup is supposed to be non-blocking but we make blocking calls in the wrapper */
int
MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request)
{

    int ret,comms,r;
    prof_attrs *communicator;
    char *buf=NULL;
    ret = PMPI_Comm_idup(comm, newcomm, request);

// #ifdef OMPI_MAJOR_VERSION
    requests_map[*request] = comm;
// #endif
// #ifdef MPICH_NAME
    // requests_map[*request] = comm;
// #endif

    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;
    communicator = get_comm_name(comm);
    buf = (char *)malloc(sizeof(char) * 8);
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 8, "_i%d", my_coms);
    //truncation error handling
    if (r < 0 || r >= 8) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    /* Append prefix+suffix and initialize the data for the new communicator */
    // sprintf(buf,"_i%d",my_coms);
    init_comm(buf, &communicator, comm, newcomm);

    free(buf);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}

/* TODO idup wrapper F77 */

int
MPI_Cart_create(MPI_Comm old_comm, int ndims, const int *dims,
                const int *periods, int reorder, MPI_Comm *comm_cart)
{
    int ret,comms,r;
    prof_attrs *communicator;
    char *buf = NULL;
    ret = PMPI_Cart_create(old_comm, ndims, dims, periods, reorder, comm_cart);
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, old_comm);
    my_coms = comms;
    /* Should we have an if condition here to check if comm_cart is null? */
    communicator = get_comm_name(old_comm);
    buf = (char *)malloc(8 * sizeof(char));
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 8, "_a%d", my_coms);
    //truncation error handling
    if (r < 0 || r >= 8) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    /* Append prefix+suffix and initialize the data for the new communicator */

    // sprintf(buf,"_a%d",my_coms);
    init_comm(buf, &communicator, old_comm, comm_cart);
    free(buf);
    PMPI_Comm_set_attr(*comm_cart, namekey(), communicator);
    return ret;
}


extern "C" {
void
F77_MPI_CART_CREATE(MPI_Fint  * comm_old, int  * ndims, const int  *dims,
                    const int  *periods, int  * reorder,
                    MPI_Fint  *comm_cart , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm_old, c_comm_cart;
    c_comm_old = MPI_Comm_f2c(*comm_old);
    ret = MPI_Cart_create(c_comm_old, *ndims, dims, periods, *reorder, &c_comm_cart);
    *ierr = ret;
    if ( ret == MPI_SUCCESS  )
        *comm_cart = MPI_Comm_c2f(c_comm_cart);
    return;
}
}


int
MPI_Cart_sub(MPI_Comm comm, const int *remain_dims, MPI_Comm *new_comm)
{
    int ret, my_rank,r;
    prof_attrs *communicator;
    char *buf=NULL;
    int id,min_rank;
    ret = PMPI_Cart_sub(comm, remain_dims, new_comm);
    PMPI_Allreduce(&my_coms, &id, 1, MPI_INT, MPI_MAX, comm);
    if ( new_comm == NULL || *new_comm == MPI_COMM_NULL ){
        return ret;
    }
    PMPI_Comm_rank(comm, &my_rank);
    PMPI_Allreduce(&my_rank, &min_rank, 1, MPI_INT , MPI_MIN, *new_comm);
    my_coms = id;
    communicator = get_comm_name(comm);
    buf = (char *)malloc(sizeof(char) * 16);

    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 16, "_b%d.%d", my_coms,min_rank);
    //truncation error handling
    if (r < 0 || r >= 16) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    // sprintf(buf,"_b%d.%d",id,min_rank);
    init_comm(buf, &communicator, comm, new_comm);
    PMPI_Comm_set_attr(*new_comm, namekey(), communicator);
    return ret;
}

extern "C" {
void
F77_MPI_CART_SUB(MPI_Fint  * comm, const int  *remain_dims,
                 MPI_Fint  *comm_new , MPI_Fint *ierr)
{
    int rc;
    MPI_Comm c_comm;
    MPI_Comm c_comm_new;
    c_comm = MPI_Comm_f2c(*comm);
    c_comm_new = MPI_Comm_f2c(*comm_new);
    rc = MPI_Cart_sub(c_comm, remain_dims, &c_comm_new);
    *comm_new = MPI_Comm_c2f(c_comm_new);
    *ierr = rc;
    return;
}
}

int
MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int *index,
                 const int *edges, int reorder, MPI_Comm *comm_graph)
{
    int ret,comms,r;
    prof_attrs *communicator;
    char *buf=NULL;

    ret = PMPI_Graph_create(comm_old, nnodes, index, edges, reorder, comm_graph);

    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm_old);
    my_coms = comms;

    communicator = get_comm_name(comm_old);
    buf = (char *)malloc(8 * sizeof(char));
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 8, "_r%d", my_coms);
    //truncation error handling
    if (r < 0 || r >= 8) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    /* Append prefix+suffix and initialize the data for the new communicator */

    // sprintf(buf,"_r%d",my_coms);
    init_comm(buf, &communicator, comm_old, comm_graph);
    free(buf);

    PMPI_Comm_set_attr(*comm_graph, namekey(), communicator);

    return ret;
}


extern "C" {
void
F77_MPI_GRAPH_CREATE(MPI_Fint  * comm_old, int  * nnodes, const int  *index,
                     const int  *edges, int  * reorder, MPI_Fint  *comm_graph,
                     MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm_old;
    MPI_Comm c_comm_graph;

    c_comm_old = MPI_Comm_f2c(*comm_old);

    ret = MPI_Graph_create(c_comm_old, *nnodes, index, edges, *reorder, &c_comm_graph);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS ) {
        *comm_graph = MPI_Comm_c2f(c_comm_graph);
    }
    return;
}

}


int
MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int *nodes,
                      const int *degrees, const int *targets,
                      const int *weights, MPI_Info info, int reorder,
                      MPI_Comm *newcomm)
{
    int ret,comms,r;
    prof_attrs *communicator;
    char *buf=NULL;
    ret = PMPI_Dist_graph_create(comm_old, n, nodes, degrees, targets, weights, info, reorder, newcomm);

    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm_old);
    my_coms = comms;

    communicator = get_comm_name(comm_old);
    buf = (char *)malloc(8 * sizeof(char));
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 8, "_g%d", my_coms);
    //truncation error handling
    if (r < 0 || r >= 8) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    /* Append prefix+suffix and initialize the data for the new communicator */

    // sprintf(buf,"_g%d",my_coms);
    init_comm(buf, &communicator, comm_old, newcomm);
    free(buf);

    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}


int
MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info,
                    MPI_Comm *newcomm){
    int ret,comms,r;
    prof_attrs *communicator;
    char *buf=NULL;
    int rank, min_rank;

    ret = PMPI_Comm_split_type(comm, split_type, key, info, newcomm);

    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    PMPI_Comm_rank(comm, &rank);
    if ( newcomm== NULL || *newcomm == MPI_COMM_NULL  ){
        return ret;
    }
    PMPI_Allreduce(&rank, &min_rank, 1, MPI_INT, MPI_MIN, *newcomm);
    communicator = get_comm_name(comm);
    buf = (char *)malloc(sizeof(char) * 16);
    if (buf == NULL) {
      mcpt_abort("Malloc failed\n");
      return 1;

    }
    /* Suffix of the new communicator with the two ids */
    r = snprintf(buf, 16, "_b%d.%d", my_coms,min_rank);
    //truncation error handling
    if (r < 0 || r >= 16) {
      mcpt_abort("Snprintf truncation error\n");
      return 2;
    }
    // sprintf(buf,"_t%d.%d",comms,min_rank);
    init_comm(buf, &communicator, comm, newcomm);
    free(buf);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}

extern "C" {
void
F77_MPI_COMM_SPLIT_TYPE(MPI_Fint  * comm, int  * split_type, int  * key,
                        MPI_Fint *info, MPI_Fint  *newcomm , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm,c_comm_out;
    MPI_Info c_info;
    c_info = PMPI_Info_f2c(*info);
    c_comm = PMPI_Comm_f2c(*comm);
    ret = MPI_Comm_split_type(c_comm, *split_type, *key, c_info,&c_comm_out);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *newcomm = MPI_Comm_c2f(c_comm_out);
    return;
}
}

int
MPI_Isend(const void *buf, int count, MPI_Datatype datatype,int dest, int tag,
          MPI_Comm comm, MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled  == 1){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm, count, datatype, Isend, t_elapsed, 0);
// #ifdef OMPI_MAJOR_VERSION
        requests_map[*request] = comm;
// #endif
// #ifdef MPICH_NAME
//         requests_map[*request] = comm;
// #endif

    }
    else{
        ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    }
    return ret;
}

extern "C" {
void
F77_MPI_ISEND(const void  *buf, int  * count, MPI_Fint  * datatype,
                          int  * dest, int  * tag, MPI_Fint  * comm,
                          MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Isend(buf, *count, c_datatype, *dest, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}


extern "C" {
int
MPI_Send(const void *buf, int count, MPI_Datatype datatype,
         int dest,int tag, MPI_Comm comm)
{
    int ret;
    double t_elapsed = 0.0;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count, datatype, Send, t_elapsed, 0);
    }
    else{
        ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    }
    return ret;
}
}

extern "C" {
void
F77_MPI_SEND(const void  *buf, int  * count, MPI_Fint  * datatype,
             int  * dest, int  * tag, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Send(buf, *count, c_datatype, *dest, *tag, c_comm);
    *ierr = ret;
    return;
}
}

extern "C"{
int
MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
          MPI_Comm comm, MPI_Request *request)
{

    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm, count,datatype,Irecv,t_elapsed,0);
// #ifdef OMPI_MAJOR_VERSION
        requests_map[*request] = comm;
// #endif
// #ifdef MPICH_NAME
//         requests_map[request] = comm;
// #endif
    }
    else{
        ret = PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
    }
    return ret;
}
}

extern "C" {
void
F77_MPI_IRECV(void  *buf, int  * count, MPI_Fint  * datatype, int  * source,
              int  * tag, MPI_Fint  * comm, MPI_Fint  *request , MPI_Fint *ierr)
{

    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Irecv(buf, *count, c_datatype, *source, *tag, c_comm, &c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Recv(buf,count,datatype,source,tag,comm,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,count,datatype,Recv,t_elapsed,0);
    }
    else{
        ret = PMPI_Recv(buf,count,datatype,source,tag,comm,status);
    }
    return ret;
}

extern "C" {
void
F77_MPI_RECV(void* buf, int* count,MPI_Fint* datatype, int* source, int* tag,
             MPI_Fint  * comm, MPI_Status  *status , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Recv(buf,*count,c_datatype,*source,*tag,c_comm,status);
    *ierr = ret;
    return;
}
}

int
MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             int dest, int sendtag, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int source, int recvtag,
             MPI_Comm comm, MPI_Status *status)
{

    int ret;
    double t_elapsed;
    int64_t sum;

    if ( prof_enabled == 1){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                            recvcount, recvtype, source, recvtag, comm, status);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = sendcount;
        sum = sum | 0x1;
        sum = sum>>1;
        profile_this(comm,sum,sendtype,Sendrecv,t_elapsed,source);
    }
    else{
        ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                        recvcount, recvtype, source, recvtag, comm, status);
    }

    return ret;
}


extern "C" {
void
F77_MPI_SENDRECV(const void  *sendbuf, int  * sendcount,MPI_Fint  * sendtype,
                 int  * dest, int  * sendtag, void  *recvbuf, int  * recvcount,
                 MPI_Fint  * recvtype, int  * source, int  * recvtag,
                 MPI_Fint  * comm, MPI_Status  *status , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;
    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Sendrecv(sendbuf, *sendcount, c_sendtype, *dest, *sendtag, recvbuf,
                        *recvcount, c_recvtype, *source, *recvtag, c_comm, status);
    *ierr = ret;
    return;
}
}


int
MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
          MPI_Comm comm)
{
    int ret,rank,sum;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Bcast(buffer, count, datatype, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            sum = 0;
        }
        else {
            sum = count;
        }
        profile_this(comm,sum,datatype,Bcast,t_elapsed,root);
    }
    else{
        ret = PMPI_Bcast(buffer, count, datatype, root, comm);
    }
    return ret;

}

extern "C" {
void
F77_MPI_BCAST(void  *buffer, int  * count, MPI_Fint  * datatype, int  * root,
              MPI_Fint  * comm , MPI_Fint *ierr)
{

    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Bcast(buffer, *count, c_datatype, *root, c_comm);
    *ierr = ret;
}
}

int
MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root,
           MPI_Comm comm, MPI_Request *request)
{
    int ret,rank,sum;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Ibcast(buffer, count, datatype, root, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            sum = 0;
        }
        else {
            sum = count;
        }
        profile_this(comm,sum,datatype,Ibcast,t_elapsed,root);
    }
    else{

        ret = PMPI_Ibcast(buffer, count, datatype, root, comm, request);
    }
    return ret;
}

extern "C" {
void
F77_MPI_IBCAST(void  *buffer, int  * count, MPI_Fint  * datatype, int  * root,
               MPI_Fint  * comm, MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Comm c_comm;
    MPI_Request c_request;
    c_datatype = MPI_Type_f2c(*datatype);
    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Ibcast(buffer, *count, c_datatype, *root, c_comm,&c_request);
    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
}
}

int
MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ret;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Allreduce,t_elapsed,0);
    }
    else{
        ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    }
    return ret;
}

extern "C" {
void
F77_MPI_ALLREDUCE(const void  *sendbuf, void  *recvbuf, int  * count,
                  MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm ,
                  MPI_Fint *ierr)
{
    int ret;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Datatype c_datatype;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Allreduce(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm);
    *ierr =ret;
    return;
}
}

int
MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
               MPI_Request *request)
{
    int ret;
    double t_elapsed;
    if (prof_enabled == 1){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm, count, datatype, Iallreduce, t_elapsed, 0);
    }
    else{
        ret = PMPI_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, request);
    }
    return ret;
}


extern "C" {
void
F77_MPI_IALLREDUCE(const void  *sendbuf, void  *recvbuf, int  * count,
                   MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm,
                   MPI_Fint  *request , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;
    MPI_Request c_request;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Iallreduce(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm, &c_request);

    *ierr = ret;
    if ( ret == MPI_SUCCESS )
        *request = MPI_Request_c2f(c_request);
    return;
}
}

int
MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,recvcount,recvtype,Allgather,t_elapsed,0);
    }
    else{
        ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                             recvtype, comm);
    }
    return ret;
}


extern "C" {
void
F77_MPI_ALLGATHER(const void *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                  void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype,
                  MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;
    MPI_Datatype c_sendtype, c_recvtype;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Allgather(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcount, c_recvtype, c_comm);
    *ierr =ret;
    return;
}
}

int
MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed =  MPI_Wtime();
        ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, comm);

        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,sendcount,sendtype,Alltoall,t_elapsed,0);
    }
    else{
        ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                            recvtype, comm);
    }
    return ret;
}

extern "C" {
void
F77_MPI_ALLTOALL(const void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                 void  *recvbuf, int  * recvcnt, MPI_Fint  * recvtype,
                 MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Alltoall(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcnt, c_recvtype, c_comm);
    *ierr = ret;
    return;
}
}


int
MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
              const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
              const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm)
{
    int ret,sum,i,sz;
    double t_elapsed;
    sum = 0;
    t_elapsed = MPI_Wtime();
    if ( prof_enabled == 1 ){
        ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts,
                             rdispls, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if ( sendcounts[i] > 0 )
                sum+=sendcounts[i];
        }
        /* We won't need this reduce just sum all in the end */
        /* PMPI_Reduce(&sum, &sum_max, 1, MPI_INT, MPI_MAX, 0, comm); */
        profile_this(comm,sum,sendtype,Alltoallv,t_elapsed,0);
    }
    else{
        ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts,
                             rdispls, recvtype, comm);
    }
    return ret;
}


extern "C" {
void
F77_MPI_ALLTOALLV(const void  *sendbuf, const int  *sendcnts, const int  *sdispls,
                  MPI_Fint  * sendtype, void  *recvbuf, const int  *recvcnts,
                  const int  *rdispls, MPI_Fint  * recvtype, MPI_Fint  * comm ,
                  MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Alltoallv(sendbuf, sendcnts, sdispls, c_sendtype, recvbuf, recvcnts, rdispls, c_recvtype, c_comm);
    *ierr = ret;
    return;
}
}

int
MPI_Alltoallw(const void *sendbuf, const int *sendcounts, const int *sdispls,
              const MPI_Datatype *sendtypes, void *recvbuf, const int *recvcounts,
              const int *rdispls, const MPI_Datatype *recvtypes, MPI_Comm comm)
{
    int ret,sum,i,sz,type_sz;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        sum = 0;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if ( sendcounts[i] > 0 ){
                PMPI_Type_size(sendtypes[i], &type_sz);
                sum+=(sendcounts[i]*type_sz);
            }
        }
        /* We won't need this reduce just sum all in the end */
        /* PMPI_Reduce(&sum, &sum_max, 1, MPI_INT, MPI_MAX, 0, comm); */
        profile_this(comm,sum,MPI_DATATYPE_NULL,Alltoallw,t_elapsed,0);
    }
    else{
        ret = PMPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
    }
    return ret;
}

/* TODO F77 Wrapper Alltoallw */

int
MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, const int *recvcounts, const int *displs,
               MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret,sum,i,sz;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                              displs, recvtype, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = 0;
        /* tmp = recvcounts; */
        /* while(tmp){ */
        /*     sum += *tmp; */
        /*     tmp++; */
        /* } */
        MPI_Comm_size(comm, &sz);
        for ( i=0; i<sz; i++ ){
            if( recvcounts[i] > 0 )
                sum+=recvcounts[i];
        }

        profile_this(comm,sum,recvtype,Allgatherv,t_elapsed,0);
    }
    else{
        ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                              displs, recvtype, comm);
    }
    return ret;
}


extern "C" {
void
F77_MPI_ALLGATHERV(const void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                   void  *recvbuf, int  *recvcounts, int  *displs,
                   MPI_Fint  * recvtype, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Allgatherv(sendbuf, *sendcount, c_sendtype, recvbuf, recvcounts, displs, c_recvtype, c_comm);
    *ierr = ret;
    return;
}
}


int
MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
           MPI_Op op, int root, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();

        ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Reduce,t_elapsed,root);
    }
    else{
        ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    }
    return ret;
}

extern "C" {
void F77_MPI_REDUCE(const void  *sendbuf, void  *recvbuf, int  * count,
                    MPI_Fint  * datatype, MPI_Fint  * op, int  * root,
                    MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Reduce(sendbuf, recvbuf, *count, c_datatype, c_op, *root, c_comm);
    *ierr = ret;
    return;

}
}

int
MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
           int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);

        t_elapsed = MPI_Wtime() - t_elapsed;
        profile_this(comm,sendcount,sendtype,Gather,t_elapsed,root);
    }
    else{
        ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    }
    return ret;

}

extern "C" {
void
F77_MPI_GATHER(const void  *sendbuf, int  * sendcnt, MPI_Fint  * sendtype,
               void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype,
               int  * root, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Gather(sendbuf, *sendcnt, c_sendtype, recvbuf, *recvcount, c_recvtype, *root, c_comm);
    *ierr = ret;
    return;

}
}


int
MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, const int *recvcounts, const int *displs,
            MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ret,sum,rank,comm_size;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                           recvtype, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        sum = 0;
        /* tmp = recvcounts; */
        /* while(tmp){ */
        /*     sum += *tmp; */
        /*     tmp++; */
        /* } */
        /* PMPI_Reduce(&sum, &max_sum, 1, MPI_INT, MPI_MAX, 0, comm); */

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            PMPI_Comm_size(comm, &comm_size);
            for (int i = 0; i < comm_size; i++){
                sum += recvcounts[i];
            }
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, recvtype, Gatherv, t_elapsed, root);
    }
    else{
        ret = PMPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                           recvtype, root, comm);
    }

    return ret;
}

extern "C" {
void
F77_MPI_GATHERV(const void  *sendbuf, int* sendcount, MPI_Fint  * sendtype,
                void  *recvbuf, const int* recvcounts, const int  *displs,
                MPI_Fint  * recvtype, int  * root, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Gatherv(sendbuf, *sendcount, c_sendtype, recvbuf, recvcounts,
                      displs, c_recvtype, *root, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}


int
MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
             MPI_Datatype sendtype, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int root, MPI_Comm comm)
{

    int ret,rank,comm_size;
    uint64_t sum;
    double t_elapsed;
    if ( prof_enabled == 1 ){
        sum = 0;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                            recvtype, root, comm);

        t_elapsed = MPI_Wtime() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            PMPI_Comm_size(comm, &comm_size);
            for (int i = 0; i < comm_size; i++){
                sum += sendcounts[i];
            }
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, sendtype, Scatterv, t_elapsed, root);
    }
    else{
        ret = PMPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                            recvtype, root, comm);
    }
    /* if ( rank == root ){ */
    /*     communicator->bytes += sum; */
    /*     communicator->prim_bytes[Scatterv] += sum; */
    /* } */
    return ret;
}

extern "C" {
void
F77_MPI_SCATTERV(const void  *sendbuf, const int  *sendcounts, const int  *displs,
                 MPI_Fint  * sendtype, void  *recvbuf, int  * recvcount,
                 MPI_Fint  * recvtype, int  * root, MPI_Fint  * comm, MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Scatterv(sendbuf, sendcounts, displs, c_sendtype, recvbuf, *recvcount,
                 c_recvtype, *root, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
            void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
            MPI_Comm comm)
{
    int ret,rank,sum;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, root, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        PMPI_Comm_rank(comm, &rank);
        if ( rank == root ){
            sum = sendcount;
        }
        else{
            sum = 0;
        }
        profile_this(comm, sum, sendtype, Scatter, t_elapsed, root);
    }
    else{
        ret = PMPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                           recvtype, root, comm);
    }
    return ret;
}

extern "C" {
void
F77_MPI_SCATTER(const void  *sendbuf, int  * sendcount, MPI_Fint  * sendtype,
                void  *recvbuf, int  * recvcount, MPI_Fint  * recvtype, int  * root,
                MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_sendtype;
    MPI_Datatype c_recvtype;
    MPI_Comm c_comm;

    c_sendtype = MPI_Type_f2c(*sendtype);
    c_recvtype = MPI_Type_f2c(*recvtype);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Scatter(sendbuf, *sendcount, c_sendtype, recvbuf, *recvcount,
                      c_recvtype, *root, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}


int
MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
         MPI_Op op, MPI_Comm comm)
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
        t_elapsed = MPI_Wtime() - t_elapsed;

        profile_this(comm,count,datatype,Scan,t_elapsed,0);
    }
    else{
        ret = PMPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
    }
    return ret;

}


extern "C" {
void
F77_MPI_SCAN(const void  *sendbuf, void  *recvbuf, int  * count, MPI_Fint  * datatype,
             MPI_Fint  * op, MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = PMPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Scan(sendbuf, recvbuf, *count, c_datatype, c_op, c_comm);

    *ierr = (MPI_Fint)ret;
    return;
}
}


int
MPI_Barrier ( MPI_Comm comm )
{
    int ret;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Barrier(comm);
        t_elapsed = MPI_Wtime()-t_elapsed;

        profile_this(comm,0,MPI_DATATYPE_NULL,Barrier,t_elapsed,0);
    }
    else{
        ret = PMPI_Barrier(comm);
    }
    return ret;
}


extern "C" {
void
F77_MPI_BARRIER(MPI_Fint  * comm , MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Barrier(c_comm);
    *ierr = (MPI_Fint)ret;
    return;
}
}

int
MPI_Wait(MPI_Request *request, MPI_Status *status)
{
    int ret;
    double t_elapsed;
    /* int i; */

// #ifndef MPICH_API_PUBLIC
    MPI_Comm comm ;
// #else
//     MPI_Comm *comm = NULL;
// #endif
    if ( prof_enabled == 1 ){
        //comm = Table_get(request_tab, request);
        comm = requests_map[*request];
        t_elapsed = MPI_Wtime();
        ret = PMPI_Wait(request, status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        //PMPI_Comm_rank(comm,&rank);
        if ( comm == NULL  ){
            fprintf(stderr, "MCPT: NULL COMMUNICATOR in MPI_Wait\n"); 
            return ret;
        }
// #ifdef OMPI_MAJOR_VERSION
       profile_this(comm, 0, MPI_DATATYPE_NULL, Wait, t_elapsed, 0);
// #endif
// #ifdef MPICH_NAME
//         profile_this(*comm, 0, MPI_DATATYPE_NULL, Wait, t_elapsed, 0);
// #endif
    }
    else{
        ret = PMPI_Wait(request, status);
    }
    return ret;
}

extern "C" {
void
F77_MPI_WAIT(MPI_Fint  *request, MPI_Status  *status , MPI_Fint *ierr)
{
   int ret;
   MPI_Request c_request;
   c_request = MPI_Request_f2c(*request);
   ret = MPI_Wait(&c_request, status);
   *ierr = ret;
}
}


int
MPI_Waitall(int count, MPI_Request array_of_requests[],
            MPI_Status array_of_statuses[])
{
    int ret,i;
    double t_elapsed;
// #ifndef MPICH_API_PUBLIC
    MPI_Comm comm = NULL;
// #else
//     MPI_Comm *comm = NULL;
// #endif

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
        t_elapsed = MPI_Wtime() - t_elapsed;
        for ( i =0; i<count; i++ ){
            if ( comm != NULL ){
                comm = requests_map[array_of_requests[i]];
                profile_this(comm, 0, MPI_DATATYPE_NULL, Waitall, t_elapsed, 0);
                return ret;
            }
        }
    }
    else{
        ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
    }
    return ret;
}


extern "C" {
void
F77_MPI_WAITALL(int  * count, MPI_Fint  *array_of_requests,
                MPI_Status  *array_of_statuses , MPI_Fint *ierr)
{

   int ret,i;
   MPI_Request *c_requests;
   c_requests = (MPI_Request*) malloc (sizeof(MPI_Request)*(*count));
   for ( i =0; i<*count; i++ ){
       c_requests[i] = MPI_Request_f2c(array_of_requests[i]);
   }
   ret = MPI_Waitall(*count, c_requests, array_of_statuses);
   *ierr = ret;
   if ( ret == MPI_SUCCESS ){
       for ( i =0; i<*count; i++ ){
           array_of_requests[i] = MPI_Request_c2f(c_requests[i]);
       }
   }
   free( c_requests );
}
}

int
MPI_Waitany(int count, MPI_Request *array_of_requests, int *index, MPI_Status *status)
{
    int ret,i;
    double t_elapsed;
// #ifndef MPICH_API_PUBLIC
    MPI_Comm comm = NULL;
// #else
//     MPI_Comm *comm = NULL;
// #endif

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitany(count, array_of_requests,index,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        for ( i =0; i<count; i++ ){
            comm = requests_map[array_of_requests[i]];
            if ( comm != NULL ){
                profile_this(comm, 0, MPI_DATATYPE_NULL, Waitany, t_elapsed, 0);
                break;
            }
        }
    }
    else{
        ret = PMPI_Waitany(count, array_of_requests,index,status);
    }
    return ret;
}

extern "C" {
void
F77_MPI_WAITANY(int  * count, MPI_Fint  *array_of_requests, int  *index,
                MPI_Status  *status , MPI_Fint *ierr)
{
    int ret,i;
    MPI_Request *c_array_of_requests = NULL;
    //c_array_of_requests = (MPI_Request*)malloc(sizeof(MPI_Request)*(*count));
    //assert(c_array_of_requests);
    for (i = 0; i < *count; i++) {
        c_array_of_requests[i] = MPI_Request_f2c(array_of_requests[i]);
    }

    ret = MPI_Waitany(*count, c_array_of_requests, index, status);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS ) {
        array_of_requests[*index] = MPI_Request_c2f(c_array_of_requests[*index]);
        if ( *index >= 0 ) (*index)++;
    }
    free(c_array_of_requests);
    return;

}
}

int
MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)
{

    int ret;
    double t_elapsed;
    /* int i; */

// #ifndef MPICH_API_PUBLIC
    MPI_Comm comm = NULL;
// #else
//     MPI_Comm *comm = NULL;
// #endif
    if ( prof_enabled == 1 ){
        comm = requests_map[*request];
        t_elapsed = MPI_Wtime();
        ret = PMPI_Test(request,flag,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        //comm = Table_get(request_tab, request);
        if ( comm == NULL ){
            /* fprintf(stderr, "MCPT: NULL COMMUNICATOR in MPI_Wait\n"); */
            return ret;
        }
// #ifndef MPICH_API_PUBLIC
        profile_this(comm, 0, MPI_DATATYPE_NULL, Test, t_elapsed, 0);
// #else
//         profile_this(*comm, 0, MPI_DATATYPE_NULL, Test, t_elapsed, 0);
// #endif
    }
    else{
        ret = PMPI_Test(request,flag,status);
    }
    return ret;
}


extern "C" {
void
F77_MPI_TEST(MPI_Fint  *request, int  *flag, MPI_Status  *status , MPI_Fint *ierr)
{
    int ret;
    MPI_Request c_request;

    c_request = MPI_Request_f2c(*request);
    ret = MPI_Test(&c_request,flag,status);
    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS ) {
        *request = MPI_Request_c2f(c_request);
    }
    return;
}
}

int
MPI_Testany(int count, MPI_Request *array_of_requests, int *index, int *flag, MPI_Status *status)
{
    int ret,i;
    double t_elapsed;
// #ifndef MPICH_API_PUBLIC
    MPI_Comm comm = NULL;
// #else
//     MPI_Comm *comm = NULL;
// #endif

    if ( prof_enabled == 1 ){
        t_elapsed = MPI_Wtime();
        ret = PMPI_Testany(count, array_of_requests, index, flag, status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        for ( i =0; i<count; i++ ){
            comm = requests_map[array_of_requests[i]];
            if ( comm != NULL ){
                profile_this(comm, 0, MPI_DATATYPE_NULL, Waitany, t_elapsed, 0);
                break;
            }
        }
    }
    else{
        ret = PMPI_Testany(count, array_of_requests, index, flag, status);
    }
    return ret;

}


extern "C" {
void
F77_MPI_TESTANY(int  * count, MPI_Fint  *array_of_requests, int  *index,
                int  *flag, MPI_Status  *status , MPI_Fint *ierr)
{

    int ret,i;
    MPI_Request *c_array_of_requests = NULL;
    //c_array_of_requests = (MPI_Request*)malloc(sizeof(MPI_Request)*(*count));
    //assert(c_array_of_requests != NULL);
    for (i = 0; i < *count; i++) {
        c_array_of_requests[i] = MPI_Request_f2c(array_of_requests[i]);
    }

    ret = MPI_Testany(*count, c_array_of_requests, index, flag, status);

    *ierr = (MPI_Fint)ret;
    if ( ret == MPI_SUCCESS ) {
        array_of_requests[*index] = MPI_Request_c2f(c_array_of_requests[*index]);
        if ( *index >= 0 ) (*index)++;
    }
    free(c_array_of_requests);
    return;

}
}

int
MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int *recvcounts,
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{
    int ret,sum,rank;
    const int *cnt;
    double t_elapsed;

    if ( prof_enabled == 1 ){
        sum = 0;
        cnt = recvcounts;
        t_elapsed = MPI_Wtime();
        ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
        t_elapsed = MPI_Wtick() - t_elapsed;
        PMPI_Comm_rank(comm, &rank);
        while(cnt){
            sum += *cnt;
        }
        sum += recvcounts[rank];
        profile_this(comm, sum, datatype, Reduce_scatter, t_elapsed, 0);
    }
    else{
        ret = PMPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm);
    }
    return ret;

}

extern "C" {
void F77_MPI_REDUCE_SCATTER(const void  *sendbuf, void  *recvbuf, const int *recvcnts,
                            MPI_Fint  * datatype, MPI_Fint  * op, MPI_Fint  * comm,
                            MPI_Fint *ierr)
{
    int ret;
    MPI_Datatype c_datatype;
    MPI_Op c_op;
    MPI_Comm c_comm;

    c_datatype = MPI_Type_f2c(*datatype);
    c_op = MPI_Op_f2c(*op);
    c_comm = MPI_Comm_f2c(*comm);

    ret = MPI_Reduce_scatter(sendbuf, recvbuf, recvcnts, c_datatype, c_op, c_comm);
    *ierr = (MPI_Fint)ret;
    return;

}
}

int
MPI_Comm_free(MPI_Comm *comm)
{
    int ret,flag;
    prof_attrs *com_info, *tmp;
    PMPI_Comm_get_attr(*comm, namekey(), &tmp, &flag);
    if (flag) {
        auto it = std::find(local_communicators.begin(), local_communicators.end(), tmp);
        if (it != local_communicators.end()) {
            // Calculate the index of the found element
            size_t index = std::distance(local_communicators.begin(), it);
            // Now "index" holds the position of com_info in the vector
            com_info = (prof_attrs *)malloc(sizeof(prof_attrs));
            memcpy(com_info, tmp, sizeof(prof_attrs));
            local_communicators[index] = com_info;
        } else {
            // com_info is not present in the vector
            mcpt_abort("Comm_free on invalid communicator\n");
        }

    }else {
        mcpt_abort("Comm free: Comm_get_attr did not find communicator\n");
    }
    ret = PMPI_Comm_free(comm);
    return ret;
}


extern "C" {
void
F77_MPI_COMM_FREE(MPI_Fint *comm, MPI_Fint *ierr)
{
    int ret;
    MPI_Comm c_comm;

    c_comm = MPI_Comm_f2c(*comm);
    ret = MPI_Comm_free(&c_comm);
    *ierr = ret;
    return;
}
}

static int
_Finalize(void) {
    prof_attrs *array = NULL;
    int rank, size;
    int i, k, j, len;
    prof_attrs *recv_buffer = NULL;
    // char **names, **unames;
    int  num_of_comms, resultlen;
    // uint64_t *bytes, *prims_bytes;
    // uint32_t *prims;
    // uint64_t *msgs;
    // double *time_info;
    // int *sizes;
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    char proc_name[MPI_MAX_PROCESSOR_NAME];
    char *proc_names = NULL;
    char *ptr;
    double *alltimes = NULL;
    std::vector<double> mpi_times;
    // double mpi_time = 0.0;
    total_time = MPI_Wtime() - total_time;

    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);
    num_of_comms = local_communicators.size();


    array = (prof_attrs *) malloc(sizeof(prof_attrs) * num_of_comms);
    if (array == NULL) {
        mcpt_abort("malloc error for send buffer Rank: %d\n", rank);
    }
    if (rank == 0) {
        recv_buffer =
                (prof_attrs *) malloc(sizeof(prof_attrs) * num_of_comms * size);
        if (recv_buffer == NULL) {
            mcpt_abort("malloc error for receive buffer Rank: %d\n", rank);
        }
    }


    MPI_Datatype profiler_data;
    MPI_Aint base, displacements[4];
    int blocklengths[4] = {NAMELEN, 1, NUM_OF_PRIMS * NUM_BUCKETS, NUM_OF_PRIMS * NUM_BUCKETS};
    MPI_Datatype types[4] = {MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_UINT64_T};
    prof_attrs dummy;

    // Create a dummy instance to calculate displacements
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.name, &displacements[0]);
    MPI_Get_address(&dummy.size, &displacements[1]);
    MPI_Get_address(&dummy.buckets_time, &displacements[2]);
    MPI_Get_address(&dummy.buckets_msgs, &displacements[3]);

    // Convert addresses to displacements
    for (int i = 0; i < 4; i++) {
        displacements[i] = MPI_Aint_diff(displacements[i], base);
    }

    MPI_Type_create_struct(4, blocklengths, displacements, types, &profiler_data);
    PMPI_Type_commit(&profiler_data);
    k = 0;

    for (i = 0; i < num_of_comms; i++) {
         memcpy(&array[i], local_communicators[i], sizeof(prof_attrs));
    }

    MPI_Get_processor_name(proc_name, &len);

    if (rank == 0) {
        proc_names = NULL;
        proc_names = (char *)malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME * size);
        if (proc_names == NULL) {
            mcpt_abort("malloc error for proc_names buffer Rank: %d\n",rank);
        }
        alltimes = NULL;
        alltimes = (double*) malloc (sizeof(double)*size);
        if (alltimes == NULL) {
            mcpt_abort("malloc error for alltimes buffer Rank: %d\n",rank);
        }
        // iterate over local_communicators to see if they have something meaningful
     }

    PMPI_Gather(proc_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, proc_names,MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0 , MPI_COMM_WORLD);

    PMPI_Gather(&total_time, 1, MPI_DOUBLE, alltimes,1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);

    PMPI_Gather(array, num_of_comms, profiler_data, recv_buffer,
                num_of_comms, profiler_data, 0, MPI_COMM_WORLD);

    PMPI_Barrier(MPI_COMM_WORLD);
    if ( rank == 0 ){
        int rc,commId,maxsize,minsize;
        int powers_of_2[NUM_BUCKETS - 1];
        sqlite3 *db = NULL;
        const char *env_var = getenv("MPISEE_OUTFILE");
        char *outfile;
        int r = -1;
        if (env_var != NULL) {
          rc = sqlite3_open(env_var, &db);
          if (rc) {
              std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
              return 1;
          } else {
              std::cout << "Opened database successfully" << std::endl;
          }
          outfile = strdup(env_var);
        }
        else{
          rc = sqlite3_open("mpisee_profile.db", &db);
          if (rc) {
              std::cerr << "Can't open database: " << sqlite3_errmsg(db) << std::endl;
              return 1;
          } else {
              std::cout << "Opened database successfully" << std::endl;
          }
          outfile = strdup("mpisee_profile.db");
        }
        PMPI_Get_library_version(version, &resultlen);
        /* Remove line breaks in MPI version string, as it may create bugs during parsing. */
        for (i = 0; i < strlen(version); i++)
        {
            if(version[i] == '\n')
            {
                version[i] = ' ';
            }
        }
        createTables(db);
        insertMetadata(db, version, size, av, ac, mpisee_major_version,
                       mpisee_minor_version, mpisee_build_date,
                       mpisee_build_time, env_var);

        for (i = 0; i<NUM_OF_PRIMS ; i++ ) {
            insertIntoOperations(db, prim_names[i]);
        }


        ptr = proc_names;
        for (i = 0; i < size; i++) {
            if ( ptr != NULL ){
              snprintf(proc_name, MPI_MAX_PROCESSOR_NAME, "%s", ptr);
            }
            insertIntoMappings(db, proc_name);
            ptr+=MPI_MAX_PROCESSOR_NAME;
        }
        // Precompute powers of 2 for each bucket
        for (int i = 0; i < NUM_BUCKETS - 1; i++) {
            powers_of_2[i] = 1 << buckets[i];
        }

        for (i = 0; i < num_of_comms*size; i++) {
            insertIntoComms(db, recv_buffer[i].name, recv_buffer[i].size);
            commId = getCommId(db, recv_buffer[i].name);
            if (i % num_of_comms == 0) {
                r++;
            }
            for (k = 0; k < NUM_OF_PRIMS; k++) {
                // Handle the first bucket separately
                minsize = 0;
                maxsize = powers_of_2[0];
                if (recv_buffer[i].buckets_msgs[k][0] > 0) {
                    insertIntoData(db, r, commId, k, maxsize, minsize,
                                   recv_buffer[i].buckets_msgs[k][0],
                                   recv_buffer[i].buckets_time[k][0]);
                }
                for (j = 1; j < NUM_BUCKETS-1; j++) {
                    minsize = powers_of_2[j - 1];
                    maxsize = (j == NUM_BUCKETS - 2) ? INT_MAX : powers_of_2[j];
                    if (recv_buffer[i].buckets_msgs[k][j] > 0) {
                        insertIntoData(db, r, commId, k, maxsize, minsize,
                                       recv_buffer[i].buckets_msgs[k][j],
                                       recv_buffer[i].buckets_time[k][j]);
                    }

                }
            }
        }

        printf("Database File Written: %s\n", outfile);

        printMetadata(db);
        //printCommsTable(db);
        //printData(db);
        sqlite3_close(db);

        free(outfile);
        free(alltimes);
        free(proc_names);
    }

    PMPI_Type_free(&profiler_data);
    MPI_Barrier(MPI_COMM_WORLD);
    free(array);

    return PMPI_Finalize();
}


int
MPI_Finalize (void)
{
  int rc = 0;

  rc = _Finalize ();

  return rc;
}

extern "C" {
void
F77_MPI_FINALIZE (int *ierr)
{
  int rc = 0;

  rc = _Finalize ();
  *ierr = rc;

  return;
}
}

/* void apply_print(const void *key,void **value,void *cl){ */
/*     prof_attrs *val; */
/*     unsigned long long bytes; */
/*     char name[32], prim[32]; */
/*     val = (prof_attrs*)value; */
/*     strcpy(name, val->name); */
/*     bytes = val->bytes; */
/*     strcpy(prim, val->prim); */
/*     printf("Name: %s created by %s bytes transferred %llu\n",name,prim,bytes); */
/* } */

/* int compare(const void *x, const void *y) { */
/*     prof_attrs *a = (prof_attrs*) x; */
/*     prof_attrs *b = (prof_attrs*) y; */
/*     char bufa[32], bufb[32]; */
/*     strcpy(bufa,a->name); */
/*     strcpy(bufb,b->name); */
/*     printf("Compare %s with %s\n",bufa,bufb); */
/*     return strcmp(a->name, b->name); */
/* } */

/* void print_entry(commtor* communicator){ */
/*     printf("Name: %s, Bytes  %ld, Prim created %s\n",communicator->name, */
/*            communicator->bytes,communicator->prim); */
/* } */

/* static void __attribute__((no_instrument_function)) */
/* log_func(const void* funcAddr, const char* action, const void* callSite){ */
/*     fprintf(stdout,"%p %d \n",callSite,__LINE__); */
/* } */

/* void __attribute__((no_instrument_function)) */
/* __cyg_profile_func_enter(void* this_fn, void* call_site){ */
/*     log_func(this_fn, "->", call_site); */
/* } */

/* void __attribute__((no_instrument_function)) */
/* __cyg_profile_func_exit(void* this_fn, void* call_site){ */
/*     log_func(this_fn, "<-", call_site); */
/* } */

/* static int */
/* compare_int(const void *a, const void *b){ */
/*     int *a0, *b0; */
/*     a0 = (int*)a; */
/*     b0 = (int*)b; */
/*     return (*a0)-(*b0); */
/* } */
