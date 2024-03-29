#include "utils.h"
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <mpi.h>
#include <ostream>
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
int mpisee_major_version = 1;
int mpisee_minor_version = 0;
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

extern "C" {
int
win_namekey(void){
    static int win_keyval = MPI_KEYVAL_INVALID;

    if (win_keyval == MPI_KEYVAL_INVALID) {
        MPI_Win_create_keyval(MPI_WIN_NULL_COPY_FN, MPI_WIN_NULL_DELETE_FN, &win_keyval, NULL);
    }
    return win_keyval;

}
}


/*
 * Create a new communicator structure and add the parent
 * communicator's name prefix to it
 */
extern "C" {
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

int
choose_bucket(int64_t bytes) {
    int index;
    int64_t tmp;
    for (index = 0; index < NUM_BUCKETS-1; index++) {
        tmp = buckets[index];
        if (tmp > bytes) {
            break;
        }
    }
    return index;
}


extern "C" {
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

std::vector<std::string> convertToArrayOfStrings(char *proc_names, int size,
                                                 int name_length) {
    std::vector<std::string> machineNames;
    for (int i = 1; i < size; ++i) {
        std::string machineName(proc_names + i * name_length);
        machineNames.push_back(machineName);
    }
    return machineNames;
}

std::vector<std::string> convertToArrayOfPrims() {

    std::vector<std::string> machineNames;
    for (int i = 1; i < NUM_OF_PRIMS; ++i) {
        std::string PrimName(prim_names[i]);
        machineNames.push_back(PrimName);
    }
    return machineNames;
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
MPI_Wait(MPI_Request *request, MPI_Status *status)
{
    int ret;
    double t_elapsed;

    MPI_Comm comm ;
    if ( prof_enabled == 1 ){
        comm = requests_map[*request];
        t_elapsed = MPI_Wtime();
        ret = PMPI_Wait(request, status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm == NULL  ){
            fprintf(stderr, "mpisee: NULL COMMUNICATOR in MPI_Wait\n");
            return ret;
        }
       profile_this(comm, 0, MPI_DATATYPE_NULL, Wait, t_elapsed, 0);
       requests_map.erase(*request);
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
    int ret,i,j;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    j = 0;
    if ( prof_enabled == 1 ){
        for ( i =0; i<count; i++ ){
            if ( j == 0 )
                comm = requests_map[array_of_requests[i]];
            j++;
            requests_map.erase(array_of_requests[i]);
        }
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm != MPI_COMM_NULL)
            profile_this(comm, 0, MPI_DATATYPE_NULL, Waitall, t_elapsed, 0);
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

    int ret,i,j;
    double t_elapsed;
    MPI_Comm comm = MPI_COMM_NULL;
    j = 0;

    if ( prof_enabled == 1 ){
        for ( i =0; i<count; i++ ){
            if ( j == 0 )
                comm = requests_map[array_of_requests[i]];
            j++;
            requests_map.erase(array_of_requests[i]);
        }
        t_elapsed = MPI_Wtime();
        ret = PMPI_Waitany(count, array_of_requests,index,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm != MPI_COMM_NULL)
            profile_this(comm, 0, MPI_DATATYPE_NULL, Waitany, t_elapsed, 0);
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

    MPI_Comm comm = NULL;
    if ( prof_enabled == 1 ){
        comm = requests_map[*request];
        t_elapsed = MPI_Wtime();
        ret = PMPI_Test(request,flag,status);
        t_elapsed = MPI_Wtime() - t_elapsed;
        if ( comm == NULL ){
            return ret;
        }
        profile_this(comm, 0, MPI_DATATYPE_NULL, Test, t_elapsed, 0);
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
    MPI_Comm comm = NULL;

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
    int  num_of_comms, resultlen;
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    char proc_name[MPI_MAX_PROCESSOR_NAME];
    char *proc_names = NULL;
    double *alltimes = NULL;
    std::vector<double> mpi_times;
    int *recvcounts = NULL;
    int *displs = NULL;
    int total_num_of_comms;
    total_time = MPI_Wtime() - total_time;

    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);
    num_of_comms = local_communicators.size();

    // Do all processes have the same number of communicators?
    recvcounts = (int *)malloc(sizeof(int) * size);
    if (recvcounts == NULL) {
        mcpt_abort("malloc error for recvcounts Rank: %d\n", rank);
    }
    if ( rank == 0 ){
        displs = (int *)malloc(sizeof(int) * size);
        if (displs == NULL) {
            mcpt_abort("malloc error for displs");
        }
    }
    PMPI_Gather(&num_of_comms, 1, MPI_INT, recvcounts, 1, MPI_INT, 0,
                MPI_COMM_WORLD);

    if (rank == 0) {
        // Compute displacements
        displs[0] = 0;
        total_num_of_comms = recvcounts[0];
        for (i = 1; i < size; ++i) {
          displs[i] = displs[i - 1] + recvcounts[i - 1];
          total_num_of_comms += recvcounts[i];
        }
        std::cout << "mpisee: total number of communicators = " << total_num_of_comms << std::endl;
        recv_buffer =
            (prof_attrs *)malloc(sizeof(prof_attrs) * total_num_of_comms );

        if (recv_buffer == NULL) {
          mcpt_abort("malloc error for receive buffer Rank: %d\n", rank);

        }
    }

    array = (prof_attrs *)malloc(sizeof(prof_attrs) * num_of_comms);
    if (array == NULL) {
        mcpt_abort("malloc error for send buffer Rank: %d\n", rank);
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
    for (i = 0; i < 4; i++) {
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
        proc_names = (char *)malloc(sizeof(char) * MPI_MAX_PROCESSOR_NAME * size);
        if (proc_names == NULL) {
          mcpt_abort("malloc error for proc_names buffer Rank: %d\n", rank);
        }
        alltimes = (double *)malloc(sizeof(double) * size);
        if (alltimes == NULL) {
          mcpt_abort("malloc error for alltimes buffer Rank: %d\n", rank);
        }
     }

    PMPI_Gather(proc_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, proc_names,
                 MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);


    PMPI_Gather(&total_time, 1, MPI_DOUBLE, alltimes, 1, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);


    PMPI_Gatherv(array, num_of_comms, profiler_data, recv_buffer, recvcounts,
                 displs, profiler_data, 0, MPI_COMM_WORLD);


    if ( rank == 0 ){
        int rc,commId,maxsize,minsize;
        //int powers_of_2[NUM_BUCKETS - 1];
        sqlite3 *db = NULL;
        char *outfile = NULL;
        int l, proc, startIdx, numElements;
        double t;
        const char *env_var = getenv("MPISEE_OUTFILE");
        if (env_var != NULL) {
          rc = sqlite3_open(env_var, &db);
          if (rc) {
              std::cerr << "mpisee: Can't open database: " << sqlite3_errmsg(db) << std::endl;
              return 1;
          } else {
              std::cout << "mpisee: Opened database successfully" << std::endl;
          }
          outfile = strdup(env_var);
          if (outfile == NULL) {
              mcpt_abort("mpisee: strdup returned NULL\n");
          }
        }
        else{
          rc = sqlite3_open("mpisee_profile.db", &db);
          if (rc) {
              std::cerr << "mpisee: Can't open database: " << sqlite3_errmsg(db) << std::endl;
              return 1;
          } else {
              std::cout << "mpisee: Opened database successfully" << std::endl;
          }
          outfile = strdup("mpisee_profile.db");
          if (outfile == NULL) {
              mcpt_abort("mpisee: strdup returned NULL\n");
          }
        }

        PMPI_Get_library_version(version, &resultlen);
        for (i = 0; i < strlen(version); i++)
        {
            if(version[i] == '\n')
            {
                version[i] = ' ';
            }
        }
        createTables(db);
        std::cout << "mpisee: Writing the metadata table" << std::endl;

        insertMetadata(db, version, size, av, ac, mpisee_major_version,
                       mpisee_minor_version, mpisee_build_date,
                       mpisee_build_time, env_var);

        std::cout << "mpisee: Writing the MPI operations table" << std::endl;


        insertIntoOperationsEmpty(db, prim_names[0]);
        std::vector<std::string> operations = convertToArrayOfPrims();
        BatchInsertIntoOperations(db, operations);
        operations.clear();
        operations.shrink_to_fit();

        std::vector<double> times;
        if (alltimes != NULL){
          std::cout << "mpisee: Writing the exectimes table" << std::endl;
          insertIntoTimes(db, alltimes[0]);
          for (i = 1; i < size; i++) {
            times.push_back(alltimes[i]);
          }
          BatchInsertIntoTimes(db, times);
          times.clear();
          times.shrink_to_fit();
          free(alltimes);
        } else {
          std::cout << "mpisee: Execution times NULL" << std::endl;
        }

        std::string machineName(proc_names);
        insertIntoMappings(db, machineName);
        std::vector<std::string> machines =
            convertToArrayOfStrings(proc_names, size, MPI_MAX_PROCESSOR_NAME);
        BatchInsertIntoMappings(db, machines);
        machines.clear();
        machines.shrink_to_fit();
        free(proc_names);

        // Precompute powers of 2 for each bucket
        // for (i = 0; i < NUM_BUCKETS - 1; i++) {
        //     powers_of_2[i] = 1 << buckets[i];
        // }

        std::vector<CommData> comms;
        std::vector<int> commIds;

        for (proc = 0; proc < size; ++proc) {
            startIdx = displs[proc];
            numElements = recvcounts[proc];

          for (j = 0; j < numElements; ++j) {
              comms.push_back({recv_buffer[startIdx + j].name,
                      recv_buffer[startIdx + j].size});
          }
        }

        commIds=CommsInsert(db, comms);
        comms.clear();
        comms.shrink_to_fit();

        std::vector<DataEntry> entries;
        std::cout << "mpisee: Writing the main data table"
                  << std::endl;
        i = 0;
        t = MPI_Wtime();
        commId = 0;
        for (proc = 0; proc < size; ++proc) {

          startIdx = displs[proc];
          numElements = recvcounts[proc];

          for (j = 0; j < numElements; ++j) {
            if (i < commIds.size()) {
                commId = commIds[i];
            } else {
              std::cout << "mpisee: index in commids (" << i
                        << ") out of bounds" << std::endl;
              mcpt_abort("commId out of bounds\n");

            }
            prof_attrs &item = recv_buffer[startIdx + j];

            for (k = 0; k < NUM_OF_PRIMS; k++) {
              minsize = 0;
              maxsize = buckets[0]; //powers_of_2[0];
              if (item.buckets_msgs[k][0] > 0) {
                  insertIntoDataEntry(entries, proc, commId, k, maxsize, minsize,
                                      item.buckets_msgs[k][0],
                                      item.buckets_time[k][0]);
              }
              for (l = 1; l < NUM_BUCKETS-1; l++) {
                  minsize = buckets[l-1]; //powers_of_2[l - 1];
                  maxsize = (l == NUM_BUCKETS - 2) ? INT_MAX : buckets[l]; //powers_of_2[l];
                  if (item.buckets_msgs[k][l] > 0) {
                    insertIntoDataEntry(entries, proc, commId, k, maxsize,
                                        minsize, item.buckets_msgs[k][l],
                                        item.buckets_time[k][l]);
                  }
              }
            }
            i++;
          }
        }

        executeBatchInsert(db, entries);
        t = MPI_Wtime() - t;
        std::cout << "mpisee: Output database file: " << outfile << ", time to write: " << t << " seconds" << std::endl;
        sqlite3_close(db);
        free(outfile);
        free(recv_buffer);

    }


    MPI_Barrier(MPI_COMM_WORLD);
    PMPI_Type_free(&profiler_data);
    free(array);
    free(displs);
    free(recvcounts);

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

