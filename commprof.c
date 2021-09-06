#include "utils.h"
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <unistd.h>
#define MAX_ARGS 32
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


MPI_Comm *communicators = NULL;
prof_attrs **local_data = NULL;
int num_of_comms = 1;
int my_coms = 1;
int num_of_local = 0;
int ac;
char *av[MAX_ARGS];
/* int line_called; */
/* char file_called[32]; */

static int
namedel(MPI_Comm comm, int keyval, void *attr, void *s)
{
  prof_attrs *com = (prof_attrs*)attr;
  if ( comm == NULL )
      return -1;
  else
      free(com);
  return 0;
}

static int
namekey()
{
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}


prof_attrs*
get_comm_parent(MPI_Comm comm)
{
    int i,flag;
    prof_attrs *communicator, *com_info;
    /* communicator = (prof_attrs*) malloc (sizeof(prof_attrs)); */
    communicator = (prof_attrs*) calloc(1,sizeof(prof_attrs));
    if ( comm != MPI_COMM_WORLD ){
        //find the mother
        for ( i = 0; i<my_coms; i++ ){
            if ( comm  == communicators[i] )
                break;
        }
        if ( i == my_coms ){
            // Should have found the mother
            fprintf(stderr, "Error: could not find the mother of communicator.\n");
            mcpt_abort("File:%s line:%d Aborting\n",__FILE__,__LINE__);

        }
        else{
            PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag);
            if ( flag ){
                strcpy(communicator->name, com_info->name);
            }
            else{
                mcpt_abort("Flag in file:%s line:%d invalid\nAborting\n",__FILE__,__LINE__);
            }
        }
    }
    else{
        strcpy(communicator->name, "W");
    }
    return communicator;
}



static int
_MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
    int ret,rank,size;
    int i;
    prof_attrs *communicator;
    ret = PMPI_Init_thread(argc, argv, required, provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    local_data = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4);
    for ( i =0 ; i<size*4; i++ ){
        communicators[i] = NULL;
        local_data[i] = NULL;
    }
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*256);
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    strcpy(communicator->name,"W");
    communicator->bytes = 0;
    communicator->size = size;
    communicator->msgs = 0;
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[0] = MPI_COMM_WORLD;
    return ret;
}

static int
_MPI_Init(int *argc, char ***argv){
    int ret,rank,size;
    int i;
    prof_attrs *communicator;
    ret = PMPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    local_data = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4);
    for ( i =0 ; i<size*4; i++ ){
        communicators[i] = NULL;
        local_data[i] = NULL;
    }
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*256);
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    strcpy(communicator->name,"W");
    communicator->bytes = 0;
    communicator->size = size;
    communicator->msgs = 0;
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[0] = MPI_COMM_WORLD;
    return ret;
}

extern int
MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
    return _MPI_Init_thread(argc, argv, required, provided);
}

extern void
F77_MPI_INIT_THREAD (int *required, int *provided, int *ierr){
    char **tmp;
    int ret;
    getProcCmdLine (&ac, av);
    tmp = av;
    printf("F77_INIT_THREAD\n");
    ret = _MPI_Init_thread(&ac, (char***)&tmp , *required, provided);
    *ierr = ret;
}

extern void
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

extern int
MPI_Init(int *argc, char ***argv)
{
    int ret,rank,size;
    int i;
    prof_attrs *communicator;
    /* unsigned long long bytes; */
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    local_data = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4);
    for ( i =0 ; i<size*4; i++ ){
        communicators[i] = NULL;
        local_data[i] = NULL;
    }
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*256);
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    strcpy(communicator->name,"W");
    communicator->bytes = 0;
    communicator->size = size;
    communicator->msgs = 0;
    /* communicator->index = num_of_comms; */
    /* communicator->size = size; */
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[0] = MPI_COMM_WORLD;

    /* Table_put(table, MPI_COMM_WORLD, communicator); */
    /* num_of_comms++; */
    /* my_coms++; */
    return ret;
}

extern int
MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    int ret;
    int i,length,comms;
    char *buf;
    prof_attrs *communicator;
    ret = PMPI_Comm_create(comm, group, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        communicators[num_of_comms] = NULL;
        num_of_comms++;
        return ret;
    }
    MPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    communicator = get_comm_parent(comm);
    buf = (char*) malloc ( sizeof(char)*8);
    sprintf(buf,"_c%d",my_coms);
    length = strlen(communicator->name);
    for ( i =0; i<strlen(buf); i++ ){
        communicator->name[length+i] = buf[i];
    }
    /* communicator->name[i]='\0'; */
    communicator->bytes = 0;
    communicator->msgs = 0;
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[my_coms] = *newcomm;
    num_of_comms++;
    my_coms++;
    return ret;
}

extern int
MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    int ret,i,length,comms;
    prof_attrs *communicator;
    char *buf;
    ret = PMPI_Comm_split(comm, color, key, newcomm);
    MPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    communicator = get_comm_parent(comm);
    buf = (char*) malloc ( sizeof(char)*16);
    sprintf(buf,"_s%d.%d",my_coms,color);
    length = strlen(communicator->name);
    for ( i =0; i<strlen(buf); i++ ){
        communicator->name[length+i] = buf[i];
    }
    /* communicator->name[i+1]='\0'; */
    communicator->bytes = 0;
    communicator->msgs = 0;
    /* printf("MPI_Comm_split comm with name %s and %c\n",communicator->name,communicator->name[i]); */
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[my_coms] = *newcomm;
    num_of_comms+=color+1;
    my_coms++;
    return ret;
}

extern int
MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
    int ret,length,i, comms;
    prof_attrs *communicator;
    char *buf;
    ret = PMPI_Comm_dup(comm, newcomm);
    /* if ( comm == MPI_COMM_WORLD ){ */
        /* Synchronize */
    MPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
    my_coms = comms;
    /* } */
    communicator = get_comm_parent(comm);
    buf = (char*) malloc ( sizeof(char)*8);
    sprintf(buf,"_d%d",my_coms);
    length = strlen(communicator->name);
    for ( i =0; i<strlen(buf); i++ ){
        communicator->name[length+i] = buf[i];
    }
    /* communicator->name[i]='\0'; */
    communicator->bytes = 0;
    communicator->msgs = 0;
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[my_coms] = *newcomm;
    num_of_comms+=1;
    my_coms++;
    return ret;
}

extern int
MPI_Cart_create(MPI_Comm old_comm, int ndims, const int *dims,
                const int *periods, int reorder, MPI_Comm *comm_cart)
{

    int ret,length,i,comms;
    prof_attrs *communicator;
    char *buf, *buffer;
    ret = PMPI_Cart_create(old_comm, ndims, dims, periods, reorder, comm_cart);
    MPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, old_comm);
    my_coms = comms;
    communicator = get_comm_parent(old_comm);
    buffer = strdup(communicator->name);
    /* buf = (char*) malloc ( sizeof(char)*8); */
    buf = (char*)calloc(8, sizeof(char));
    /* snprintf(buf, 8, "_a%d",my_coms); */
    sprintf(buf,"_a%d",my_coms);
    length = strlen(communicator->name);
    for ( i =0; i<strlen(buf); i++ ){
        communicator->name[length+i] = buf[i];
    }
    /* communicator->name[i+1]='\0'; */
    communicator->bytes = 0;
    communicator->msgs = 0;
    PMPI_Comm_set_attr(*comm_cart, namekey(), communicator);
    /* printf("MPI Cart comm with name %s and parent %s\n",communicator->name,buffer); */
    communicators[my_coms] = *comm_cart;
    num_of_comms+=1;
    my_coms++;
    free(buffer);
    return ret;
}

extern int
MPI_Isend(const void *buf, int count, MPI_Datatype datatype,int dest, int tag,
          MPI_Comm comm, MPI_Request *request)
{
    int ret;
    int size,flag,i;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }

    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(datatype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + count * size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
    }
    return ret;
}

extern int
MPI_Send(const void *buf, int count, MPI_Datatype datatype,
         int dest,int tag, MPI_Comm comm)
{
    int ret;
    int size,flag,i;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }

    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(datatype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + count * size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
    }
    return ret;
}

extern int
MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             int dest, int sendtag, void *recvbuf, int recvcount,
             MPI_Datatype recvtype, int source, int recvtag,
             MPI_Comm comm, MPI_Status *status)
{

    int ret;
    int size,flag,i;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf,
                        recvcount, recvtype, source, recvtag, comm, status);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }

    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(sendtype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sendcount * size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
    }
    return ret;
}

extern int
MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root,
          MPI_Comm comm)
{
    int ret,i,flag,size,comm_size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Bcast(buffer, count, datatype, root, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(datatype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = (sum + count * size);
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
    }
    return ret;

}

extern int
MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{

    int ret,i,flag,size,comm_size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(datatype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + (count * size)*comm_size + (count * size);
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
    }
    return ret;
}

extern int
MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
              void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret,i,flag,size,comm_size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                         recvtype, comm);

    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(sendtype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + (sendcount * size)*comm_size + (recvcount*size);
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
    }
    return ret;
}

extern int
MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
             void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
{
    int ret,i,flag,size,comm_size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                         recvtype, comm);

    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(sendtype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sendcount*size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
    }
    return ret;
}

extern int
MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
              const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
              const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm)
{

    int ret,i,j,flag,size,comm_size,sent,recvd;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts,
                   rdispls, recvtype, comm);

    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(sendtype, &size);
    PMPI_Comm_size(comm, &comm_size);
    sent = 0;
    for ( i = 0; i<comm_size; i++ ){
        sent += sendcounts[i];
    }
    j = sent;
    sent = sent*size;
    recvd = 0;
    for ( i = 0; i<comm_size; i++ ){
        recvd += recvcounts[i];
    }
    recvd = recvd*size;
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sent;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + j;
    }
    return ret;
}

extern int
MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, const int *recvcounts, const int *displs,
               MPI_Datatype recvtype, MPI_Comm comm)
{

    int ret,i,j,flag,size,comm_size,recvd;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                          displs, recvtype, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(sendtype, &size);
    PMPI_Comm_size(comm, &comm_size);
    recvd = 0;
    j = 0;
    for ( i = 0; i<comm_size; i++ ){
        recvd += recvcounts[i];
        j++;
    }
    recvd = recvd*size;
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + recvd;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + j;
    }
    return ret;
}


extern int
MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
           MPI_Op op, int root, MPI_Comm comm)
{
    int ret,i,flag,size,comm_size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;

    ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    PMPI_Type_size(datatype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + count*size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
    }
    return ret;
}

extern int
MPI_Comm_free(MPI_Comm *comm){
    int ret,flag,i;
    prof_attrs *com_info;
    PMPI_Comm_get_attr(*comm, namekey(), &com_info, &flag);
    if ( flag ){
        local_data[num_of_local] = (prof_attrs*) malloc (sizeof(prof_attrs));
        local_data[num_of_local]->bytes = com_info->bytes;
        local_data[num_of_local]->msgs = com_info->msgs;
        local_data[num_of_local]->size = com_info->size;
        strcpy(local_data[num_of_local]->name,com_info->name);
        num_of_local++;
    }
    for ( i = 0; i<num_of_comms; i++ ){
        if ( *comm  == communicators[i])
            break;
    }
    communicators[i]=NULL;
    ret = PMPI_Comm_free(comm);
    return ret;
}

static int
_Finalize(){
    FILE *fp;
    prof_attrs *array;
    int rank,size;
    int i,j,k,flag,found;
    prof_attrs *com_info;
    prof_attrs *recv_buffer;
    prof_attrs dummy;
    char **names, **unames;
    int total_comms,total;
    uint64_t *bytes, *ubytes;
    uint32_t *msgs, *umsgs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    PMPI_Allreduce(&num_of_comms, &total_comms, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    num_of_comms = total_comms;


    array =(prof_attrs*) malloc(sizeof(prof_attrs)*num_of_comms);
    recv_buffer = (prof_attrs*) malloc (sizeof(prof_attrs)*num_of_comms*size);
    MPI_Datatype types[4] = { MPI_CHAR,MPI_UINT64_T, MPI_UINT32_T, MPI_INT };
    int blocklen[4] = {NAMELEN,1,1,1};
    MPI_Aint displacements[4];
    MPI_Aint base_address;
    MPI_Datatype profiler_data;
    MPI_Get_address(&dummy, &base_address);
    MPI_Get_address(&dummy.name[0], &displacements[0]);
    /* MPI_Get_address(&dummy.parent[0], &displacements[1]); */
    /* MPI_Get_address(&dummy.prim[0], &displacements[2]); */
    MPI_Get_address(&dummy.bytes, &displacements[1]);
    MPI_Get_address(&dummy.msgs, &displacements[2]);
    MPI_Get_address(&dummy.size, &displacements[3]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    displacements[3] = MPI_Aint_diff(displacements[3], base_address);
    /* displacements[4] = MPI_Aint_diff(displacements[4], base_address); */

    MPI_Type_create_struct(4, blocklen, displacements, types, &profiler_data);
    MPI_Type_commit(&profiler_data);
    k = 0;
    for ( i = 0; i < num_of_comms; i++ ){
        if ( communicators[i] != NULL ){
            PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag);
            if ( flag ){
                strcpy(array[i].name, com_info->name);
                /* strcpy(array[i].parent, com_info->parent); */
                /* strcpy(array[i].prim, com_info->prim); */
                array[i].bytes = com_info->bytes;
                array[i].msgs= com_info->msgs;
                array[i].size = com_info->size;
            }
        }
        else{
            if ( num_of_local > 0 && k < num_of_local){
                strcpy(array[i].name, local_data[k]->name);
                /* strcpy(array[i].parent, local_data[k]->parent); */
                /* strcpy(array[i].prim, local_data[k]->prim); */
                array[i].bytes = local_data[k]->bytes;
                array[i].size = local_data[k]->size;
                array[i].msgs = local_data[k]->msgs;
                k++;
            }
            else{
                strcpy(array[i].name, "NULL");
                /* strcpy(array[i].parent, "NULL"); */
                /* strcpy(array[i].prim, "NULL"); */
                array[i].msgs = 0;
                array[i].bytes = 0;
                array[i].size = 0;
            }
        }
    }
    fp = fopen("profiler_stats.txt","w");
    PMPI_Gather(array, num_of_comms*sizeof(prof_attrs), MPI_BYTE, recv_buffer,
                num_of_comms*sizeof(prof_attrs), MPI_BYTE, 0, MPI_COMM_WORLD);

    if ( rank == 0 ){
        names = ( char**)malloc(sizeof(char*)*num_of_comms*size);
        /* parents = ( char**)malloc(sizeof(char*)*num_of_comms*size); */
        unames = (char **) malloc (sizeof(char*)*num_of_comms*size);
        /* uparents =(char **) malloc (sizeof(char*)*num_of_comms*size); */
        bytes = (uint64_t *) malloc (sizeof(uint64_t )
                                               *num_of_comms*size);
        msgs = (uint32_t *) malloc (sizeof(uint32_t )
                                               *num_of_comms*size);
        j = 0;
        for ( i =0; i<size*num_of_comms; i++ ){
            if ( strcmp(recv_buffer[i].name, "NULL") != 0 ){
                names[j]=strdup(recv_buffer[i].name);
                /* parents[j]=strdup(recv_buffer[i].parent); */
                unames[j] = strdup("NULL");
                /* uparents[j] = strdup("NULL"); */
                bytes[j] = recv_buffer[i].bytes;
                msgs[j] = recv_buffer[i].msgs;
                j++;
            }
        }
        total = j;
        ubytes = (uint64_t *) malloc (sizeof(uint64_t )
                                      *total);
        umsgs = (uint32_t *) malloc (sizeof(uint32_t )
                                                *total);
        memset(ubytes, 0, sizeof(uint64_t )*total);
        memset(umsgs, 0, sizeof(uint32_t )*total);
        num_of_comms = 1;
        j = 0;
        for ( i=0; i<total; i++ ){
            /* Build the global communicator tree */
            if ( strcmp(names[i], "W") != 0 ){
                found = 0;
                for ( k =0; k<total; k++ ){
                    if ( strcmp(names[i], unames[k] ) == 0 ){
                         /* strcmp(parents[i], uparents[k] ) == 0 ) */
                        found = 1;
                    }
                }
                if ( !found ){
                    strcpy(unames[j], names[i]);
                    /* strcpy(uparents[j], parents[i]); */
                    j++;
                    num_of_comms++;
                }
            }
        }
        strcpy(unames[j], "W");
        /* strcpy(uparents[j],"NULL"); */
        for ( i = 0; i<num_of_comms; i++){
            /* printf("Searching for %s parent %s\n",unames[i],uparents[i]); */
            for ( j=0; j<total; j++ ){
                /* printf("Found Comm: %s parent: %s Bytes= %llu\n",names[j], */
                /*        parents[j],bytes[j]); */
                if ( strcmp(unames[i], names[j]) == 0 ){
                     /* strcmp(uparents[i], parents[j]) == 0){ */
                    ubytes[i]+= bytes[j];
                    umsgs[i]+= msgs[j];

                    /* printf("Match!\n"); */
                }
            }
        }

        for ( i =0; i<total; i++ ){
            free(names[i]);
            /* free(parents[i]); */
        }
        free(names);
        /* free(parents); */
        free(bytes);
        free(msgs);

        if (fp == NULL){
            fprintf(stderr, "Failed to open output file: profiler_stats.txt\n");
            mcpt_abort("Aborting\n");
        }
        fprintf(fp,"Num of REAL comms = %d\n",num_of_comms);
        for ( i =0; i<num_of_comms; i++ )
            fprintf(fp,"Comm: %s Bytes = %lu Msgs = %u\n",unames[i],ubytes[i],
                   umsgs[i]);

        for ( i =0; i<total; i++ ){
            free(unames[i]);
            /* free(uparents[i]); */
        }
        free(unames);
        /* free(uparents); */
        free(ubytes);
        free(umsgs);
    }
    fclose(fp);
    /* if (communicators) */
    /*     FREE(communicators); */
    return PMPI_Finalize();
}


extern int
MPI_Finalize (void)
{
  int rc = 0;

  rc = _Finalize ();

  return rc;
}

extern void
F77_MPI_FINALIZE (int *ierr)
{
  int rc = 0;

  rc = _Finalize ();
  *ierr = rc;

  return;
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
