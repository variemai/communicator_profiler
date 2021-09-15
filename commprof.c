#include "utils.h"
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <unistd.h>
#include "symbols.h"
#define MAX_ARGS 1024
#define MAX_DIMS 8

MPI_Comm *communicators = NULL;
prof_attrs **local_data = NULL;
prof_attrs **local_comms = NULL;
int local_cid= 0;
/* int num_of_comms = 1; */
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
  /* printf("Namedel %s\n",com->name); */
  /* fflush(stdout); */
  free(com);
  return MPI_SUCCESS;
}

static int
namekey()
{
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    PMPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}

/* static int */
/* compare_int(const void *a, const void *b){ */
/*     int *a0, *b0; */
/*     a0 = (int*)a; */
/*     b0 = (int*)b; */
/*     return (*a0)-(*b0); */
/* } */

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
            fprintf(stderr, "Error: could not find the parent of communicator.\n");
            mcpt_abort("File:%s line:%d Aborting\n",__FILE__,__LINE__);

        }
        else{
            PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag);
            if ( flag ){
                strcpy(communicator->name, com_info->name);
                /* printf("Parent = %s\n",communicator->name); */
                /* fflush(stdout); */
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

void
new_comm(char *buf, prof_attrs** communicator, MPI_Comm comm, MPI_Comm* newcomm){
    size_t length;
    int comm_size;
    if ( buf == NULL || communicator == NULL ||
         comm == MPI_COMM_NULL || newcomm == NULL){
        mcpt_abort("Newcomm called with NULL\n");
    }
    PMPI_Comm_size(*newcomm, &comm_size);
    length = strlen((*communicator)->name);
    strcpy(&(*communicator)->name[length], buf);
    (*communicator)->bytes = 0;
    (*communicator)->msgs = 0;
    (*communicator)->size = comm_size;
    local_comms[local_cid] = *communicator;
    communicators[my_coms] = *newcomm;
    my_coms++;
    local_cid++;
    return;
}

void
profile_this(prof_attrs *communicator,int index){
    return;
}

int
_MPI_Init(int *argc, char ***argv){
    int ret,rank,size;
    int i,rc;
    prof_attrs *communicator;
    ret = PMPI_Init(argc, argv);
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    /* local_data = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4); */
    local_comms = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4);
    for ( i =0 ; i<size*4; i++ ){
        communicators[i] = MPI_COMM_NULL;
        /* local_data[i] = NULL; */
        local_comms[i] = NULL;
    }
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*1024);
        appname = get_appname();
        printf("MPI_Init: MPI Communicator Profiling Tool\nProfiling application\
 %s\n",appname);
        fflush(stdout);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    if ( communicator == NULL ){
        mcpt_abort("malloc failed at line %s\n",__LINE__);
    }
    strcpy(communicator->name,"W");
    communicator->bytes = 0;
    communicator->size = size;
    communicator->msgs = 0;
    rc = PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    local_comms[local_cid] = communicator;
    local_cid++;
    if ( rc != MPI_SUCCESS ){
        mcpt_abort("Comm_set_attr failed at line %s\n",__LINE__);
    }
    communicators[0] = MPI_COMM_WORLD;
    /* printf("Rank: %d Return from Init\n",rank); */
    /* fflush(stdout); */
    return ret;
}

static int
_MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
    int ret,rank,size;
    int i,rc;
    prof_attrs *communicator;
    ret = PMPI_Init_thread(argc, argv, required, provided);
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    local_data = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4);
    local_comms = (prof_attrs**) malloc (sizeof(prof_attrs*)*size*4);
    for ( i =0 ; i<size*4; i++ ){
        communicators[i] = MPI_COMM_NULL;
        local_data[i] = NULL;
        local_comms[i] = NULL;
    }
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*1024);
        appname = get_appname();
        printf("MPI_Init_thread: MPI Communicator Profiling Tool\nProfiling\
 application %s\n",appname);
        fflush(stdout);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    if ( communicator == NULL ){
        mcpt_abort("malloc failed at line %s\n",__LINE__);
    }
    strcpy(communicator->name,"W");
    communicator->bytes = 0;
    communicator->size = size;
    communicator->msgs = 0;
    local_comms[local_cid] = communicator;
    local_cid++;
    rc = PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    if ( rc != MPI_SUCCESS ){
        mcpt_abort("Comm_set_attr failed at line %s\n",__LINE__);
    }
    communicators[0] = MPI_COMM_WORLD;
    /* printf("Rank: %d Return from Init_thread\n",rank); */
    /* fflush(stdout); */
    return ret;
}


int
MPI_Init_thread(int *argc, char ***argv, int required, int *provided){
    return _MPI_Init_thread(argc, argv, required, provided);
}

void
F77_MPI_INIT_THREAD (int *required, int *provided, int *ierr){
    char **tmp;
    int ret;
    getProcCmdLine (&ac, av);
    tmp = av;
    ret = _MPI_Init_thread(&ac, (char***)&tmp , *required, provided);
    *ierr = ret;
    return;
}

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

int
MPI_Init(int *argc, char ***argv)
{
    return _MPI_Init(argc, argv);
}

int
MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    int ret;
    int i,length,comms;
    char *buf;
    prof_attrs *communicator;
    ret = PMPI_Comm_create(comm, group, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        /* communicators[my_coms] = MPI_COMM_NULL; */
        /* my_coms++; */
        return ret;
    }
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
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
    communicator->size = 0;
    /* int rank; */
    /* PMPI_Comm_rank(MPI_COMM_WORLD, &rank); */
    /* printf("Rank %d, local_cid = %d Func %s\n",rank,local_cid,__FUNCTION__); */
    /* fflush(stdout); */
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    local_data[local_cid] = communicator;
    local_cid++;
    communicators[my_coms] = *newcomm;
    /* num_of_comms++; */
    my_coms++;
    return ret;
}

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

/*
 * The naming scheme for the MPI_Comm_split requires two numbers
 * 1. current number of communicators
 * 2. a unique id for each new communicator via split
 * The prefix of this name is always the parent's name
 */
int
MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    int ret,comms;
    prof_attrs *communicator;
    char *buf;
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
    PMPI_Comm_rank(comm, &rank);
    if ( newcomm== NULL || *newcomm == MPI_COMM_NULL  ){
        return ret;
    }
    /*
     * Call Allreduce in the new communicators and get the minimum rank
     * Note that the rank is the rank from the parent communicator
     * Every new communicator will have an id based on the minimum rank
     * of a process in the parent communicator
     */
    PMPI_Allreduce(&rank, &min_rank, 1, MPI_INT, MPI_MIN, *newcomm);
    /* Use the parent's name as a prefix for the newly created communicator */
    communicator = get_comm_parent(comm);
    buf = (char*) malloc ( sizeof(char)*16);
    /* Suffix of the new communicator with the two ids */
    sprintf(buf,"_s%d.%d",my_coms,min_rank);
    /* Append prefix+suffix and initialize the data for the new communicator */
    new_comm(buf, &communicator, comm, newcomm);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    return ret;
}

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


int
MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
{
    int ret,length,i, comms;
    prof_attrs *communicator;
    char *buf;
    ret = PMPI_Comm_dup(comm, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL )
        return ret;
    /* if ( comm == MPI_COMM_WORLD ){ */
        /* Synchronize */
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, comm);
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
    local_comms[local_cid] = communicator;
    local_cid++;
    communicators[my_coms] = *newcomm;
    /* num_of_comms+=1; */
    my_coms++;
    return ret;
}

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


int
MPI_Cart_create(MPI_Comm old_comm, int ndims, const int *dims,
                const int *periods, int reorder, MPI_Comm *comm_cart)
{

    int ret,length,i,comms;
    prof_attrs *communicator;
    char *buf, *buffer;
    ret = PMPI_Cart_create(old_comm, ndims, dims, periods, reorder, comm_cart);
    PMPI_Allreduce(&my_coms, &comms, 1, MPI_INT, MPI_MAX, old_comm);
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
    /* printf("MPI_cart_create\n"); */
    /* fflush(stdout); */
    PMPI_Comm_set_attr(*comm_cart, namekey(), communicator);
    local_comms[local_cid] = communicator;
    local_cid++;
    printf("MPI Cart comm with name %s and parent %s\n",communicator->name,buffer);
    communicators[my_coms] = *comm_cart;
    /* num_of_comms+=1; */
    my_coms++;
    free(buffer);
    return ret;
}


void
F77_MPI_CART_CREATE(MPI_Fint  * comm_old, int  * ndims, mpip_const_int_t  *dims,
                    mpip_const_int_t  *periods, int  * reorder,
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


int
MPI_Cart_sub(MPI_Comm comm, const int *remain_dims, MPI_Comm *new_comm){
    int ret,length,i, my_rank, newcoms,size;
    prof_attrs *communicator;
    char *buf;
    int dims[MAX_DIMS],id,min_rank;
    ret = PMPI_Cart_sub(comm, remain_dims, new_comm);
    PMPI_Comm_rank(comm, &my_rank);
    PMPI_Comm_size(*new_comm, &size);
    PMPI_Allreduce(&my_coms, &id, 1, MPI_INT, MPI_MAX, comm);
    memset(dims, 0, sizeof(int)*MAX_DIMS);
    PMPI_Cartdim_get(*new_comm, dims);
    newcoms = 1;
    for ( i =0; i<MAX_DIMS; i++ ){
        if ( dims[i] != 0 ){
            newcoms = newcoms*dims[i];
        }
    }
    /* recvbuf = (int*) malloc ( size*sizeof(int) ); */
    /* PMPI_Allgather(&my_rank, 1, MPI_INT, recvbuf, 1, MPI_INT, *new_comm); */
    /* my_coms = newcoms; */
    PMPI_Allreduce(&my_rank, &min_rank, 1, MPI_INT , MPI_MIN, *new_comm);
    communicator = get_comm_parent(comm);
    /* buffer = strdup(communicator->name); */
    buf = (char*) malloc ( sizeof(char)*16);
    sprintf(buf,"_b%d.%d",id,min_rank);
    length = strlen(communicator->name);
    for ( i =0; i<strlen(buf); i++ ){
        communicator->name[length+i] = buf[i];
    }
    communicator->bytes = 0;
    communicator->msgs = 0;
    PMPI_Comm_set_attr(*new_comm, namekey(), communicator);
    local_comms[local_cid] = communicator;
    local_cid++;
    communicators[my_coms] = *new_comm;
    /* num_of_comms+=newcoms; */
    my_coms++;
    return ret;
}

void
F77_MPI_CART_SUB(MPI_Fint  * comm, mpip_const_int_t  *remain_dims,
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


int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Isend on wrong communicator\n");
    PMPI_Type_size(datatype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + count * size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs + 1;
    }
    return ret;
}

void
F77_MPI_ISEND(mpip_const_void_t  *buf, int  * count, MPI_Fint  * datatype,
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



int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Send on wrong communicator\n");
    PMPI_Type_size(datatype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + count * size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}

void
F77_MPI_SEND(mpip_const_void_t  *buf, int  * count, MPI_Fint  * datatype,
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



int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Sendrecv on wrong communicator\n");
    PMPI_Type_size(sendtype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sendcount * size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}


void
F77_MPI_SENDRECV(mpip_const_void_t  *sendbuf, int  * sendcount,MPI_Fint  * sendtype,
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


int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Bcast on wrong communicator\n");
    PMPI_Type_size(datatype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = (sum + count * size);
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;

}

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


int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Allreduce on wrong communicator\n");
    PMPI_Type_size(datatype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + (count * size)*comm_size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}

void
F77_MPI_ALLREDUCE(mpip_const_void_t  *sendbuf, void  *recvbuf, int  * count,
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




int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Allgather on wrong communicator\n");
    PMPI_Type_size(sendtype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + (sendcount * size)*comm_size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
        /* printf("Send bytes %d to %s\n",(sendcount * size)*comm_size,communicator->name); */
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}

int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Alltoall on wrong communicator\n");
    PMPI_Type_size(sendtype, &size);
    PMPI_Comm_size(comm, &comm_size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sendcount*size*comm_size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + comm_size;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}

int
MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
              const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
              const int *recvcounts, const int *rdispls, MPI_Datatype recvtype,
              MPI_Comm comm)
{

    int ret,i,j,flag,size,comm_size,sent;
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
    /* recvd = 0; */
    /* for ( i = 0; i<comm_size; i++ ){ */
    /*     recvd += recvcounts[i]; */
    /* } */
    /* recvd = recvd*size; */

    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Alltoallv on wrong communicator\n");
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sent;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + j;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}

int
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
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Allgatherv on wrong communicator\n");
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + recvd;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + j;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}


int
MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
           MPI_Op op, int root, MPI_Comm comm)
{
    int ret,i,flag,size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;

    ret = PMPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Reduce on wrong communicator\n");
    PMPI_Type_size(datatype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + count*size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;
}

int
MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
           int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
    int ret,i,flag,size;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    flag = 0;
    for ( i=0; i< my_coms; i++){
        if ( comm == communicators[i] )
            break;
    }
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    for ( i =0; i< local_cid; i++ ){
        if ( strcmp(communicator->name, local_comms[i]->name) == 0 )
            break;
    }
    if ( i == local_cid  )
        mcpt_abort("Gather on wrong communicator\n");
    PMPI_Type_size(sendtype, &size);
    if ( flag ){
        sum = communicator->bytes;
        sum = sum + sendcount*size;
        communicator->bytes = sum;
        communicator->msgs = communicator->msgs + 1;
        local_comms[i]->bytes = sum;
        local_comms[i]->msgs = communicator->msgs;
    }
    return ret;

}

int
MPI_Comm_free(MPI_Comm *comm){
    int ret,flag,i;
    prof_attrs *com_info;
    /* printf("CALLING SPLIT\n"); */
    /* fflush(stdout); */
    PMPI_Comm_get_attr(*comm, namekey(), &com_info, &flag);
    if ( flag ){
        for ( i = 0; i<local_cid; i++ ){
            if ( strcmp(com_info->name, local_comms[i]->name) == 0 )
            break;
        }
        if ( i == local_cid  )
            mcpt_abort("Comm_free on wrong communicator\n");
        /* local_data[num_of_local] = (prof_attrs*) malloc (sizeof(prof_attrs)); */
        /* local_data[num_of_local]->bytes = com_info->bytes; */
        /* local_data[num_of_local]->msgs = com_info->msgs; */
        /* local_data[num_of_local]->size = com_info->size; */
        /* strcpy(local_data[num_of_local]->name,com_info->name); */
        local_comms[i] = (prof_attrs*) malloc (sizeof(prof_attrs));
        local_comms[i]->bytes = com_info->bytes;
        local_comms[i]->msgs = com_info->msgs;
        local_comms[i]->size = com_info->size;
        strcpy(local_comms[i]->name,com_info->name);
        /* local_comms[i] = local_data[num_of_local]; */
        num_of_local++;
    }
    /* for ( i = 0; i<num_of_comms; i++ ){ */
    /*     if ( *comm  == communicators[i]) */
    /*         break; */
    /* } */
    /* communicators[i]=MPI_COMM_NULL; */
    /* free(local_comms[i]); */
    ret = PMPI_Comm_free(comm);
    return ret;
}

static int
_Finalize(){
    FILE *fp;
    prof_attrs *array;
    int rank,size;
    int i,j,k,found;
    prof_attrs *recv_buffer;
    prof_attrs dummy;
    char **names, **unames;
    int total_comms,total,num_of_comms;
    uint64_t *bytes, *ubytes;
    uint32_t *msgs, *umsgs;
    PMPI_Comm_rank(MPI_COMM_WORLD, &rank);
    PMPI_Comm_size(MPI_COMM_WORLD, &size);
    num_of_comms = my_coms;
    PMPI_Allreduce(&num_of_comms, &total_comms, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    num_of_comms = total_comms;
    /* if ( rank == 0 ){ */
        /* printf("Num of REAL comms = %d\n",num_of_comms); */
        /* fflush(stdout); */
    /* } */
    array =(prof_attrs*) malloc(sizeof(prof_attrs)*num_of_comms);
    recv_buffer = (prof_attrs*) malloc (sizeof(prof_attrs)*num_of_comms*size);
    MPI_Datatype types[4] = { MPI_CHAR,MPI_UINT64_T, MPI_UINT32_T, MPI_INT };
    int blocklen[4] = {NAMELEN,1,1,1};
    MPI_Aint displacements[4];
    MPI_Aint base_address;
    MPI_Datatype profiler_data;
    PMPI_Get_address(&dummy, &base_address);
    PMPI_Get_address(&dummy.name[0], &displacements[0]);
    /* MPI_Get_address(&dummy.parent[0], &displacements[1]); */
    /* MPI_Get_address(&dummy.prim[0], &displacements[2]); */
    PMPI_Get_address(&dummy.bytes, &displacements[1]);
    PMPI_Get_address(&dummy.msgs, &displacements[2]);
    PMPI_Get_address(&dummy.size, &displacements[3]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    displacements[3] = MPI_Aint_diff(displacements[3], base_address);
    /* displacements[4] = MPI_Aint_diff(displacements[4], base_address); */

    PMPI_Type_create_struct(4, blocklen, displacements, types, &profiler_data);
    PMPI_Type_commit(&profiler_data);
    k = 0;
    for ( i = 0; i < num_of_comms; i++ ){
        if ( local_comms[i] != NULL ){
            /* printf("RANK %d %p %d\n",rank,communicators[i],namekey()); */
            /* fflush(stdout); */
            /* PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag); */
            /* com_info = local_comms[i]; */
            /* if ( rank == 1 || rank == 3 || rank == 5){ */
            /*     printf("RANK %d : i = %d COMM %s bytes = %lu, Msgs = %u\n",rank,i, */
            /*            local_comms[i]->name,local_comms[i]->bytes,local_comms[i]->msgs); */
            /*     fflush(stdout); */
            /* } */

            strcpy(array[i].name, local_comms[i]->name);
            array[i].bytes = local_comms[i]->bytes;
            array[i].msgs= local_comms[i]->msgs;
            array[i].size = local_comms[i]->size;
            /* if ( flag ){ */
                /* strcpy(array[i].name, com_info->name); */
                /* array[i].bytes = com_info->bytes; */
                /* array[i].msgs= com_info->msgs; */
                /* array[i].size = com_info->size; */
            /* } */
        }
        else{
            /* if ( num_of_local > 0 && k < num_of_local){ */
            /*     strcpy(array[i].name, local_data[k]->name); */
            /*     array[i].bytes = local_data[k]->bytes; */
            /*     array[i].size = local_data[k]->size; */
            /*     array[i].msgs = local_data[k]->msgs; */
            /*     k++; */
            /* } */
            /* else{ */
                strcpy(array[i].name, "NULL");
                array[i].msgs = 0;
                array[i].bytes = 0;
                array[i].size = 0;
            /* } */
        }
    }
    /* printf ( "RANK %d i = %d\n",rank,i ); */
    /* fflush(stdout); */
    /* for ( i =0; i< num_of_comms; i++ ){ */
    /*     printf("RANK %d: comm =%s bytes =%lu\n",rank,array[i].name,array[i].bytes); */
    /* } */
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
            /* if ( strcmp(names[i], "W") != 0 ){ */
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
            /* } */
        }
        /* strcpy(unames[j], "W"); */
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


 int
MPI_Finalize (void)
{
  int rc = 0;

  rc = _Finalize ();

  return rc;
}

 void
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
