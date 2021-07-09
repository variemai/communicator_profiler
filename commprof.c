#include "atom.h"
#include "table.h"
#include "utils.h"
#include <mpi.h>
#include <mem.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>


typedef struct _node{
    char *name;
    struct _node *next;
}node;


node* head;

const char *comm_name = "comm";

void newNode(const char *n){
    node *ptr;
    ptr = ALLOC(sizeof(node));
    ptr->name = ALLOC(80);
    ptr->next = head->next;
    head->next = ptr;
    ptr->name=strcpy(ptr->name,n);
}


void apply_print(const void *key,void **value,void *cl){
    commtor *val;
    val = (commtor*)value;
    char name[32], prim[32];
    strcpy(name, val->name);
    int64_t bytes = val->bytes;
    strcpy(prim, val->prim);
    printf("Name: %s created by %s bytes transferred %ld\n",name,prim,bytes);
}

int compare(const void *x, const void *y) {
    return strcmp(*(char **)x, *(char **)y);
}

void print_entry(commtor* communicator){
    printf("Name: %s, Bytes  %ld, Prim created %s\n",communicator->name,
           communicator->bytes,communicator->prim);
}

int MPI_Init(int *argc, char ***argv){
    int ret,rank;
    char *buf;
    commtor *communicator;
    int64_t *bytes;
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    comm_table = ALLOC(sizeof(commtor)*64);
    /* MPI_Comm_size(MPI_COMM_WORLD, &size); */
    table = Table_new(128, NULL, NULL);
    if ( rank == 0 ){
        appname = (malloc(sizeof(char)*256));
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    buf = ALLOC(32);
    sprintf(buf,"%s%d",comm_name,num_of_comms);
    NEW(communicator);
    NEW(bytes);
    *bytes = 0;
    communicator->bytes = bytes;
    communicator->name = buf;
    communicator->prim = strdup("INIT");
    communicator->comm = NULL;
    comm_table[num_of_comms]=communicator;
    Table_put(table, MPI_COMM_WORLD, communicator);
    num_of_comms++;
    return ret;
}

int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
    int ret;
    char *buf;
    commtor *communicator;
    int64_t *bytes;
    ret = PMPI_Comm_create(comm, group, newcomm);
    buf = ALLOC(32);
    sprintf(buf,"%s%d",comm_name,num_of_comms);
    /* strcpy(buf, comm_name); */
    NEW(communicator);
    NEW(bytes);
    *bytes = 0;
    communicator->bytes = bytes;
    communicator->name = buf;
    communicator->prim = strdup("COMM_CREATE");
    /* MPI_Comm_set_name(*newcomm, buf); */
    communicator->comm = *newcomm;
    printf("New Comm Entry %s\n",communicator->name);
    comm_table[num_of_comms]=communicator;
    Table_put(table, *newcomm, communicator);
    /* MPI_Comm_get_name(*newcomm, buffer, &len); */
    /* if ( len && buffer){ */
    /*     printf ( "Name = %.*s\n",len,buffer ); */
    /* } */
    num_of_comms++;
    return ret;
}

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm){
    int ret,size;
    commtor *communicator;
    int64_t sum = 0;
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    communicator = Table_get(table,comm);
    MPI_Type_size(datatype, &size);
    if ( communicator ){
        printf("MPI_Send to communicator %s\n",communicator->name);
        sum = *(communicator->bytes);
        sum = sum + count * size;
        *(communicator->bytes) = sum;
        Table_put(table, comm, communicator);
    }
    /* Table_put(table, comm, void *value) */
    /* for ( i=0; i<num_of_comms; i++ ){ */
    /*     if ( comm == (comm_table[i]->comm)){ */
    /*         comm_table[i]->bytes += count* sizeof(datatype); */
    /*         printf("MPI_Send to communicator %s\n",comm_table[i]->name); */
    /*         break; */
    /*     } */
    /* } */
    return ret;
}

extern int MPI_Finalize(){
    int rank,i;
    void **array;
    commtor *com;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( rank == 0 ){
        array = Table_toArray(table, NULL);
        /* Table_map(table, apply_print, NULL); */
        for (i = 0; array[i]; i+=2) {
            com = (commtor*)array[i+1];
            printf("%ld %s ",*(int64_t*)array[i],(char*)com->name);
            printf("%ld\n",*(int64_t*)com->bytes);
        }
        /* for ( i=0; i<num_of_comms; i++ ){ */
        /*     printf("Communicator %s: %ld Bytes\n",comm_table[i]->name,comm_table[i]->bytes); */
        /* } */
    }
    return PMPI_Finalize();
}
