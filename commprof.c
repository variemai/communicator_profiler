#include "utils.h"
#include <mpi.h>
#include <list.h>
#include <mem.h>
#include <string.h>


typedef struct _node{
    char *name;
    struct _node *next;
}node;

node* head;

void newNode(const char *n){
    node *ptr;
    ptr = ALLOC(sizeof(node));
    ptr->name = ALLOC(80);
    ptr->next = head->next;
    head->next = ptr;
    ptr->name=strcpy(ptr->name,n);
}

void applyPrint(void **ptr,void *cl){
    printf("%s\n",(char*)*ptr);
}

extern int MPI_Init(int *argc, char ***argv){
    int ret,rank;
    char **buf;
    List_T lista;
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( rank == 0 ){
        appname = (malloc(sizeof(char)*256));
        appname = getProcExeLink();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
        buf=ALLOC(10*sizeof(char*));
        buf[0]="PEOS";
        lista=List_list("arxodoa", "malakies", "kala", NULL);
        printf("%d\n",List_length(lista));
        List_map(lista,applyPrint,NULL);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return ret;
}
