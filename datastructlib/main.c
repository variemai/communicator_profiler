#include <stdio.h>
#include <assert.h>
#include "arith.h"
#include <stdlib.h>
#include <string.h>
#include "except.h"
#include "mem.h"
#include "stack.h"
#include "queue.h"
#include "list.h"
#include "atom.h"
#include <time.h>
#include "table.h"
#include <mpi.h>


Except_T FileOpen_failed = {"FILE OPEN FAILED"};

typedef struct _node{
    char *name;
    struct _node *next;
}node;

typedef struct fuck_my_life{
    const char* key;
    int* value;
}fml;

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

void apply_print(const void *key,void **value,void *cl){
    int *val = *value;
    printf("%s %d\n",(char*)key,*val);
}

int compare(const void *x, const void *y) {
    return strcmp(*(char **)x, *(char **)y);
}

int main(int32_t argc, char const *argv[]){
	/* uint32_t i; */
	/* const char *strn; */
	/* clock_t begin,end; */
	/* double time_spent; */
    /* char **buf; */
    int length,i;
    const char *com0, *com1, *com2, *com3;
    int* b_com0, *b_com1, *b_com2, *b_com3;
    void** array;
    /* List_T lista; */
    Table_T table;
    MPI_Comm newcomm,split,third;
    MPI_Request req0, req1, req2;
    int rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    table = Table_new(5,NULL,NULL);
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &newcomm);
    MPI_Comm_rank(newcomm, &rank);
    NEW(b_com0);
    NEW(b_com1);
    NEW(b_com2);
    NEW(b_com3);
    if ( rank == 0 ){
        MPI_Isend(&rank, 1, MPI_INT, 1, 0, newcomm, &req0);
        Table_put(table, &req0, newcomm);
    }
    if ( rank == 1 ){
        MPI_Irecv(b_com3, 1, MPI_INT, 0, 0, newcomm, &req0);
        Table_put(table, &req0, newcomm);
    }
    array = Table_toArray(table, NULL);
    for (i = 0; array[i]; i += 2){
        printf("%d: %p\t%p\n", rank,(MPI_Request)array[i],(MPI_Comm)array[i+1]);
        fflush(stdout);
    }
    MPI_Wait(&req0, MPI_STATUS_IGNORE);
    com0=Atom_string("com0");
    com1=Atom_string("com1");
    com2=Atom_string("com2");
    com3=Atom_string("com3");
    *b_com0 = 100;
    *b_com1 = 10;
    *b_com2 = 33;
    *b_com3 = 20;
    MPI_Comm_split(MPI_COMM_WORLD, rank %2, rank, &split);
    MPI_Comm_rank(split, &rank);
    if ( rank == 0 ){
        MPI_Isend(&rank, 1, MPI_INT, 1, 0, split, &req1);
    }
    if ( rank == 1 ){
        MPI_Irecv(b_com3, 1, MPI_INT, 0, 0, split, &req1);
    }
    MPI_Wait(&req1, MPI_STATUS_IGNORE);
    Table_put(table, &req1, split);
    /* count = Table_put(table,com0,b_com0); */
    /* if ( count  ) */
    /*     *count = (*count) + (*b_com0); */
    /* count = Table_put(table,com1,b_com1); */
    /* if ( count  ) */
    /*     *count = (*count) + (*b_com1); */
    /* count = Table_put(table,com2,b_com2); */
    /* if ( count  ) */
    /*     *count = (*count) + (*b_com2); */
    /* count =Table_put(table,com3,b_com3); */
    /* if ( count  ) */
    /*     *count = (*count) + (*b_com3); */
    /* count = Table_put(table,com0,b_com0); */
    /* if ( count  ) */
    /*     *count = (*count) + (*b_com0); */
    /* count = Table_put(table,com0,b_com0); */
    /* if ( count  ) */
    /*     *count = (*count) + (*b_com0); */
    /* Table_remove(table, com2); */
    /* Table_map(table,apply_print,NULL); */
    /* printf("%s\n",com2); */
    req2 = MPI_REQUEST_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &third);
    Table_put(table, &req2, third);
    Table_remove(table, &req1);
    length = Table_length(table);
    printf("Table length = %d\n",Table_length(table));
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* qsort(array, Table_length(table), 2*sizeof (*array), */
    /*    compare); */
    if ( rank == 0 )
      printf("%p\n", MPI_REQUEST_NULL);
    MPI_Barrier(MPI_COMM_WORLD);
	/* strn="paparia"; */
	/* begin=clock(); */
	/* for(i=0; i<1000; i++){ */
	/* 	Atom_new(strn,7); */
	/* } */
	/* end=clock(); */
	/* time_spent=(double)(end-begin)/CLOCKS_PER_SEC; */
	/* printf("%f\n",time_spent); */
    /* buf=ALLOC(10*sizeof(char*)); */
    /* buf[0]="PEOS"; */
    /* lista=List_list("arxodoa", "malakies", "kala", NULL); */
    /* printf("%d\n",List_length(lista)); */
    /* List_map(lista,applyPrint,NULL); */
   FREE(array);
   /* FREE(b_com0); */
   /* FREE(b_com1); */
   /* FREE(b_com2); */
   /* FREE(b_com3); */
   MPI_Finalize();
   return 0;

}
