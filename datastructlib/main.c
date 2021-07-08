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

Except_T FileOpen_failed = {"FILE OPEN FAILED"};

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

void apply_print(const void *key,void **value,void *cl){
    int *val = *value;
    printf("%s %d\n",(char*)key,*val);
}

int compare(const void *x, const void *y) {
    return strcmp(*(char **)x, *(char **)y);
}

int main(int32_t argc, char const *argv[]){
	uint32_t i;
	const char *strn;
	clock_t begin,end;
	double time_spent;
    char **buf;
    const char *com0, *com1, *com2, *com3;
    int* b_com0, *b_com1, *b_com2, *b_com3, *count;
    void** array;
    List_T lista;
    Table_T table;
    table = Table_new(5,NULL,NULL);
    com0=Atom_string("com0");
    com1=Atom_string("com1");
    com2=Atom_string("com2");
    com3=Atom_string("com3");
    NEW(b_com0);
    NEW(b_com1);
    NEW(b_com2);
    NEW(b_com3);
    *b_com0 = 100;
    *b_com1 = 10;
    *b_com2 = 33;
    *b_com3 = 20;
    count = Table_put(table,com0,b_com0);
    if ( count  )
        *count = (*count) + (*b_com0);
    count = Table_put(table,com1,b_com1);
    if ( count  )
        *count = (*count) + (*b_com1);
    count = Table_put(table,com2,b_com2);
    if ( count  )
        *count = (*count) + (*b_com2);
    count =Table_put(table,com3,b_com3);
    if ( count  )
        *count = (*count) + (*b_com3);
    count = Table_put(table,com0,b_com0);
    if ( count  )
        *count = (*count) + (*b_com0);
    count = Table_put(table,com0,b_com0);
    if ( count  )
        *count = (*count) + (*b_com0);
    Table_map(table,apply_print,NULL);
    array = Table_toArray(table, NULL);
    printf("Table length = %d\n",Table_length(table));
    qsort(array, Table_length(table), 2*sizeof (*array),
       compare);
   for (i = 0; array[i]; i += 2)
       printf("%d\t%s\n", *(int *)array[i+1],
              (char *)array[i]);
	strn="paparia";
	begin=clock();
	for(i=0; i<1000; i++){
		Atom_new(strn,7);
	}
	end=clock();
	time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
	/* printf("%f\n",time_spent); */
    buf=ALLOC(10*sizeof(char*));
    buf[0]="PEOS";
    lista=List_list("arxodoa", "malakies", "kala", NULL);
    /* printf("%d\n",List_length(lista)); */
    /* List_map(lista,applyPrint,NULL); */
    return 0;

}
