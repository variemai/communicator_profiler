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


int main(int32_t argc, char const *argv[]){
	uint32_t i;
	const char *strn;
	clock_t begin,end;
	double time_spent;
    char **buf;
    List_T lista;
	strn="paparia";
	begin=clock();
	for(i=0; i<1000; i++){
		Atom_new(strn,7);
	}
	end=clock();
	time_spent=(double)(end-begin)/CLOCKS_PER_SEC;
	printf("%f\n",time_spent);
    buf=ALLOC(10*sizeof(char*));
    buf[0]="PEOS";
    lista=List_list("arxodoa", "malakies", "kala", NULL);
    printf("%d\n",List_length(lista));
    List_map(lista,applyPrint,NULL);
    return 0;

}
