#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 56
typedef struct _tnode{
    char *name;
    struct _tnode *parent;
    int bytes;
    int msgs;
}node;
char names[N][16] = {"WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.0","WORLD","s_1.1","s_2.1","s_3.0","WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.1","WORLD","s_1.1","s_2.1","s_3.1","WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.0","WORLD","s_1.1","s_2.1","s_3.0","WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.1","WORLD","s_1.1","s_2.1","s_3.1"};
char parents[N][16] = {"NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1","NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1","NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1","NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1"};
uint64_t bytes[N] = {32,64,8,0,16,4,32,0,0,4,0,0,4,0,32,0,0,0,16,0,32,0,0,4,0,0,0,8,32,0,0,0,16,4,32,0,0,0,0,0,4,0,32,0,0,0,16,0,32,0,0,0,0,0,0,0};
uint32_t msgs[N]= {1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0};
int i,j,k;
node **head = NULL;

int main(void){
    int i;
    for ( i = 0; i < N; i++){
        printf("%s\n",names[i]);
    }
    /* if ( strcmp(names[i],"WORLD") == 0 && head == NULL){ */
    /* } */
    return 0;
}
