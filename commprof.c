#include "atom.h"
#include "table.h"
#include "utils.h"
#include <mpi.h>
#include <mem.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>


/* typedef struct _node{ */
/*     char *name; */
/*     struct _node *next; */
/* }node; */


/* node* head; */

const char *comm_name = "comm";
MPI_Comm *communicators = NULL;
int num_of_comms;

/* void newNode(const char *n){ */
/*     node *ptr; */
/*     ptr = ALLOC(sizeof(node)); */
/*     ptr->name = ALLOC(80); */
/*     ptr->next = head->next; */
/*     head->next = ptr; */
/*     ptr->name=strcpy(ptr->name,n); */
/* } */

static int namedel(MPI_Comm comm, int keyval, void *attr, void *s) {
  prof_attrs *com = (prof_attrs*)attr;
  if ( comm == NULL )
    return -1;
  /* if ( com->name ) */
  /*   free(com->name); */
  /* if(com->prim) */
  /*   free(com->prim); */
  /* free(com); */
  /* if ( com->name ) */
  /*     FREE(com->name); */
  /* if (com->prim) */
  /*     FREE(com->prim); */

  FREE(com);
  return 0;
}

static int namekey() {
  // hidden key value for type attributes
  static int namekeyval = MPI_KEYVAL_INVALID;

  if (namekeyval == MPI_KEYVAL_INVALID) {
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN,namedel,&namekeyval,NULL);
  }

  return namekeyval;
}

void apply_print(const void *key,void **value,void *cl){
    prof_attrs *val;
    unsigned long long bytes;
    char name[32], prim[32];
    val = (prof_attrs*)value;
    strcpy(name, val->name);
    bytes = val->bytes;
    strcpy(prim, val->prim);
    printf("Name: %s created by %s bytes transferred %llu\n",name,prim,bytes);
}

int compare(const void *x, const void *y) {
    prof_attrs *a = (prof_attrs*) x;
    prof_attrs *b = (prof_attrs*) y;
    char bufa[32], bufb[32];
    strcpy(bufa,a->name);
    strcpy(bufb,b->name);
    printf("Compare %s with %s\n",bufa,bufb);
    return strcmp(a->name, b->name);
}

/* void print_entry(commtor* communicator){ */
/*     printf("Name: %s, Bytes  %ld, Prim created %s\n",communicator->name, */
/*            communicator->bytes,communicator->prim); */
/* } */

int MPI_Init(int *argc, char ***argv){
    int ret,rank,size;
    prof_attrs *communicator;
    /* unsigned long long bytes; */
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* comm_table = ALLOC(sizeof(prof_attrs)*64); */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators = ALLOC(sizeof(MPI_Comm)*size);
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (malloc(sizeof(char)*256));
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    NEW(communicator);
    strcpy(communicator->name,"WORLD");
    communicator->bytes = 0;
    strcpy(communicator->prim , "INIT");
    /* communicator->index = num_of_comms; */
    communicator->size = size;
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[num_of_comms] = MPI_COMM_WORLD;

    /* comm_table[num_of_comms]=communicator; */
    /* Table_put(table, MPI_COMM_WORLD, communicator); */
    num_of_comms++;
    return ret;
}

int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
    int ret,size;
    prof_attrs *communicator;
    ret = PMPI_Comm_create(comm, group, newcomm);
    NEW(communicator);
    sprintf(communicator->name,"%s%d",comm_name,num_of_comms);
    /* strcpy(buf, comm_name); */
    communicator->bytes = 0;
    strcpy(communicator->prim,"COMM_CREATE");
    /* communicator->name = buf; */
    /* communicator->prim = prim; */
    MPI_Comm_size(*newcomm, &size);
    communicator->size = size;
    /* MPI_Comm_set_name(*newcomm, buf); */
    printf("New Comm Entry %s\n",communicator->name);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[num_of_comms] = *newcomm;
    /* comm_table[num_of_comms]=communicator; */
    /* Table_put(table, *newcomm, communicator); */
    /* MPI_Comm_get_name(*newcomm, buffer, &len); */
    /* if ( len && buffer){ */
    /*     printf ( "Name = %.*s\n",len,buffer ); */
    /* } */
    num_of_comms++;
    return ret;
}

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm){
    int ret,size,flag;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    flag = 0;
    /* communicator = Table_get(table,comm); */
    PMPI_Type_size(datatype, &size);
    PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag);
    if ( flag ){
        printf("MPI_Send to communicator %s\n",communicator->name);
        sum = communicator->bytes;
        sum = sum + count * size;
        communicator->bytes = sum;
        /* Table_put(table, comm, communicator); */
    }
    /* We need to check in case this is invalid */

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
    /* FILE *fp; */
    /* fp = fopen("profiler_stats.txt","w"); */
    prof_attrs *array;
    int rank,i,j,index,size,flag;
    prof_attrs *com_info;
    prof_attrs *recv_buffer;
    prof_attrs dummy;
    /* char **names, **names_buf; */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* names = (char**) malloc ( sizeof (char*) *num_of_comms ); */
    /* for (i = 0; i < num_of_comms; i++) { */
    /*     names[i] = (char*) malloc(32); */
    /* } */
    /* names_buf = (char**) malloc (sizeof(char*)*num_of_comms*size); */

    /* for (i = 0; i < num_of_comms*size; i++) { */
    /*     names_buf[i] = (char*) malloc(32); */
    /* } */
    array =(prof_attrs*) malloc(sizeof(prof_attrs)*num_of_comms);
    recv_buffer = (prof_attrs*) malloc (sizeof(prof_attrs)*num_of_comms*size);
    /* for ( i =0; i<num_of_comms*size; i++ ){ */
    /*     recv_buffer[i]=(prof_attrs*) malloc ( sizeof(prof_attrs) ); */
    /* } */
    MPI_Datatype types[4] = { MPI_CHAR, MPI_CHAR, MPI_UNSIGNED_LONG_LONG, MPI_INT };
    int blocklen[4] = {32,32,1,1};
    MPI_Aint displacements[4];
    MPI_Aint base_address;
    MPI_Datatype profiler_data;
    MPI_Get_address(&dummy, &base_address);
    MPI_Get_address(&dummy.name[0], &displacements[0]);
    MPI_Get_address(&dummy.prim[0], &displacements[1]);
    MPI_Get_address(&dummy.bytes, &displacements[2]);
    MPI_Get_address(&dummy.size, &displacements[3]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    displacements[3] = MPI_Aint_diff(displacements[3], base_address);

    /* int len = Table_length(table); */
    /* for ( i = 0; i < num_of_comms; i++ ){ */
    /*     } */
    /* } */
    /* array = (prof_attrs**)Table_toArray(table, NULL); */
    /* qsort(array, len, sizeof(*array), compare); */
    /* bytes = ALLOC(sizeof(unsigned long long)*len); */
    /* recv_buf_bytes = ALLOC(sizeof(unsigned long long)*len); */
    /* sizes = ALLOC(sizeof(int)*len); */
    /* recv_buf_sizes = ALLOC(sizeof(int)*len*size); */

    /* Table_map(table, apply_print, NULL); */
    /* printf("Num of communicators = %d\n",num_of_comms); */
    MPI_Type_create_struct(4, blocklen, displacements, types, &profiler_data);
    MPI_Type_commit(&profiler_data);
    for ( i = 0; i < num_of_comms; i++ ){
        /* array = Table_get(table, communicators[i]); */
        PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag);
        if ( flag ){
            /* printf("Rank %d Communicator %s created by %s Bytes %llu\n", */
            /*        rank,com_info->name,com_info->prim,com_info->bytes); */
            /* strcpy(names[i], com_info->name); */
            strcpy(array[i].name, com_info->name);
            strcpy(array[i].prim, com_info->prim);
            array[i].bytes = com_info->bytes;
            array[i].size = com_info->size;
            /* printf("Rank %d Communicator %s created by %s Bytes %llu\n", */
            /*        rank,array[i].name,array[i].prim,array[i].bytes); */
            /* printf("Rank %d Comm Name = %s\n",rank,names[i]); */
        }
    }
    /* PMPI_Gather(names, num_of_comms*32, MPI_CHAR, names_buf, num_of_comms*32, MPI_CHAR, 0, MPI_COMM_WORLD); */
    PMPI_Gather(array, num_of_comms*sizeof(prof_attrs), MPI_BYTE, recv_buffer,
                num_of_comms*sizeof(prof_attrs), MPI_BYTE, 0, MPI_COMM_WORLD);
    if ( rank == 0 ){
        for ( i =0; i<size; i++ ){
            for ( j = 0; j < num_of_comms; j++ ){
                index = i*num_of_comms + j;
                printf("Communicator: %s create from %s bytes = %llu, size = %d\n",
                   recv_buffer[index].name,recv_buffer[index].prim,recv_buffer[index].bytes,
                   recv_buffer[index].size);
            }
        }
    }
    /* for (i = 0; array[i]; i+=1) { */
    /*     com_info = (prof_attrs*)array[i]; */
        /* printf("Rank :%d Communicator %s ",rank,(char*)com_info->name); */
        /* printf("Bytes %llu\n",(unsigned long long)com_info->bytes); */
        /* bytes[i] = (unsigned long long)com_info->bytes; */
        /* sizes[i] = (int)(com_info)->size; */
        /* printf("Rank %d: Communicator %s: Bytes = %llu Size = %d, i = %d\n",rank,com_info->name,bytes[i],sizes[i],i); */
    /* } */
    /* PMPI_Reduce(bytes, recv_buf_bytes, len, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD); */
    /* PMPI_Gather(sizes, 1, MPI_INT, recv_buf_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD); */
    /* if ( rank == 0 ){ */
    /*     for ( i =0; array[i]; i+=1 ){ */
    /*         com_info = (prof_attrs*)array[i]; */
    /*         printf("Communicator %s: Bytes = %llu Size = %d\n",com_info->name,recv_buf_bytes[i],recv_buf_sizes[i]); */
    /*     } */
    /* } */
    /* for ( i=0; i<num_of_comms; i++ ){ */
    /*     printf("Communicator %s: %ld Bytes\n",comm_table[i]->name,comm_table[i]->bytes); */
    /* } */
    /* commtor *com; */
    /* if ( rank == 0 ){ */
        /* array = Table_toArray(table, NULL); */
        /* for (i = 0; array[i]; i+=2) { */
        /*     com = (commtor*)array[i+1]; */
        /*     printf("%ld %s ",*(int64_t*)array[i],(char*)com->name); */
        /*     printf("%ld\n",*(int64_t*)com->bytes); */
        /* } */
        /* for ( i=0; i<num_of_comms; i++ ){ */
        /*     printf("Communicator %s: %ld Bytes\n",comm_table[i]->name,comm_table[i]->bytes); */
        /* } */
    /* } */
    /* fclose(fp); */
    /* if (communicators) */
    /*     FREE(communicators); */
    return PMPI_Finalize();
}
