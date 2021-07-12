#include "atom.h"
#include "table.h"
#include "utils.h"
#include <mpi.h>
#include <mem.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>


/* typedef struct _node{ */
/*     char *name; */
/*     struct _node *next; */
/* }node; */


/* node* head; */

const char *comm_name = "comm";
MPI_Comm *communicators = NULL;

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
  if ( com->name )
      FREE(com->name);
  if (com->prim)
      FREE(com->prim);

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
    return strcmp(*(char **)x, *(char **)y);
}

/* void print_entry(commtor* communicator){ */
/*     printf("Name: %s, Bytes  %ld, Prim created %s\n",communicator->name, */
/*            communicator->bytes,communicator->prim); */
/* } */

int MPI_Init(int *argc, char ***argv){
    int ret,rank;
    char *buf, *prim;
    prof_attrs *communicator;
    /* unsigned long long bytes; */
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* comm_table = ALLOC(sizeof(prof_attrs)*64); */
    communicators = ALLOC(sizeof(MPI_Comm)*64);
    /* MPI_Comm_size(MPI_COMM_WORLD, &size); */
    table = Table_new(128, NULL, NULL);
    if ( rank == 0 ){
        appname = (malloc(sizeof(char)*256));
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    buf = ALLOC(32);
    strcpy(buf, "WORLD");
    NEW(communicator);
    communicator->bytes = 0;
    communicator->name = buf;
    prim = ALLOC(32);
    strcpy(prim, "INIT");
    communicator->prim = prim;
    communicator->comm = MPI_COMM_WORLD;
    communicator->index = num_of_comms;
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[num_of_comms] = MPI_COMM_WORLD;

    /* comm_table[num_of_comms]=communicator; */
    Table_put(table, MPI_COMM_WORLD, communicator);
    num_of_comms++;
    return ret;
}

int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
    int ret;
    char *buf,*prim;
    prof_attrs *communicator;
    unsigned long long bytes;
    ret = PMPI_Comm_create(comm, group, newcomm);
    buf = ALLOC(32);
    sprintf(buf,"%s%d",comm_name,num_of_comms);
    /* strcpy(buf, comm_name); */
    NEW(communicator);
    bytes = 0;
    prim = ALLOC(32);
    strcpy(prim, "COMM_CREATE");
    communicator->bytes = bytes;
    communicator->name = buf;
    communicator->prim = prim;
    /* MPI_Comm_set_name(*newcomm, buf); */
    communicator->comm = *newcomm;
    printf("New Comm Entry %s\n",communicator->name);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[num_of_comms] = *newcomm;
    /* comm_table[num_of_comms]=communicator; */
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
    void **array;
    int rank,i,flag;
    prof_attrs *com_info;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* for ( i = 0; i < num_of_comms; i++ ){ */
    /*     PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag); */
    /*     if ( flag ){ */
    /*         printf("Rank %d Communicator %s Bytes %llu\n",rank,com_info->name,com_info->bytes); */
    /*     } */
    /* } */
    array = Table_toArray(table, NULL);
    for (i = 0; array[i]; i+=2) {
        com_info = (prof_attrs*)array[i+1];
        printf("Rank :%d Communicator %s ",rank,(char*)com_info->name);
        printf("Bytes %llu\n",(unsigned long long)com_info->bytes);
    }
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
    if (communicators)
        FREE(communicators);
    return PMPI_Finalize();
}
