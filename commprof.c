#include "utils.h"
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <ctype.h>
#include <unistd.h>

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
int num_of_comms;
int line_called;

static int namedel(MPI_Comm comm, int keyval, void *attr, void *s) {
  prof_attrs *com = (prof_attrs*)attr;
  if ( comm == NULL )
      return -1;
  else
      free(com);
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

extern int MPI_Init(int *argc, char ***argv){
    int ret,rank,size,i;
    prof_attrs *communicator;
    /* unsigned long long bytes; */
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    /* table = Table_new(128, NULL, NULL); */
    if ( rank == 0 ){
        appname = (char*)malloc(sizeof(char)*256);
        appname = get_appname();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    strcpy(communicator->name,"WORLD");
    communicator->bytes = 0;
    strcpy(communicator->prim , "INIT");
    /* communicator->index = num_of_comms; */
    /* communicator->size = size; */
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[num_of_comms] = MPI_COMM_WORLD;

    /* Table_put(table, MPI_COMM_WORLD, communicator); */
    num_of_comms++;
    return ret;
}

extern int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
    int ret,size,i,flag,rank;
    prof_attrs *communicator, *com_info;
    char *name;
    ret = PMPI_Comm_create(comm, group, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        /* communicators[num_of_comms] = NULL; */
        /* num_of_comms++; */
        return ret;
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    /* //append mother's name if mother != WORLD */
    /* if ( comm != MPI_COMM_WORLD ){ */
    /*     //find the mother */
    /*     for ( i = 0; i<num_of_comms; i++ ){ */
    /*         if ( comm  == communicators[i] ) */
    /*             break; */
    /*     } */
    /*     if ( i == num_of_comms ){ */
    /*         // Should have found the mother */
    /*         fprintf(stderr, "Error: could not find the mother of communicator.\nAborting\n"); */

    /*     } */
    /*     else{ */
    /*         PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag); */
    /*         if ( flag ){ */
    /*             strncpy(name, com_info->name, 32); */
    /*             printf("Mother: %s\n",name); */
    /*             sprintf(communicator->name, "%s-%d",name,__LINE__); */
    /*             printf("New comm with name %s\n",communicator->name); */
    /*         } */
    /*         else{ */
    /*             fprintf(stderr, "Flag invalid\nAborting\n"); */
    /*         } */
    /*     } */

    /* } */
    /* else */
    sprintf(communicator->name,"%c%d",'c',__LINE__);
    /* /\* strcpy(buf, comm_name); *\/ */
    communicator->bytes = 0;
    strcpy(communicator->prim,"COMM_CREATE");
    /* /\* communicator->name = buf; *\/ */
    /* /\* communicator->prim = prim; *\/ */
    /* MPI_Comm_size(*newcomm, &size); */
    /* communicator->size = size; */
    /* /\* MPI_Comm_set_name(*newcomm, buf); *\/ */
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
    /* printf("%p",__builtin_return_address(0)); */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
    /* printf("RANK %d: Num of Comms After create = %d\n",rank,num_of_comms); */
    return ret;
}

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm){
    int ret,size,rank;
    prof_attrs *communicator;
    ret = PMPI_Comm_split(comm, color, key, newcomm);
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    sprintf(communicator->name,"%c%d_%d",'s',line_called,color);
    /* /\* strcpy(buf, comm_name); *\/ */
    /* communicator->bytes = 0; */
    strcpy(communicator->prim,"COMM_SPLIT");
    /* /\* communicator->name = buf; *\/ */
    /* /\* communicator->prim = prim; *\/ */
    /* MPI_Comm_size(*newcomm, &size); */
    /* communicator->size = size; */
    /* /\* MPI_Comm_set_name(*newcomm, buf); *\/ */
    printf("New Comm Entry %s\n",communicator->name);
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[num_of_comms] = *newcomm;
    /* comm_table[num_of_comms]=communicator; */
    /* Table_put(table, *newcomm, communicator); */
    /* MPI_Comm_get_name(*newcomm, buffer, &len); */
    /* if ( len && buffer){ */
    /*     printf ( "Name = %.*s\n",len,buffer ); */
    /* } */
    num_of_comms+=color+1;
    /* printf("%p",__builtin_return_address(0)); */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
    /* printf("RANK %d: Num of Comms After split = %d\n",rank,num_of_comms); */
    return ret;
}


extern int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm){
    int ret,size,flag;
    prof_attrs  *communicator;
    unsigned long long  sum = 0;
    ret = PMPI_Send(buf, count, datatype, dest, tag, comm);
    flag = 0;
    /* communicator = Table_get(table,comm); */
    /* PMPI_Type_size(datatype, &size); */
    /* PMPI_Comm_get_attr(comm, namekey(), &communicator, &flag); */
    /* if ( flag ){ */
    /*     /\* printf("MPI_Send to communicator %s\n",communicator->name); *\/ */
    /*     sum = communicator->bytes; */
    /*     sum = sum + count * size; */
    /*     communicator->bytes = sum; */
    /*     /\* Table_put(table, comm, communicator); *\/ */
    /* } */
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
    /* prof_attrs *array; */
    int rank,i,j,k,size,flag;
    prof_attrs *com_info;
    /* prof_attrs *recv_buffer; */
    /* prof_attrs dummy; */
    /* char **names, **snames; */
    int found, comm_num,mycom;
    int *total_comms;
    /* unsigned long long *bytes; */
    char **names, **names_buf;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* free(communicators); */

    /* Before we do anything we should first learn the number of communicators */
    total_comms = (int*)malloc(sizeof(int)*size);
    mycom = num_of_comms-1;
    PMPI_Gather(&num_of_comms, 1, MPI_INT, total_comms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ( rank == 0 ){
    for ( i = 0; i< size; i++ ){
        if ( total_comms[i] > num_of_comms ){
            num_of_comms = total_comms[i];
        }
    }
    printf( "Num of comms = %d\n",num_of_comms);
    }

    /* names = (char**)malloc(sizeof(char*)*num_of_comms); */
   /* snames = (char**)malloc(sizeof(char*)*num_of_comms*size); */
   /* for ( i =0; i<num_of_comms; i++ ){ */
   /*     names[i] = (char*)malloc(32); */
   /* /\*     snames[i] = (char*)malloc(32); *\/ */
   /* } */
   /* for ( i =0; i<mycom; i++ ){ */
   /*     if ( communicators[i] != NULL ){ */
   /*          PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag); */
   /*          if ( flag ){ */
   /*              printf("%s\n",com_info->name); */
   /*              /\* strcpy(com_info->name, names[i]); *\/ */
   /*          } */
   /*     } */
   /* } */
   /* for ( i =0; i<num_of_comms; i++ ){ */
   /*     free(names[i]); */
   /* } */
   /* free(names); */
   /* free(total_comms); */
   /* free(communicators); */
   /* names_buf = (char**)malloc(sizeof(char*)*num_of_comms*size); */
   /* for ( i =0; i<num_of_comms*size; i++ ){ */
   /*     names_buf[i] = (char*)malloc(32); */
   /* } */


   /*  array =(prof_attrs*) malloc(sizeof(prof_attrs)*num_of_comms); */
   /*  recv_buffer = (prof_attrs*) malloc (sizeof(prof_attrs)*num_of_comms*size); */
   /*  MPI_Datatype types[4] = { MPI_CHAR, MPI_CHAR, MPI_UNSIGNED_LONG_LONG, MPI_INT }; */
   /*  int blocklen[4] = {32,32,1,1}; */
   /*  MPI_Aint displacements[4]; */
   /*  MPI_Aint base_address; */
   /*  MPI_Datatype profiler_data; */
   /*  MPI_Get_address(&dummy, &base_address); */
   /*  MPI_Get_address(&dummy.name[0], &displacements[0]); */
   /*  MPI_Get_address(&dummy.prim[0], &displacements[1]); */
   /*  MPI_Get_address(&dummy.bytes, &displacements[2]); */
   /*  MPI_Get_address(&dummy.size, &displacements[3]); */
   /*  displacements[0] = MPI_Aint_diff(displacements[0], base_address); */
   /*  displacements[1] = MPI_Aint_diff(displacements[1], base_address); */
   /*  displacements[2] = MPI_Aint_diff(displacements[2], base_address); */
   /*  displacements[3] = MPI_Aint_diff(displacements[3], base_address); */

   /*  MPI_Type_create_struct(4, blocklen, displacements, types, &profiler_data); */
   /*  MPI_Type_commit(&profiler_data); */
   /*  for ( i = 0; i < mycom; i++ ){ */
   /*      if ( communicators[i] != NULL ){ */
   /*          PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag); */
   /*          if ( flag ){ */
   /*              strcpy(array[i].name, com_info->name); */
   /*              strcpy(array[i].prim, com_info->prim); */
   /*              array[i].bytes = com_info->bytes; */
   /*              array[i].size = com_info->size; */
   /*          } */
   /*      } */
   /*      else{ */
   /*          strcpy(array[i].name, "NULL"); */
   /*          strcpy(array[i].prim, "NULL"); */
   /*          array[i].bytes = 0; */
   /*          array[i].size = 0; */
   /*      } */
   /*  } */
   /* PMPI_Gather(names, num_of_comms*32, MPI_CHAR, names_buf, num_of_comms*32, MPI_CHAR, 0, MPI_COMM_WORLD); */
    /* PMPI_Gather(array, num_of_comms*sizeof(prof_attrs), MPI_BYTE, recv_buffer, */
    /*             num_of_comms*sizeof(prof_attrs), MPI_BYTE, 0, MPI_COMM_WORLD); */
    /* if ( rank == 0 ){ */
    /*     for ( i =0; i<size; i++ ){ */
    /*         for ( j = 0; j<num_of_comms; j++ ){ */
    /*             printf("%s\n",names_buf[i*size+j]); */
    /*         } */
    /*     } */
    /* } */

        /* j = 0; */
        /* for ( i =0; i<size; i++ ){ */
        /*     for ( j = 0; j<num_of_comms; j++ ){ */
        /*         strncpy(names[i*size+j], recv_buffer[i*size+j].name,32); */
        /*         /\* bytes[i] = recv_buffer[i].bytes; *\/ */
        /*         /\* printf("Copy %s",names[j]); *\/ */
        /*     } */
        /* } */
        /* /\* Clear duplicates *\/ */
        /* j = 0; */
        /* k = 0; */
        /* for ( i=0; i<size*num_of_comms; i++ ){ */
        /*     found = 0; */
        /*     for ( j=0; j<num_of_comms*size; j++ ){ */
        /*         if ( i == 0  ){ */
        /*             strncpy(snames[j], names[i],32); */
        /*             k++; */
        /*             found = 1; */
        /*             break; */
        /*         } */
        /*         else{ */
        /*             if ( strncmp(snames[j], names[i], 32) == 0 ){ */
        /*                 found = 1; */
        /*                 break; */
        /*             } */
        /*         } */
        /*     } */
        /*     if ( !found ){ */
        /*         strncpy(snames[k], names[i],32); */
        /*         k++; */
        /*     } */

        /* } */
        /* comm_num = k; */
        /* printf ( "Actual communicators = %d\n",k ); */
        /* bytes = (unsigned long long*)malloc(sizeof(unsigned long long)*comm_num); */
        /* for ( i = 0; i < comm_num; i++ ){ */
        /*     bytes[i] = 0; */
        /* } */
        /* printf("Rank 0 Communicators found:\n"); */
        /* for ( i = 0; i < comm_num; i++ ){ */
        /*     for ( j =0; j < num_of_comms*size; j++ ){ */
        /*         if ( strcmp(snames[i], names[j]) == 0 ){ */
        /*             bytes[i] = bytes[i] + recv_buffer[j].bytes; */
        /*         } */
        /*     } */
        /* } */

        /* for ( i = 0; i < comm_num; i++ ){ */
        /*     printf("Communicator %s Bytes = %llu\n",snames[i],bytes[i]); */
        /* } */
        /* free(bytes); */
    /* } */




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

