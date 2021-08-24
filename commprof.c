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
int num_of_comms = 1;
int my_coms = 1;
/* int line_called; */
/* char file_called[32]; */

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
    int ret,rank,size;
    int i;
    prof_attrs *communicator;
    /* unsigned long long bytes; */
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    communicators =(MPI_Comm*) malloc(sizeof(MPI_Comm)*size*4);
    for ( i =0 ; i<size*4; i++ ){
        communicators[i] = NULL;
    }
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
    strcpy(communicator->parent, "NULL");
    /* communicator->index = num_of_comms; */
    /* communicator->size = size; */
    PMPI_Comm_set_attr(MPI_COMM_WORLD, namekey(), communicator);
    communicators[0] = MPI_COMM_WORLD;

    /* Table_put(table, MPI_COMM_WORLD, communicator); */
    /* num_of_comms++; */
    /* my_coms++; */
    return ret;
}

extern int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm){
    int ret;
    int size,i,flag,rank;
    prof_attrs *communicator, *com_info;
    char *name;
    ret = PMPI_Comm_create(comm, group, newcomm);
    if ( newcomm == NULL || *newcomm == MPI_COMM_NULL ){
        /* communicators[num_of_comms] = NULL; */
        num_of_comms++;
        return ret;
    }
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    //append mother's name if mother != WORLD
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
    /*             sprintf(communicator->name, "%s_c%d",com_info->name,my_coms); */
    /*             printf("New comm with name %s\n",communicator->name); */
    /*         } */
    /*         else{ */
    /*             mcpt_abort("Flag in file:%s line:%d invalid\nAborting\n",__FILE__,__LINE__); */
    /*         } */
    /*     } */

    /* } */
    /* else */
    sprintf(communicator->name,"c%d",my_coms);
    communicator->bytes = 0;
    strcpy(communicator->prim,"COMM_CREATE");
    /* /\* communicator->name = buf; *\/ */
    /* /\* communicator->prim = prim; *\/ */
    /* MPI_Comm_size(*newcomm, &size); */
    /* communicator->size = size; */
    /* /\* MPI_Comm_set_name(*newcomm, buf); *\/ */
    /* printf("New Comm Entry %s\n",communicator->name); */
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[my_coms] = *newcomm;
    /* comm_table[num_of_comms]=communicator; */
    /* Table_put(table, *newcomm, communicator); */
    /* MPI_Comm_get_name(*newcomm, buffer, &len); */
    /* if ( len && buffer){ */
    /*     printf ( "Name = %.*s\n",len,buffer ); */
    /* } */
    num_of_comms++;
    my_coms++;
    /* printf("%p",__builtin_return_address(0)); */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
    /* printf("RANK %d: Num of Comms After create = %d\n",rank,num_of_comms); */
    return ret;
}

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm){
    int ret,i,flag;
    prof_attrs *communicator,*com_info;
    ret = PMPI_Comm_split(comm, color, key, newcomm);
    communicator = (prof_attrs*) malloc (sizeof(prof_attrs));
    sprintf(communicator->name,"s_%d.%d",my_coms,color);
    //append mother's name if mother != WORLD
    if ( comm != MPI_COMM_WORLD ){
    /*     //find the mother */
        for ( i = 0; i<num_of_comms; i++ ){
            if ( comm  == communicators[i] )
                break;
        }
        if ( i == num_of_comms ){
            // Should have found the mother
            fprintf(stderr, "Error: could not find the mother of communicator.\nAborting\n");

        }
        else{
            PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag);
            if ( flag ){
                strcpy(communicator->parent, com_info->name);
            }
            else{
                /* fprintf(stderr, "Flag invalid\nAborting\n"); */
                mcpt_abort("Flag in file:%s line:%d invalid\nAborting\n",__FILE__,__LINE__);
            }
        }
    }
    else{
        strcpy(communicator->parent, "WORLD");
    }
    /* printf("New comm with name %s and parent %s\n",communicator->name,communicator->parent); */
    communicator->bytes = 0;
    strcpy(communicator->prim,"COMM_SPLIT");
    PMPI_Comm_set_attr(*newcomm, namekey(), communicator);
    communicators[my_coms] = *newcomm;
    num_of_comms+=color+1;
    my_coms++;
    /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
    /* printf("RANK %d: Num of Comms After split = %d\n",rank,num_of_comms); */
    return ret;
}


extern int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest,
                    int tag, MPI_Comm comm){
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
    PMPI_Type_size(datatype, &size);
    if ( flag ){
        /* printf("MPI_Send to communicator %s\n",communicator->name); */
        sum = communicator->bytes;
        sum = sum + count * size;
        communicator->bytes = sum;
    /*     /\* Table_put(table, comm, communicator); *\/ */
    }
    /* communicator = Table_get(table,comm); */
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
    prof_attrs *array;
    int rank,size;
    int i,j,k,flag,found;
    prof_attrs *com_info;
    prof_attrs *recv_buffer;
    prof_attrs dummy;
    char **names, **parents, **unames, **uparents;
    /* int found, comm_num; */
    int *total_comms,total;
    /* unsigned long long *bytes; */
    /* char **names, **names_buf; */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* free(communicators); */

    /* Before we do anything we should first learn the number of communicators */
    total_comms = (int*)malloc(sizeof(int)*size*num_of_comms);
    PMPI_Gather(&num_of_comms, 1, MPI_INT, total_comms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* PMPI_Reduce(&num_of_comms, total_comms, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD); */
    if ( rank == 0 ){
        for ( i = 0; i< size; i++ ){
            if ( total_comms[i] > num_of_comms ){
                num_of_comms = total_comms[i];
            }
        }
    }
    PMPI_Bcast(&num_of_comms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* names = (char**)malloc(sizeof(char*)*num_of_comms); */
    /* for ( i =0; i<num_of_comms*size; i++ ){ */
    /*     snames[i] = (char*)malloc(64); */
    /* } */
    /* j = 0; */
    /* for ( i =0; i<num_of_comms; i++ ){ */
    /*     if ( communicators[i] != NULL ){ */
    /*         PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag); */
    /*         if ( flag ){ */
    /*             names[j]=strndup(com_info->name,64); */
    /*             j++; */
    /*         } */
    /*     } */
    /*     else{ */
    /*         names[i]=strdup(names[j-1]); */
    /*     } */
    /* } */
    /* names[i]="\0"; */
    /* for ( i =0; i<num_of_comms; i++ ){ */
    /*     printf("Rank %d %s\n",rank,names[i]); */
    /* } */
    /* PMPI_Gather(names, num_of_comms*64, MPI_CHAR, snames, num_of_comms*64, MPI_CHAR, 0, MPI_COMM_WORLD); */

    /* if ( rank == 0 ){ */
    /*     for ( i =0; i<num_of_comms*size; i++ ){ */
    /*         printf("%s\n",snames[i]); */
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


    array =(prof_attrs*) malloc(sizeof(prof_attrs)*num_of_comms);
    recv_buffer = (prof_attrs*) malloc (sizeof(prof_attrs)*num_of_comms*size);
    MPI_Datatype types[5] = { MPI_CHAR, MPI_CHAR, MPI_CHAR,MPI_UNSIGNED_LONG_LONG, MPI_INT };
    int blocklen[5] = {NAMELEN,NAMELEN,PRIMLEN,1,1};
    MPI_Aint displacements[5];
    MPI_Aint base_address;
    MPI_Datatype profiler_data;
    MPI_Get_address(&dummy, &base_address);
    MPI_Get_address(&dummy.name[0], &displacements[0]);
    MPI_Get_address(&dummy.parent[0], &displacements[1]);
    MPI_Get_address(&dummy.prim[0], &displacements[2]);
    MPI_Get_address(&dummy.bytes, &displacements[3]);
    MPI_Get_address(&dummy.size, &displacements[4]);
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    displacements[3] = MPI_Aint_diff(displacements[3], base_address);
    displacements[4] = MPI_Aint_diff(displacements[4], base_address);

    MPI_Type_create_struct(5, blocklen, displacements, types, &profiler_data);
    MPI_Type_commit(&profiler_data);
    for ( i = 0; i < num_of_comms; i++ ){
        if ( communicators[i] != NULL ){
            PMPI_Comm_get_attr(communicators[i], namekey(), &com_info, &flag);
            if ( flag ){
                strcpy(array[i].name, com_info->name);
                strcpy(array[i].parent, com_info->parent);
                strcpy(array[i].prim, com_info->prim);
                array[i].bytes = com_info->bytes;
                array[i].size = com_info->size;
            }
        }
        else{
            strcpy(array[i].name, "NULL");
            strcpy(array[i].parent, "NULL");
            strcpy(array[i].prim, "NULL");
            array[i].bytes = 0;
            array[i].size = 0;
        }
    }
    PMPI_Gather(array, num_of_comms*sizeof(prof_attrs), MPI_BYTE, recv_buffer,
                num_of_comms*sizeof(prof_attrs), MPI_BYTE, 0, MPI_COMM_WORLD);

    if ( rank == 0 ){
        names = ( char**)malloc(sizeof(char*)*num_of_comms*size);
        parents = ( char**)malloc(sizeof(char*)*num_of_comms*size);
        unames = (char **) malloc (sizeof(char*)*num_of_comms*size);
        uparents =(char **) malloc (sizeof(char*)*num_of_comms*size);
        j = 0;
        for ( i =0; i<size*num_of_comms; i++ ){
            if ( strcmp(recv_buffer[i].name, "NULL") != 0 ){
                names[j]=strdup(recv_buffer[i].name);
                parents[j]=strdup(recv_buffer[i].parent);
                unames[j] = strdup("NULL");
                uparents[j] = strdup("NULL");
                j++;
            }
        }
        total = j;
        num_of_comms = 1;
        j = 0;
        printf("TOTAL %d\n",total);
        for ( i=0; i<total; i++ ){
                if ( strcmp(names[i], "WORLD") != 0 ){
                    found = 0;
                    for ( k =0; k<total; k++ ){
                        if ( strcmp(names[i], unames[k] ) == 0 &&
                             strcmp(parents[i], uparents[k] ) == 0 ){
                            found = 1;
                        }
                    }
                    if ( !found ){
                        strcpy(unames[j], names[i]);
                        strcpy(uparents[j], parents[i]);
                        j++;
                        num_of_comms++;
                    }
                }
        }
        printf( "Num of comms = %d\n",num_of_comms);
        for ( i =0; i<num_of_comms-1; i++ )
            printf("Comm: %s Parent: %s\n",unames[i],uparents[i]);
        for ( i =0; i<total; i++ ){
            free(names[i]);
            free(unames[i]);
            free(parents[i]);
            free(uparents[i]);
        }
    }

    /* fclose(fp); */
    /* if (communicators) */
    /*     FREE(communicators); */
    return PMPI_Finalize();
}

