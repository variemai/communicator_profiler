#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>

int main (int argc, char *argv[]){
    int size;
    MPI_Comm newcomm,splitcomm,subcomm;
    MPI_Group group;
    int rank;
    char *buffer;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* printf("Hello from all %d ranks\n",size); */
    /* MPI_Comm_get_name(MPI_COMM_WORLD, comm_name, &len); */
    MPI_Comm_group(MPI_COMM_WORLD, &group);
    /* wrapper_shit(MPI_COMM_WORLD, group, &newcomm); */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm) */
    /* if ( rank == 0 || rank == 1 ){ */
    MPI_Comm_create(MPI_COMM_WORLD, group,&newcomm);
    MPI_Comm_rank(newcomm, &rank);
    /* } */
    /* MPI_Comm_set_name(newcomm, "BOURDELO"); */
    /* MPI_Comm_get_name(newcomm, comm_name, &len); */
    /* set_attr_wrapper(newcomm, 0, val); */
    /* MPI_Comm_set_attr(newcomm, namekey(), com_info); */
    /* comm_create_wrapper(MPI_COMM_WORLD, group, &newcomm); */
    /* bullshit(&newcomm) ; */
    /* sprintf(comm_name, "TRIBOURDELO"); */
    /* MPI_Comm_set_attr(newcomm, 0, comm_name); */
    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm);
    buffer = malloc(64);
    if ( rank  % 2 == 0  ){
        MPI_Comm_group(splitcomm, &group);
        MPI_Comm_create(splitcomm, group, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        /* MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); */
        /* printf("Rank = %d, World = %d\n",rank,world_rank); */
        if ( rank % 2 == 0 ){
            MPI_Send(&size,1,MPI_INT,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(&size, 1, MPI_INT, rank-1, 0, subcomm, MPI_STATUSES_IGNORE);
        }
        /* MPI_Recv(buffer, 64, MPI_CHAR, 1, 1, newcomm, MPI_STATUS_IGNORE); */
        /* MPI_Send(&size,1,MPI_INT,1,0,splitcomm); */
        /* printf("I received %s from rank %d\n",buffer,rank); */
        /* MPI_Send(&size,1,MPI_INT,1,2,MPI_COMM_WORLD); */
        /* send_wrapper(&size, 1, MPI_INT, 1, 0, newcomm); */
        /* MPI_Comm_get_attr(newcomm, namekey(), &com_info, &flag); */
        /* if ( flag ){ */
        /*     printf("Rank %d: %s HAS %llu\n",rank,com_info->name,com_info->bytes); */
        /* } */
    }
    else{

    }
    /* if ( rank == 1 ){ */
    /*     MPI_Recv(&size, 1, MPI_INT, 0, 0, newcomm, MPI_STATUS_IGNORE); */
    /*     MPI_Send(buffer,64,MPI_CHAR,0,1,newcomm); */
    /* } */
    /* if ( rank == 2 ){ */
    /*     MPI_Send(buffer,64,MPI_CHAR,3,3,MPI_COMM_WORLD); */
    /*     MPI_Recv(&size, 1, MPI_INT, 0, 0, splitcomm, MPI_STATUS_IGNORE); */
    /*     printf("I received size = %d from rank 0\n",size); */
    /* } */
    /* if ( rank == 3 ){ */
    /*     MPI_Recv(buffer, 64, MPI_CHAR, 2, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
    /* } */
    MPI_Finalize();
}
