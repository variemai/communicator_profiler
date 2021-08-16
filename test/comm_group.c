#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>

int main (int argc, char *argv[]){
    int size;
    MPI_Comm newcomm,splitcomm,subcomm;
    MPI_Comm comm;
    MPI_Group group, newgroup, world_group;
    int rank, world_rank, world_size;
    char *buffer;
    const int ranks[2] = {0,1};
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /* printf("Hello from all %d ranks\n",size); */
    /* MPI_Comm_get_name(MPI_COMM_WORLD, comm_name, &len); */
    /* wrapper_shit(MPI_COMM_WORLD, group, &newcomm); */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* if ( rank == 0 || rank == 1 ){ */
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
        /* MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); */
        /* MPI_Comm_rank(splitcomm, &rank); */
        /* printf("Rank = %d, World = %d\n",rank,world_rank); */
        MPI_Comm_group(splitcomm, &group);
        MPI_Comm_create(splitcomm, group, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank % 2 == 0 ){
            MPI_Send(&size,1,MPI_INT,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(&size, 1, MPI_INT, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
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
        MPI_Comm_group(splitcomm, &group);
        MPI_Comm_create(splitcomm, group, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        /* MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); */
        /* printf("Rank = %d, World = %d\n",rank,world_rank); */
        if ( rank % 2 == 0 ){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
        }
    }
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if ( MPI_Group_excl(world_group, 2, ranks, &newgroup) != MPI_SUCCESS ){
        fprintf(stderr,"group_include_failed");
        exit(EXIT_FAILURE);
    }
    if ( MPI_Group_incl(world_group, 2, ranks, &group) != MPI_SUCCESS ){
        fprintf(stderr,"group_include_failed");
        exit(EXIT_FAILURE);
    }
    if ( newgroup == MPI_GROUP_NULL ){
            fprintf(stderr,"NULL NEWGROUP\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }
    if ( group == MPI_GROUP_NULL ){
            fprintf(stderr,"NULL GROUP\n");
            MPI_Finalize();
            exit(EXIT_FAILURE);
    }
    // id = 0
    if ( MPI_Comm_create(MPI_COMM_WORLD, newgroup, &newcomm) != MPI_SUCCESS ){
        fprintf(stderr,"comm_create_failed");
        exit(EXIT_FAILURE);
    }
    // id = 1
    if ( MPI_Comm_create(MPI_COMM_WORLD, group, &comm) != MPI_SUCCESS ){

        fprintf(stderr,"comm_create_failed");
        exit(EXIT_FAILURE);
    }
    if ( newcomm != MPI_COMM_NULL && newcomm != NULL){
            /* fprintf(stderr,"NULL COMMUNICATOR\n"); */
            /* MPI_Finalize(); */
            /* exit(EXIT_FAILURE); */
        MPI_Comm_rank(newcomm,&rank);
        MPI_Comm_size(newcomm,&size);
        if ( rank == 0 ){
            MPI_Send(&size,1,MPI_INT,rank+1,0,newcomm);
        }
        if ( rank == 1 ){
            MPI_Recv(&size, 1, MPI_INT, rank-1, 0, newcomm, MPI_STATUS_IGNORE);
        }
        printf("Communicator X: World rank = %d, newcomm rank = %d, newcomm size = %d\n",world_rank,rank,size);
        // id = 0_0.0
        // id = 0_0.1
        /* MPI_Comm_split(newcomm, rank % 2, rank / 2, &splitcomm); */
        // id = 0_1
        /* MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) */
    }
    if ( comm != MPI_COMM_NULL && comm != NULL){
        /* fprintf(stderr,"NULL COMMUNICATOR\n"); */
        /* MPI_Finalize(); */
        /* exit(EXIT_FAILURE); */
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&size);
        printf("Communicator Y: World rank = %d, comm rank = %d, comm size = %d\n",world_rank,rank,size);
        if ( rank == 0 ){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,comm);
        }
        if ( rank == 1 ){
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, comm, MPI_STATUS_IGNORE);
        }
        // id = 1_0.0
        // id = 1_0.0
        /* MPI_Comm_split(comm, rank % 2, rank / 2, &splitcomm); */
        // id = 1_1
        /* MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) */
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
