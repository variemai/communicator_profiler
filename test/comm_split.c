#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#define NUMPROCS 8

int main (int argc, char *argv[]){
    int size;
    MPI_Comm splitcomm,subcomm;
    int rank,color;
    char *buffer;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( size != NUMPROCS ){
        if ( rank == 0 ){
            fprintf(stderr, "Error: Application runs only with %d processes\n",NUMPROCS);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank / 2, &splitcomm);
    buffer = malloc(64);
    if ( rank  % 2 == 0  ){
        MPI_Send(buffer, 32, MPI_BYTE, rank+1, 0, MPI_COMM_WORLD);
        MPI_Comm_rank(splitcomm, &rank);
        color = rank %2;
        MPI_Comm_split(splitcomm, color, rank / 2, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank < 1 ){
            MPI_Send(buffer,64,MPI_BYTE,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(buffer, 64, MPI_BYTE, rank-1, 0, subcomm, MPI_STATUS_IGNORE);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            printf("Process %d recvd %d bytes in subcomm color: %d\n",rank,64,color);
        }

    }
    else{
        MPI_Recv(buffer, 32, MPI_BYTE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Comm_rank(splitcomm, &rank);
        color = rank %2;
        MPI_Comm_split(splitcomm, color, rank / 2, &subcomm);
        MPI_Comm_rank(subcomm, &rank);
        if ( rank % 2 == 0 ){
            MPI_Send(buffer,16,MPI_BYTE,rank+1,0,subcomm);
        }
        else{
            MPI_Recv(buffer, 16, MPI_BYTE, rank-1, 0,subcomm, MPI_STATUS_IGNORE);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            printf("Process %d recvd %d bytes in subcomm color: %d\n",rank,16,color);
        }
    }
    free(buffer);
    MPI_Finalize();
}

