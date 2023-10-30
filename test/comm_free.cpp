#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "../utils.h"
#include "../commprof.h"
#include <assert.h>

int main (int argc, char *argv[]){
    int flag,rank,i;
    MPI_Comm comm0,comm1,comm2;
    MPI_Group commgrp;
    prof_attrs *communicator0 = NULL;
    prof_attrs *communicator1 = NULL;
    prof_attrs *communicator2 = NULL;
    int buffer[2];
    unsigned long long bytes[3] = {0,0,0};
    flag = 0;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    MPI_Comm_group(MPI_COMM_WORLD,&commgrp);
    MPI_Comm_create(MPI_COMM_WORLD, commgrp, &comm0);
    MPI_Comm_create(MPI_COMM_WORLD, commgrp, &comm1);
    MPI_Comm_create(MPI_COMM_WORLD, commgrp, &comm2);
    buffer[0] = rank;
    buffer[1] = rank;
    MPI_Bcast(&buffer, 2, MPI_INT, 0, comm0);
    MPI_Bcast(&buffer, 1, MPI_INT, 0, comm1);
    MPI_Bcast(&buffer, 2, MPI_INT, 0, comm1);

    PMPI_Comm_get_attr(comm0, namekey(), &communicator0, &flag);
    if ( flag )
        bytes[0] = communicator0->prim_bytes[Bcast];

    MPI_Comm_free(&comm0);

    PMPI_Comm_get_attr(comm1, namekey(), &communicator1, &flag);
    if ( flag )
        bytes[1] = communicator1->prim_bytes[Bcast];
    MPI_Comm_free(&comm1);


    PMPI_Comm_get_attr(comm2, namekey(), &communicator2, &flag);
    if ( flag )
        bytes[2] = communicator2->prim_bytes[Bcast];
    MPI_Comm_free(&comm2);

    printf("%d\n",local_cid);
    for ( i = 1; i<local_cid; i++ ){
        assert(bytes[i-1] == local_comms[i]->prim_bytes[Bcast]);
    }
    MPI_Finalize();
    return 0;
}
