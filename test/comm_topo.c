/* virt_topo.c
* https://computing.llnl.gov/tutorials/mpi/#Virtual_Topologies
*/
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#define SIZE 16
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int main (int argc, char *argv[]){
    int rank = 0, source, dest, outbuf, i;
    MPI_Comm cartcomm;
    int world_rank, world_size;
    int periods[2] = {1,1};
    int dims[2] = {4,4};
    MPI_Request reqs[8];
    MPI_Status stats[8];
    int inbuf[4]={MPI_PROC_NULL,MPI_PROC_NULL, MPI_PROC_NULL,MPI_PROC_NULL},
        nbrs[4], reorder=0, coords[2];
    int tag = 1;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if ( world_size == SIZE ){
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods,reorder, &cartcomm);
        MPI_Comm_rank(cartcomm, &rank);
        MPI_Cart_coords(cartcomm, rank, 2, coords);
        MPI_Cart_shift(cartcomm, 0, 1, &nbrs[UP], &nbrs[DOWN]);
        MPI_Cart_shift(cartcomm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);

    /* printf("P[%d]: coords= %d %d \n", */
    /*        rank,coords[0],coords[1]); */
    /* printf("P[%d]: neighbors(u,d,l,r)=%d %d %d %d\n", */
       /* rank,nbrs[UP],nbrs[DOWN],nbrs[LEFT], nbrs[RIGHT]); */
    outbuf = rank;
    for (i=0; i<4; i++) {
        dest = nbrs[i];
        source = nbrs[i];
        /* MPI_Sendrecv(&outbuf, 1, MPI_INT, dest, tag, &inbuf[i], 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
        MPI_Isend(&outbuf, 1, MPI_INT, dest, tag,
                  cartcomm, &reqs[i]);
        MPI_Irecv(&inbuf[i], 1, MPI_INT, source, tag,
                  cartcomm, &reqs[i+4]);
    }
    MPI_Waitall(8, reqs, stats);
    /* printf("rank= %d inbuf(u,d,l,r)= %d %d %d %d\n", */
    /*        rank,inbuf[UP],inbuf[DOWN],inbuf[LEFT],inbuf[RIGHT]); */
    }
    else{
        printf("Must specify %d processors. Terminating.\n",SIZE);
    }

    MPI_Finalize();
}

