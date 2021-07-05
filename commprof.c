#include "utils.h"
#include <mpi.h>

extern int MPI_Init(int *argc, char ***argv){
    int ret,rank;
    /* Call Init before profiler initializations */
    ret = PMPI_Init(argc,argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( rank == 0 ){
        appname = (malloc(sizeof(char)*256));
        appname = getProcExeLink();
        printf("MPI Communicator Profiling Tool\nProfiling application %s\n",appname);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return ret;
}
