#ifndef UTILS_H_
#define UTILS_H_

#include <cstdint>
#include <mpi.h>
#include <unistd.h>
#include <string>
#define NAMELEN 256
#define NUM_OF_PRIMS 58
#define MAX_ARG_STRING_SIZE 4096
#define MAX_ARGS 1024
#define MAX_ARG_SIZE 64
#define NUM_BUCKETS 7
// Buckets for message sizes
const uint32_t buckets[NUM_BUCKETS-1] = {7,10,13,16,20,25};
// This is used to index the prim_names array
// Warning: This has to be in the same order as the prim_names array
enum primitives{
Send,
Recv,
Isend,
Irecv,
Ssend,
Issend,
Probe,
Iprobe,
Waitall,
Wait,
Waitany,
Test,
Testany,
Testall,
Sendrecv,
Bcast,
Barrier,
Allreduce,
Allgather,
Allgatherv,
Alltoall,
Alltoallv,
Alltoallw,
Reduce,
Gather,
Gatherv,
Scan,
Exscan,
Scatter,
Scatterv,
Reduce_scatter,
Iallreduce,
Ibcast,
Ialltoall,
Iscatter,
Ibarrier,
Iallgather,
Iallgatherv,
Ialltoallv,
Ialltoallw,
Ireduce,
Igather,
Igatherv,
Iscan,
Iexscan,
Iscatterv,
Ireduce_scatter,
Ireduce_scatter_block,
Neighbor_allgather,
Neighbor_allgatherv,
Neighbor_alltoall,
Neighbor_alltoallv,
Neighbor_alltoallw,
Neighbor_iallgather,
Neighbor_iallgatherv,
Neighbor_ialltoall,
Neighbor_ialltoallv,
Neighbor_ialltoallw
};

struct CommData {
    std::string name;
    int size;
};

struct DataEntry {
    int rank;
    int commId;
    int operationId;
    int bufferSizeMax;
    int bufferSizeMin;
    int calls;
    double time;
};

typedef struct profiler_attributes{
    char name[NAMELEN];
    int size;
    double buckets_time[NUM_OF_PRIMS][NUM_BUCKETS];
    uint64_t buckets_msgs[NUM_OF_PRIMS][NUM_BUCKETS];
}prof_attrs;

typedef struct request_data{
    MPI_Request *req;
    MPI_Comm comm;
}rq;

/* Globals */
extern char *appname;
extern const char prim_names[][NUM_OF_PRIMS];
extern int ac;
extern char *av[MAX_ARGS];

void mcpt_abort (const char *fmt, ...);

void getProcCmdLine (int *ac, char **av);
char * get_appname(void);

void getRunCmd(int argc, char **argv);

#endif // UTILS_H_
