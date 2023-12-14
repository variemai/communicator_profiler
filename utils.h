#ifndef UTILS_H_
#define UTILS_H_

#include <cstdint>
#include <mpi.h>
#include <stdint.h>
#include <unistd.h>
#define NAMELEN 256
#define PRIMLEN 16
#define NUM_OF_PRIMS 36
#define MAX_ARG_STRING_SIZE 4096
#define MAX_ARGS 1024
#define MAX_ARG_SIZE 64
#define NUM_BUCKETS 7
// Buckets for message sizes
const uint32_t buckets[NUM_BUCKETS-1] = {7,10,13,16,20,25};

enum primitives{
Send,           /* DO NOT ADD ANYHTING BEFORE THIS */
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
Ibarrier   /* DO NOT ADD ANYTHING AFTER THIS */
};

typedef struct communicator_attributes_msgsizes{
    char name[NAMELEN];
    int size;
    uint32_t prims[NUM_OF_PRIMS];
    uint64_t prim_bytes[NUM_OF_PRIMS];
    double time_info[NUM_OF_PRIMS];
    double buckets_time[NUM_OF_PRIMS][NUM_BUCKETS];
    uint64_t buckets_msgs[NUM_OF_PRIMS][NUM_BUCKETS];
}comm_attrs_msgsizes;

typedef struct profiler_attributes{
    char name[NAMELEN];
    int size;
    uint64_t bytes;
    uint64_t msgs;
    uint32_t prims[NUM_OF_PRIMS];
    uint64_t prim_bytes[NUM_OF_PRIMS];
    double time_info[NUM_OF_PRIMS];
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
