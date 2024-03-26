#ifndef COMMPROF_H_
#define COMMPROF_H_

#include "mpi.h"
#include <unordered_map>
#define MAX_DIMS 8
#include "utils.h"
#include <vector>
extern int prof_enabled;
extern prof_attrs **local_data;
// extern prof_attrs **local_comms;
extern std::vector<prof_attrs*> local_communicators;
extern std::unordered_map<MPI_Request, MPI_Comm> requests_map;
extern int local_cid;
extern int my_coms;

extern "C" {
int namekey();
}
extern "C" {
int namedel(MPI_Comm comm, int keyval, void *attr, void *s);
}

extern "C" {
prof_attrs*
profile_this(MPI_Comm comm, int64_t count,MPI_Datatype datatype,int prim,
             double t_elapsed,int root);
}
#endif
