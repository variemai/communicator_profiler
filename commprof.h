#ifndef COMMPROF_H_
#define COMMPROF_H_

#define MAX_DIMS 8
#include "utils.h"

extern int prof_enabled;
extern prof_attrs **local_data;
extern prof_attrs **local_comms;
extern int local_cid;
extern int my_coms;

extern "C" {
int namekey();
}
extern "C" {
int namedel(MPI_Comm comm, int keyval, void *attr, void *s);
}
#endif
