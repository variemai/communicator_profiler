#ifndef UTILS_H_
#define UTILS_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>
#include "table.h"

typedef struct _communicator{
    char *name;
    int64_t *bytes;
    MPI_Comm comm;
    char *prim;
}commtor;

/* Globals */
extern char *appname;
extern int num_of_comms;
extern Table_T table;
extern commtor **comm_table;

void mcpt_abort (char *fmt, ...);

char * get_appname(void);

#endif // UTILS_H_
