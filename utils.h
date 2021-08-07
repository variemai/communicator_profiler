#ifndef UTILS_H_
#define UTILS_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>

typedef struct profiler_attributes{
    char name[32];
    char prim[32];
    unsigned long long bytes;
    int size;
}prof_attrs;

/* Globals */
extern char *appname;
extern prof_attrs **comm_table;

void mcpt_abort (char *fmt, ...);

char * get_appname(void);

#endif // UTILS_H_
