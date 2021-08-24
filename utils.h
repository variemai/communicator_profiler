#ifndef UTILS_H_
#define UTILS_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>
#define NAMELEN 64
#define PRIMLEN 16

typedef struct profiler_attributes{
    char name[NAMELEN];
    char parent[NAMELEN];
    char prim[PRIMLEN];
    unsigned long long bytes;
    int size;
}prof_attrs;

/* Globals */
extern char *appname;
extern prof_attrs **comm_table;

void mcpt_abort (char *fmt, ...);

char * get_appname(void);

#endif // UTILS_H_
