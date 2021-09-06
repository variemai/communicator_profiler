#ifndef UTILS_H_
#define UTILS_H_

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#define NAMELEN 128
#define PRIMLEN 16

typedef struct profiler_attributes{
    char name[NAMELEN];
    /* char parent[NAMELEN]; */
    /* char prim[PRIMLEN]; */
    uint64_t bytes;
    uint32_t msgs;
    int size;
}prof_attrs;

/* Globals */
extern char *appname;

void mcpt_abort (char *fmt, ...);

void getProcCmdLine (int *ac, char **av);
char * get_appname(void);

#endif // UTILS_H_
