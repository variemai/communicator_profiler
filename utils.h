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
#define NAMELEN 512
#define PRIMLEN 16

#ifdef USE_MPI3_CONSTS
typedef const void mpip_const_void_t;
typedef const int mpip_const_int_t;
typedef const char mpip_const_char_t;
typedef const MPI_Datatype mpip_const_datatype_t;
#else
typedef void mpip_const_void_t;
typedef int mpip_const_int_t;
typedef char mpip_const_char_t;
#endif

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
