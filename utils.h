#ifndef UTILS_H_
#define UTILS_H_

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <stdarg.h>

/* Globals */

extern char *appname;

void
mpiPi_abort (char *fmt, ...);

char *
getProcExeLink ();

#endif // UTILS_H_
