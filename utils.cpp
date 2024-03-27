#include "utils.h"
#include "symbols.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>

#define MAX_ARG_STRING_SIZE 4096
#define ENUM_TO_STRING(s) #s

const char prim_names[][NUM_OF_PRIMS]={
    ENUM_TO_STRING(Send),
    ENUM_TO_STRING(Recv),
    ENUM_TO_STRING(Isend),
    ENUM_TO_STRING(Irecv),
    ENUM_TO_STRING(Sendrecv),
    ENUM_TO_STRING(Isendrecv),
    ENUM_TO_STRING(Ssend),
    ENUM_TO_STRING(Issend),
    ENUM_TO_STRING(Rsend),
    ENUM_TO_STRING(Irsend),
    ENUM_TO_STRING(Bsend),
    ENUM_TO_STRING(Ibsend),
    ENUM_TO_STRING(Waitall),
    ENUM_TO_STRING(Wait),
    ENUM_TO_STRING(Waitany),
    ENUM_TO_STRING(Test),
    ENUM_TO_STRING(Testany),
    ENUM_TO_STRING(Testall),
    ENUM_TO_STRING(Put),
    ENUM_TO_STRING(Rput),
    ENUM_TO_STRING(Get),
    ENUM_TO_STRING(Rget),
    ENUM_TO_STRING(Accumulate),
    ENUM_TO_STRING(Raccumulate),
    ENUM_TO_STRING(Fence),
    ENUM_TO_STRING(Win_start),
    ENUM_TO_STRING(Win_complete),
    ENUM_TO_STRING(Win_post),
    ENUM_TO_STRING(Win_wait),
    ENUM_TO_STRING(Win_test),
    ENUM_TO_STRING(Bcast),
    ENUM_TO_STRING(Barrier),
    ENUM_TO_STRING(Allreduce),
    ENUM_TO_STRING(Allgather),
    ENUM_TO_STRING(Allgatherv),
    ENUM_TO_STRING(Alltoall),
    ENUM_TO_STRING(Alltoallv),
    ENUM_TO_STRING(Alltoallw),
    ENUM_TO_STRING(Reduce),
    ENUM_TO_STRING(Gather),
    ENUM_TO_STRING(Gatherv),
    ENUM_TO_STRING(Scan),
    ENUM_TO_STRING(Exscan),
    ENUM_TO_STRING(Scatter),
    ENUM_TO_STRING(Scatterv),
    ENUM_TO_STRING(Reduce_scatter),
    ENUM_TO_STRING(Reduce_scatter_block),
    ENUM_TO_STRING(Iallreduce),
    ENUM_TO_STRING(Ibcast),
    ENUM_TO_STRING(Ialltoall),
    ENUM_TO_STRING(Iscatter),
    ENUM_TO_STRING(Ibarrier),
    ENUM_TO_STRING(Iallgather),
    ENUM_TO_STRING(Iallgatherv),
    ENUM_TO_STRING(Ialltoallv),
    ENUM_TO_STRING(Ialltoallw),
    ENUM_TO_STRING(Ireduce),
    ENUM_TO_STRING(Igather),
    ENUM_TO_STRING(Igatherv),
    ENUM_TO_STRING(Iscan),
    ENUM_TO_STRING(Iexscan),
    ENUM_TO_STRING(Iscatterv),
    ENUM_TO_STRING(Ireduce_scatter),
    ENUM_TO_STRING(Ireduce_scatter_block),
    ENUM_TO_STRING(Neighbor_allgather),
    ENUM_TO_STRING(Neighbor_allgatherv),
    ENUM_TO_STRING(Neighbor_alltoall),
    ENUM_TO_STRING(Neighbor_alltoallv),
    ENUM_TO_STRING(Neighbor_alltoallw),
    ENUM_TO_STRING(Ineighbor_allgather),
    ENUM_TO_STRING(Ineighbor_allgatherv),
    ENUM_TO_STRING(Ineighbor_alltoall),
    ENUM_TO_STRING(Ineighbor_alltoallv),
    ENUM_TO_STRING(Ineighbor_alltoallw)
};

char *appname = NULL;
prof_attrs **comm_table = NULL;

void mcpt_abort (const char *fmt, ...){
  va_list args;
  va_start (args, fmt);
  fprintf (stderr, "\n\n: mpisee ABORTING ");
  vfprintf (stderr, fmt, args);
  va_end (args);
  fflush (stderr);
  PMPI_Abort(MPI_COMM_WORLD, -1);
}

char *get_appname (void){
  int pid, exelen, insize = 256;
  char *inbuf = NULL, file[256];

  pid = getpid ();
  snprintf (file, 256, "/proc/%d/exe", pid);
  inbuf = (char*) malloc (insize);
  if (inbuf == NULL){
      mcpt_abort ("unable to allocate space for full executable path.\n");
  }

  exelen = readlink (file, inbuf, 256);
  if (exelen == -1){
      if (errno != ENOENT){
          while (exelen == -1 && errno == ENAMETOOLONG){
              insize += 256;
              inbuf =(char*) realloc (inbuf, insize);
              exelen = readlink (file, inbuf, insize);
          }
          inbuf[exelen] = '\0';
          return inbuf;
      }
      else
        free (inbuf);
  }
  else{
      inbuf[exelen] = '\0';
      return inbuf;
  }
  return NULL;
}

void getProcCmdLine(int *ac, char **av) {
#ifdef __linux__
  int i = 0, pid;
  char *inbuf, file[256];
  FILE *infile;
  char *arg_ptr;

  *ac = 0;
  *av = NULL;

  pid = getpid ();
  snprintf (file, 256, "/proc/%d/cmdline", pid);
  infile = fopen (file, "r");

  if (infile != NULL){
    while (!feof (infile)){
      inbuf = (char*) malloc (MAX_ARG_STRING_SIZE);
      memset(inbuf, 0, MAX_ARG_STRING_SIZE);
      if (fread (inbuf, 1, MAX_ARG_STRING_SIZE, infile) > 0){
        arg_ptr = inbuf;
        while (*arg_ptr != 0){
          av[i] = strdup (arg_ptr);
          arg_ptr += strlen (av[i]) + 1;
          i++;
        }
      }
      free(inbuf);
    }
    *ac = i;
    fclose (infile);
  }
  else{
    mcpt_abort("Error opening file %s FILE:LINE = %d",file,__FILE__,__LINE__);
  }

#elif __APPLE__
  *ac = 0;
  *av = NULL;
#else
  mcpt_abort("Unsupported platform");
#endif
}

// void
// getRunCmd(int argc, char **argv){
//   int i;
//   ac = argc;
//   av = (char**) malloc ( sizeof(char*)*ac );
//   for ( i =0; i<ac; i++ ){
//       av[i] = strdup(argv[i]);
//   }
// }
