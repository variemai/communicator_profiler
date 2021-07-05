#include "utils.h"

char *appname = NULL;

void
mcpt_abort (char *fmt, ...)
{
  va_list args;
  va_start (args, fmt);
  fprintf (stderr, "\n\n: MCPT ABORTING ");
  vfprintf (stderr, fmt, args);
  va_end (args);
  fflush (stderr);
  PMPI_Abort(MPI_COMM_WORLD, -1);
}

char *
getProcExeLink ()
{
  int pid, exelen, insize = 256;
  char *inbuf = NULL, file[256];

  pid = getpid ();
  snprintf (file, 256, "/proc/%d/exe", pid);
  inbuf = malloc (insize);
  if (inbuf == NULL)
    {
      mcpt_abort ("unable to allocate space for full executable path.\n");
    }

  exelen = readlink (file, inbuf, 256);
  /* printf("EXELEN %d\n",exelen); */
  if (exelen == -1)
    {
      if (errno != ENOENT)
        {
          while (exelen == -1 && errno == ENAMETOOLONG)
            {
              insize += 256;
              inbuf = realloc (inbuf, insize);
              exelen = readlink (file, inbuf, insize);
            }
          inbuf[exelen] = '\0';
          return inbuf;
        }
      else
        free (inbuf);
    }
  else
    {
      inbuf[exelen] = '\0';
      return inbuf;
    }
  return NULL;
}
