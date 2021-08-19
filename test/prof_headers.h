#ifndef PROF_HEADERS_H_
#define PROF_HEADERS_H_

extern int line_called;

#define MPI_Comm_split(comm,color,key,newcomm)                          \
    do { line_called = __LINE__;                                        \
        MPI_Comm_split(comm,color,key,newcomm);  }                       \
    while (0)

#endif // PROF_HEADERS_H_
