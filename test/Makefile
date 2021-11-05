##
# Testfile
#
# @file
# @version 0.1

LIBF=../lib
LIBPATH=$(abspath $(LIBF))

all:
	mpicc -g -L$(LIBPATH) hello_mpi.c -lmcpt
	mpicc -g -L$(LIBPATH) comm_group.c -o comm_group.out -lmcpt
	mpicc -g -L$(LIBPATH) comm_topo.c -o comm_topo.out -lmcpt
	mpicc -g -L$(LIBPATH) comm_split.c -o comm_split.out -lmcpt
	mpicc -g -L$(LIBPATH) comm_split_v2.c -o comm_split_v2.out -lmcpt
	mpicc -g -L$(LIBPATH) alltoallv.c -o alltoallv.out  -lmcpt
	mpicc -g -L$(LIBPATH) cart_sub.c -o cart_sub.out  -lmcpt
	mpicc -g -L$(LIBPATH) tuwqsort.c -o tuwqsort.out  -lmcpt
	mpicc -g -L$(LIBPATH) wait.c -o wait.out  -lmcpt
	mpicc -g -L$(LIBPATH) waitall.c -o waitall.out  -lmcpt
	mpicc -g -L$(LIBPATH) sendrecv.c -o sendrecv.out  -lmcpt

# end
