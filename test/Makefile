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
	mpicc -g -L$(LIBPATH) comm_topo.c -o comm_top.out -lmcpt
	mpicc -g -L$(LIBPATH) comm_split.c -o comm_split.out -lmcpt

# end
