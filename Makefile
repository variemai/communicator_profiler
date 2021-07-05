##
# MPI Communicator profiler
#
# @file
# @version 0.1

srcdir=.

.PHONY: default all clean

SRCS = commprof.c utils.c

OBJS=$(SRCS:.c=.o)
LIBF=./lib

CC=mpicc
CFLAGS=-O3 -Wall -pedantic -fPIC
SHARED=-shared

all:default

default: $(SRCS)
	$(CC) $(CFLAGS) $(SHARED) $^ -o $(LIBF)/libmcpt.so

#$(OBJS): $(SRCS)
#	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o *.out $(LIBF)/*

# end
