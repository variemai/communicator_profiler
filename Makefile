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
DSTRUCTF=datastructlib
LIBPATH=$(abspath $(DSTRUCTF))

CC=mpicc
CFLAGS=-O3 -Wall -pedantic -fPIC -I$(LIBPATH)
SHARED=-shared

#all:
#	@echo PATH IS  $(LIBPATH)

all: default $(DSTRUCTF)
	@echo MCPT compiled

default: libciface
	$(CC) $(CFLAGS) $(SHARED) $(SRCS) -L$(LIBPATH) -lciface -o $(LIBF)/libmcpt.so

libciface: $(DSTRUCTF)
	$(MAKE) -C $<

#$(OBJS): $(SRCS)
#	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o *.out $(LIBF)/* $(DSTRUCTF)/*.o $(DSTRUCTF)/*.so

# end
