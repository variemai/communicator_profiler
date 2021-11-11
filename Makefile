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

#CC=mpicc
CFLAGS=-O3 -flto=auto -Wall -std=gnu99 -pedantic -fPIC -I$(LIBPATH)
SHARED=-shared

#all:
#	@echo PATH IS  $(LIBPATH)

all: default $(DSTRUCTF) tests
	@echo MCPT compiled

default: libciface
#	@echo PATH IS $(LIBPATH)
	$(CC) $(CFLAGS) $(SHARED) $(SRCS) -L$(LIBPATH) -lciface -o $(LIBF)/libmcpt.so

libciface: $(DSTRUCTF)
	$(MAKE) -C $<

tests: test
	$(MAKE) -C $<

#$(OBJS): $(SRCS)
#	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o *.out $(LIBF)/* $(DSTRUCTF)/*.o $(DSTRUCTF)/*.so test/*.out

# end
