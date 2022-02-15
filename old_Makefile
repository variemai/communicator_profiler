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

CFLAGS=-O3 -flto -Wall -std=gnu99 -pedantic -fPIC -I$(LIBPATH)
## Use the line to run withiout lto flag if your compiler does not support it
#CFLAGS= -O3 -Wall -march=native -std=gnu99 -pedantic -fPIC -I$(LIBPATH)
SHARED=-shared
export CFLAGS
export CC

#all:
#	@echo PATH IS  $(LIBPATH)

all: default $(DSTRUCTF)
	@echo mpisee compiled

default: libciface
#	@echo PATH IS $(LIBPATH)
	mkdir -p $(LIBF)
	$(CC) $(CFLAGS) $(SHARED) $(SRCS) -L$(LIBPATH) -l$(DSTRUCTF) -o $(LIBF)/libmpisee.so

libciface: $(DSTRUCTF)
	$(MAKE) -C $<

tests: test default
	$(MAKE) -C $<

#$(OBJS): $(SRCS)
#	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o *.out $(LIBF)/* $(DSTRUCTF)/*.o $(DSTRUCTF)/*.so test/*.out

# end
