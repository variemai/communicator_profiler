EXECS=exec
CC = gcc -std=gnu99 -ansi -pedantic
CFLAGS = -g -Wall -fno-inline -O3
SRCS = arith.c except.c assert.c memchk.c stack.c queue.c list.c atom.c table.c
OBJS = $(SRCS:.c=.o)
MAIN=exec
LIBFLAGS= -shared -ldl
AR=ar rcs

.PHONY: depend clean

all: SHLIB
	@echo  Library has been compiled

static: LIB
	@echo  Shared Library has been compiled

SHLIB: $(OBJS)
	$(CC) $(CFLAGS) $(LIBFLAGS) $(OBJS) -o libcintfaces.so

$(MAIN): $(OBJS) main.c
	$(CC) $(CFLAGS) main.c -o $(MAIN) $(OBJS)

LIB: $(OBJS)
	$(AR) libmyclib.a $(OBJS)

.c.o:
	$(CC) -fPIC $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN) *.a *.so

depend: $(SRCS)
	makedepend $(INCLUDES) $^

#DO NOT DELETE
