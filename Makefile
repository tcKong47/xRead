CC = gcc
CFLAGS = -g -Wall -O2 -Wno-unused-function
# -Wc++-compat
OBJS = main.o ktime.o bseq.o thread.o bit_operation.o overlapping.o index.o graph.o paf.o
HEADERS = main.h ktime.h bseq.h thread.h ksort.h bit_operation.h overlapping.h index.h graph.h paf.h
LIBS = -lm -lz -lpthread
PROG = xRead

.SUFFIXES:.c .o .cc

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

${OBJS}:${HEADERS}

clean:
	@rm -f *.o $(PROG)
	
depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) -- *.c )
