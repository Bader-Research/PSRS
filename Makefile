EXECS     = psrs

OPT         = -O2
# OPT       = -g
CFLAGS      = $(OPT)
CC          = mpicc

OBJS	    = psrs.o mpi-printf.o seq_sort.o main.o

$(EXECS): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS) 

.SUFFIXES: .c

.c.o : 
	$(CC) $(CFLAGS) -c $< 

clean:
	rm $(OBJS) $(EXECS) *~