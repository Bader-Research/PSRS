#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "psrs.h"

#define N 1024

main(int argc, char **argv) {
  int A[N];
  int B[2*N];
  int MYPROC;
  int i, Bsize;

  MPI_Init(&argc,&argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &MYPROC);
  srandom(MYPROC);

  for (i=0 ; i<N ; i++)
    A[i] = random();
    
  Bsize = all_sort_psrs_i(A, N, B);
  all_sort_psrs_check_i(B, Bsize);
  
  MPI_Finalize();
}
