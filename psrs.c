#include "psrs.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include "mpi-printf.h"
#include "seq_sort.h"

void assert_malloc(void *ptr, int _MYPROC) {
    if (ptr==NULL) {
        fprintf(stderr,"ERROR: PE%2d cannot malloc\n",_MYPROC);
        fflush(stderr);
        exit(1);
    }
}

#define OVERSMP 1

void all_sort_psrs_check_i(int *A, int A_size) {
  int i;
  int _PROCS, _MYPROC;

  MPI_Comm_size(MPI_COMM_WORLD, &_PROCS);
  MPI_Comm_rank(MPI_COMM_WORLD, &_MYPROC);
  
  for (i=1 ; i<A_size ; i++)
    if (A[i] < A[i-1])
      fprintf(stderr,"(PE%3d)ERROR: A[%8d] < A[%8d]  (%12d %12d)\n",
	      _MYPROC, i, i-1, A[i], A[i-1]);

  MPI_fprintf(stderr,"HERE CHECKED\n");
  
}

int *findpos (int *base, int size, int key) {
/*  find position of key (or next larger) in sorted array */

  int done = 0, *l, *u, *m, *end;

  l = base;
  u = end = base + size - 1;

  while (! done && size > 0) {
    m = l + size/2;
    if (key == *m)
      done = 1;
    else {
      size /= 2;
      if (key > *m)
	l = m+1;
      else
	u = m-1;
    }
  }
  
  while (m < end && key > *m)
    ++m;

  return (m);
}


int all_sort_psrs_i(int *A, int A_size, int *B) {


  int i, B_size;

  int 
      *recv_cnt,     /* number of keys to receive from each PN */
      *send_cnt,     /* number of keys to send to PNs */
      *recv_off,
      *send_off;

  int no_samples,    /* number of samples to take from sorted data */
      *pivbuffer,    /* array of pivots */
      pivdist,       /* distance between pivots in set of samples */
      possmp,        /* position of pivot in set of all samples */
      **prtbound,    /* boundaries for partitioning local data */
      *recv_buf,     /* incoming data to merge */
      *smpbuffer,    /* array of local samples */
      *smpallbuf,    /* array of global samples */
      smpdist,       /* distance between consecutive samples */
      trecv;

  MPI_Status stat;

  int _PROCS, _MYPROC;

  MPI_Comm_size(MPI_COMM_WORLD, &_PROCS);
  MPI_Comm_rank(MPI_COMM_WORLD, &_MYPROC);

  prtbound = (int **) malloc((_PROCS+1)*sizeof(int *));
  assert_malloc(prtbound, _MYPROC);
  
  send_cnt = (int *)malloc(_PROCS*sizeof(int));
  assert_malloc(send_cnt, _MYPROC);
  
  send_off = (int *)malloc(_PROCS*sizeof(int));
  assert_malloc(send_off, _MYPROC);
  
  recv_off = (int *)malloc(_PROCS*sizeof(int));
  assert_malloc(recv_off, _MYPROC);
  
  recv_cnt = (int *)malloc(_PROCS*sizeof(int));
  assert_malloc(recv_cnt, _MYPROC);

  pivbuffer = (int *)malloc((_PROCS-1)*sizeof(int));
  assert_malloc(pivbuffer, _MYPROC);

  /* STEP 1 */
  fastsort(A, A_size);

  /* STEP 2 */
  no_samples = OVERSMP*_PROCS - 1;
  smpdist = A_size / (no_samples + 1);

  smpbuffer = (int *)malloc(no_samples*sizeof(int));
  assert_malloc(smpbuffer, _MYPROC);

  for (i=0; i<no_samples; i++)
    smpbuffer[i] = A[(i+1)*smpdist];

  /* STEP 3 */
  smpallbuf = (int *)malloc(_PROCS*no_samples*sizeof(int));
  assert_malloc(smpallbuf, _MYPROC);

  MPI_Gather(smpbuffer, no_samples, MPI_INT,
	     smpallbuf, no_samples, MPI_INT,
	     0, MPI_COMM_WORLD);

  /* STEP 4 */
  fastsort(smpallbuf, _PROCS*no_samples);

  if (_MYPROC==0) {
    pivdist = _PROCS*no_samples/(_PROCS-1);
    possmp = pivdist/2;
    for (i=0; i<_PROCS-1; i++) {
      pivbuffer[i] = smpallbuf[possmp];
      possmp += pivdist;
    }
  }

  /* STEP 5 */
  MPI_Bcast(pivbuffer, _PROCS-1, MPI_INT, 0, MPI_COMM_WORLD);

  /* STEP 6 */
  prtbound[0] = A;
  for (i=1; i<_PROCS; i++)
    prtbound[i] = findpos (A, A_size, pivbuffer[i-1]);

  prtbound[_PROCS] = A + A_size;

  free (smpbuffer);
  free (smpallbuf);
  free (pivbuffer);

  /* STEP 7 */
  for (i=0; i<_PROCS; i++)
    send_cnt[i] = max(0, prtbound[i+1] - 1 - prtbound[i]);

  MPI_Alltoall(send_cnt, 1, MPI_INT,
	       recv_cnt, 1, MPI_INT,
	       MPI_COMM_WORLD);
   	
  send_off[0] = 0;
  recv_off[0] = 0;
  for (i=1 ; i<_PROCS ; i++) {
    send_off[i] = send_off[i-1] + send_cnt[i-1];
    recv_off[i] = recv_off[i-1] + recv_cnt[i-1];
  }

  trecv = recv_off[_PROCS-1] + recv_cnt[_PROCS-1];

  if (trecv>0) {
    recv_buf = (int *)malloc(trecv*sizeof(int));
    assert_malloc(recv_buf, _MYPROC);
  }
  else {
    recv_buf = NULL;
  }

  MPI_Alltoallv(A,        send_cnt, send_off, MPI_INT,
		recv_buf, recv_cnt, recv_off, MPI_INT,
		MPI_COMM_WORLD);
  

  /* STEP 8 */
#if 0
  B_size = trecv;
  bcopy(recv_buf, B, B_size*sizeof(int));

  fastsort(B, B_size);
#else
  {
    int **listptr, *listcnt;
    int elem, minelem, minlist;
    int i, j;
    listptr = (int **)malloc(_PROCS * sizeof(int *));
    assert_malloc(listptr,_MYPROC);
    listcnt = (int *)malloc(_PROCS * sizeof(int));
    assert_malloc(listcnt,_MYPROC);

    for (i=0 ; i<_PROCS ; i++) {
      listptr[i] = recv_buf + recv_off[i];
      listcnt[i] = recv_cnt[i];
    }

    B_size = 0;
    while (B_size < trecv) {
      minlist = -1;
      for (j=0 ; j<_PROCS ; j++) {
	if ((listcnt[j] > 0) &&
	    ((minlist<0) || ((elem = *(listptr[j])) < minelem))) {
	  minlist = j;
	  minelem = elem;
	}
      }
      B[B_size++] = minelem;
      listptr[minlist]++;
      listcnt[minlist]--;
    }

    free(listcnt);
    free(listptr);
  }
#endif
  
  MPI_Barrier(MPI_COMM_WORLD);

  free(recv_buf);
  free(recv_off);
  free(send_off);
  free(send_cnt);
  free(recv_cnt);
  free(prtbound);
  return (B_size);
}

