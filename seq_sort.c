#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "seq_sort.h"

#define RADIXSORT_INT_BREAKPT 100

int intcompare(int *i, int *j)
{
    return(*i - *j);
}

unsigned bits(unsigned x, int k, int j) {
/* Return the j bits which appear k bits from the right in x */
    return (x>>k) & ~(~0<<j);
}

void seq_radixsort(int *a, int n) {
/* Radix sort a list of n integers, a[], between 0 and M-1,
   where M = 2^m, and n = 2^w */
    register int
	i,
	j,
	m,
	w,
	M,
	pass;

    int	nextpass,
	*count,
	*bitArr,
	*b;
    bitArr = (int *)malloc(n*sizeof(int));
    if (bitArr==NULL) {
        fprintf(stderr,"ERROR: malloc failed\n");
        fflush(stderr);
        exit(1);
    }
    b = (int *)malloc(n*sizeof(int));
    if (b==NULL) {
        fprintf(stderr,"ERROR: malloc failed\n");
        fflush(stderr);
        exit(1);
    }

    w = sizeof(int) << 3;   /* The number of bits in the key */
    m = w >> 2;             /* m = number of bits per pass   */
    M = 1 << m;             /* The range of each pass */
    
    if ((count = (int*)malloc(M*sizeof(int)))==NULL)
	fprintf(stderr,"ERROR: seq_radixsort count could not be malloc'ed\n");

/* Note that the loop below runs two passes each time so that
   a[] and b[] don't need to be swapped */
    
    for (pass=0 ; pass<(w/m) ; pass+=2) {
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(a[i],pass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) b[--count[bitArr[i]]] = a[i]; 

	nextpass = pass+1;
	for (j=0 ; j<M ; j++) count[j] = 0;
	for (i=0 ; i<n ; i++) count[bitArr[i] = bits(b[i],nextpass*m,m)]++;
	for (j=1 ; j<M ; j++) count[j] += count[j-1];
	for (i=n-1 ; i>=0 ; i--) a[--count[bitArr[i]]] = b[i];
    }
    free(count);
    free(b);
    free(bitArr);
}





void insertsort_i(int *A, int n) {

    register int item;
    register int i,j;
    
    for (i=1 ; i<n ; i++) {
	item = A[i];
	j = i-1;
	while ((j>=0)&&(item < A[j])) {
	    A[j+1] = A[j];
	    j--;
	}
	A[j+1] = item;
    }

}


void fastsort(int* arr, int nel) {

    if (nel>=RADIXSORT_INT_BREAKPT)
	seq_radixsort(arr,nel); 
    else
	insertsort_i(arr,nel);
}


