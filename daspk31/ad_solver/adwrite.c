#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "cfortran.h"

static FILE *fp=NULL;

void adopen()
{
  if ((fp=tmpfile()) == NULL) {
    printf(" Cannot open temporary work file.\n");
    exit(1);
  }
}

void adclose()
{
  fclose(fp);
}

FCALLSCSUB0(adopen, ADOPEN, adopen)
FCALLSCSUB0(adclose, ADCLOSE, adclose)

void getsizeof(int label, 
               int* ret)
{
  size_t size;
  if (label == 1) size = sizeof(double);
  else if (label==0) size = sizeof(int);
  *ret = (int) size;
}

FCALLSCSUB2(getsizeof, GETSIZEOF, getsizeof, INT, PINT)

void adwrite(void *buf, 
	     int label, 
	     size_t count,
	     int *irec)
{
  int    i;
  size_t size;

  if (!fp) {
      printf(" Cannot open file.\n");
      exit(1);
  }

  if (label == 1) size = sizeof(double);
  else if (label==0) size = sizeof(int);

  *irec = ftell(fp);
  fwrite(buf, size, count, fp);
}

FCALLSCSUB4(adwrite, ADWRITE, adwrite, PVOID, INT, INT, PINT)

  
void adread (void *buf, 
	     int  label, 
	     size_t count,
	     int irec)
{
  int    i;
  size_t size;
  
  if (!fp) {
      printf(" Cannot open file.\n");
      exit(1);
  }

  if (label == 1) size = sizeof(double);
  else if (label==0) size = sizeof(int);

  fseek(fp, irec, 0);
  fread(buf, size, count, fp);
}
     
FCALLSCSUB4(adread, ADREAD, adread, PVOID, INT, INT, INT)
  
void transpose(const int n,
	      const int m,
	      double *A)
{
  double *B;
  int i, j;
  long dim;

  dim = m*n*sizeof(double);
  if ((B = (double *) malloc(dim)) == NULL) {
    printf(" Memory allocation error in transpose routine \n");
    exit(0);
  } else {
    for (j = 0; j<m; j++) {
      for (i = 0; i<n; i++) {
	B[j+i*m] = A[i+j*n];
      }
    }
    memcpy(A,B,dim);
    free(B);
  }
}
  
FCALLSCSUB3(transpose, TRANSPOSE, transpose, INT, INT, PDOUBLE)  
  
