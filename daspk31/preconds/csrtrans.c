#include <stdio.h>
#include "cfortran.h"

PROTOCCALLSFSUB9(CSRCSC, csrcsc, INT, INT, INT, PDOUBLE, PINT, PINT, \
		 PDOUBLE, PINT, PINT)

#define CSRCSC(n,jb,ipos,a,ja,ia,ao,jao,iao) \
        CCALLSFSUB9(CSRCSC, csrcsc, INT, INT, INT, PDOUBLE, PINT, PINT, \
		 PDOUBLE, PINT, PINT, n,jb,ipos,a,ja,ia,ao,jao,iao)

PROTOCCALLSFSUB7(ATOB, atob, INT, PDOUBLE, PINT, PINT, PDOUBLE, PINT, PINT)

#define ATOB(n,a,ja,ia,ao,jao,iao) \
     CCALLSFSUB7(ATOB, atob, INT, PDOUBLE, PINT, PINT, PDOUBLE, PINT, PINT)

void csrcscbyc(int neq, 
	       double *a, 
	       int    *ja, 
	       int    *ia)
{
  int nnz;
  double *ao;
  int *jao, *iao;

  nnz = ia[neq];
  ao = (double *)malloc(nnz*sizeof(double));
  jao  = (int *)malloc(nnz*sizeof(int));
  iao  = (int *)malloc((neq+1)*sizeof(int));

  if (jaco == NULL || jao == NULL || iao == NULL) {
    printf(" Memory allocation erroy in csccsrbyc routine.\n");
    exit(-1);
  }

  CSRCSC (neq,1,1,a,ja,ia,ao,jao,iao);

  ATOB (neq, ao, jao, iao, a, ja, ia);
}

FCALLSCSUB4(csrcscbyc, CSRCSCBYC, csrcscbyc, INT, PDOUBLE, PINT, PINT)
    
  
