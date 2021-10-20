      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
C
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
C
      implicit double precision(a-h,o-z)
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      qsen(1) = 0.0d0
      do i = 1, ns*my*mx
         qsen(1) = qsen(1) + y(i)*y(i)
c         if (ires .eq. 1) qsen(1+i) = 2.0d0*y(i)
      end do
      qsen(2) = qsen(1)
C
      return
      end
