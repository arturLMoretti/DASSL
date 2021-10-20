      SUBROUTINE QRES (T,Y,YPRIME,QSEN,ires,RPAR,IPAR,senpar)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(*), YPRIME(*), RPAR(*), IPAR(*),senpar(*), Qsen(*)
C
C Set problem constants using IPAR and RPAR.
C
      qsen(1) = value(20, Y)

      RETURN
      END
      double precision function value(n,v)
      implicit none
      integer n, i, j, k
      double precision v(3,n), prod, dij, dsqrt
      
      prod = 1.0d0
      do i = 1, n
         do j = i+1, n
            dij = 0.0d0
            do k = 1, 3
               dij = dij + (v(k,i)-v(k,j))**2
            end do
            prod = prod*dsqrt(dij)
         end do
      end do
      value = dlog(prod)
      return
      end
