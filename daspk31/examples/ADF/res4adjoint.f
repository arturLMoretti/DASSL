      subroutine res4adjoint(
     *     n, v, 
     *     T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR,
     *     sum)
      implicit none
      integer  n, ires, ipar(*)
      double precision v(n), cj, t, Y(*), yprime(*), delta(*), 
     *     rpar(*), senpar(*), sum
      external RES

      call res(T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)

      sum = 0.0d0
      do i = 1, n
         sum = sum + delta(i)*v(i)
      end do

      return
      end
         
     
