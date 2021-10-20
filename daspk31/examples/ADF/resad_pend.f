      subroutine res_ady(t, y, ady, yp, cj, delta, addelta, 
     &     ires, rpar, ipar, senpar)
      IMPLICIT NONE
      integer ires, ipar(*)
      double precision t, y(*), yp(*), cj, rpar(*)
      double precision senpar(*), ady(*), addelta(*), delta(*)

      call adyres( y, cj, ady, addelta )

      return
      end
