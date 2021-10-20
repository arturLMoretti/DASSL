      subroutine res_ady(t, y, ady, yp, cj, delta, addelta, 
     &     ires, rpar, ipar, senpar)
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     ady(*), addelta(*), delta(*)
      integer ires, ipar(*)

      call adyres( y, cj, ady, addelta )
      return
      end
