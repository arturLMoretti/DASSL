      subroutine res_ady(t, y, yp, cj, ires, rpar, ipar, senpar, 
     *     ady, addelta)
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     ady(*), addelta(*)
      integer ires, ipar(*)

      call adyres( y, yp, cj, ipar, senpar, ady, addelta )
      return
      end

      subroutine res_adyp(t, y, yp, cj, ires, rpar, ipar, senpar, 
     *     adyp, addelta)
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     adyp(*), addelta(*)
      integer ires, ipar(*)

      call adypres( y, cj, ipar, adyp, addelta )

      return
      end


