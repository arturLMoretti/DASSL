      subroutine res_ady(t, y, ady, yp, cj, delta, addelta, ires, 
     *     rpar, ipar, senpar)
      implicit none
      integer ires, ipar(*), ad_q_, adj_op
      real*8  t, y(*), yp(*), cj, delta(*), addelta(*), rpar(*),
     *     senpar(*), ady(*)
      include 'adlib.h'
c...  call ADIfor3.0 generated routine
c ** fwd pass of adjoint
      ad_q_ = 1
      adj_op = adlib_fwd
      call ady_res_proc(adj_op, ad_q_, t, y, ady, yp, cj, delta,
     +     addelta, ires, rpar, ipar, senpar)
c ** rev pass of adjoint
      adj_op = adlib_rev
      call ady_res_proc(adj_op, ad_q_, t, y, ady, yp, cj, delta,
     +     addelta, ires, rpar, ipar, senpar)

      return
      end

      subroutine res_adp(t, y, yp, cj, delta, addelta, ires, rpar,
     *     ipar, senpar,adsenpar)
      implicit none
      integer ires, ipar(*), ad_q_, adj_op
      real*8  t, y(*), yp(*), cj, delta(*), addelta(*), rpar(*),
     *     senpar(*), adsenpar(*)
      include 'adlib.h'
c...  call ADIfor3.0 generated routine
c ** fwd pass of adjoint
      ad_q_ = 1
      adj_op = adlib_fwd
      call adp_res_proc(adj_op, ad_q_, t, y, yp, cj, delta, addelta,
     +     ires, rpar, ipar, senpar, adsenpar)
c ** rev pass of adjoint
      adj_op = adlib_rev
      call adp_res_proc(adj_op, ad_q_, t, y, yp, cj, delta, addelta,
     +     ires, rpar, ipar, senpar, adsenpar)
      
      return
      end
      
      subroutine res_adyp(t, y, yp, adyp, cj, delta, addelta, ires, 
     *     rpar, ipar, senpar)
      implicit none
      integer ires, ipar(*), ad_q_, adj_op
      real*8  t, y(*), yp(*), cj, delta(*), addelta(*), rpar(*),
     *     senpar(*), adyp(*)
      include 'adlib.h'
c...  call ADIfor3.0 generated routine
c ** fwd pass of adjoint
      ad_q_ = 1
      adj_op = adlib_fwd
c$$$      call adyp_res_proc(adj_op, ad_q_, t, y, yp, adyp, cj, delta,
c$$$     +     addelta, ires, rpar, ipar, senpar)
c$$$c ** rev pass of adjoint
c$$$      adj_op = adlib_rev
c$$$      call adyp_res_proc(adj_op, ad_q_, t, y, yp, adyp, cj, delta,
c$$$     +     addelta, ires, rpar, ipar, senpar)

      return
      end
      
      
