       program fekete_test
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      DIMENSION Y(160*9), YPRIME(160*9), RWORK(230000),IWORK(10000)
      DIMENSION INFO(30),RPAR(1),IPAR(1), senpar(8), YTRUE(160)
      dimension ysen(8*160), ypsen(160*8), tmp(8*8), qsen(200),b2(160)
      dimension ieopt(160*2), infob(10), qsenb(8), b1(160), v(160)
      dimension delta(200)

      logical consis
      real*4   t1, t2, t3, tt(2), secnds
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-4
      ATOL = 1.0d-4
      LRW = 230000
      LIW = 10000
      T = 0.0d0
      TOUT = 1000.0d0
      neq = 160
*     
*     Initialize y, yprime, senpar
*     
      senpar(1) = 1d0/13d0
      senpar(2) = 3d0/8d0
      senpar(3) = 1d0/29d0
      senpar(4) = 1d0/8d0
      senpar(5) = 1d0/7d0
      senpar(6) = -2d0/15d0
      senpar(7) = 1d0/17d0
      senpar(8) = -0.3d0
      
      call init(neq,t,y,yprime,consis,senpar)
      do i = 1, 64
         tmp(i) = 0.0d0
      end do
      do i = 1, 8
         tmp(i+(i-1)*8) = 1.0d0
      end do
*     
*     Initialize the INFO vector to 0.
*     
      DO 10 I = 1,30
         INFO(I) = 0
 10   CONTINUE
*     
*     Get solution at intermediate steps.
*     
      INFO(3)  = 0
*     
*     Adifor evaluation of the Jacobian. (1 = no,0 = yes)
*     
      info(5) = 3
C     
C     Here set INFO(14) = 1 to get the computed initial values.
      INFO(14) = 1
*     
*     Set #12 #13 as index-2 algebraic variables, does not do the error test on it.
*     
      info(16) = 1
      do i = 1, 120
         iwork(40+i) = 1
         iwork(40+neq+i) = 0
      end do
c
c...  index-2 variables and constraints
      do i = 121, 160
         iwork(40+i) = -2
         iwork(40+i+neq) = 1
      end do

      iwork(38) = 0
      Ny = neq
c
c...  initialization options
      INFO(11) = 5              ! 4 or 5 
      if (consis) info(11) = 0     
      neq = Ny
c     
c...  IEOPT 
      do i = 1, 120
         ieopt (i) = 1
      end do
      do i = 121, 160
         ieopt(i) = -2
      end do
      do i = 1, neq
         ieopt(i+neq) = iwork(40+i) ! forward ID information
      end do
c
      do i = 1, 10
         infob(i) = 0
      end do
      infob(1) = 0              ! error tolerance
      infob(2) = 2              ! 0 -- index-0 initialization
                                ! 1 -- index-1 initialization
                                ! 2 -- index-2 initialization for ODE
                                ! 3 -- user input initialization
      infob(3) = 1              ! 0 -- input RES_ADY  in G_RES
                                ! 1 -- input RES_ADYP in ADRES
                                ! 2 -- Input G_RES_ADYP in ADRES
                                ! 3 -- my adres
      nquad    = 1              ! number of derived functions
      infob(4) = 0              ! take it as variables or do the quadrature
                                ! independently
      infob(5) = 1              ! index-2 indicator
      info(18) = 0              ! error control for the quadrature variables
      info(19) = 8              ! number of the parameters
      info(20) = 0              ! finite difference
      info(22) = 0              ! no parameter in the RES routine
      info(30) = 0              ! number of quadrature variables
      neq = neq + info(30)
      nbuf = 60*neq
      t = 0.0d0
      call xspini

      do i = 1, neq
         v(i) = 1.0d0
      end do

      do i = 1, 3*neq
         iwork(i) = 0
      end do
      cj = 1.0d0

      t1 = dtime(tt)
c$$$      t1 = secnds(0.0e0)
c$$$      do i = 1, 500
c$$$         call vmulFx3(neq,v,t,y,yprime,cj, delta, ires, 
c$$$     *     rpar, ipar, senpar, b1)
c$$$         print *, i, secnds(t1)
c$$$      end do
c$$$      t2 = dtime(tt)
c$$$      write(6,*) '   Time it took:', tt(1)+tt(2)
      t1 = secnds(0.0e0)
      do i = 1, 500         
         call vmulFx2(neq,v,t,y,yprime,cj, delta, ires, 
     *     rpar, ipar, senpar, b2, iwork)
         print *, i, secnds(t1)
      end do 
      t2 = dtime(tt)
      write(6,*) '   Time it took:', tt(1)+tt(2)

      STOP
      END


      subroutine vmulFx2(n,v,t,y,yprime, cj, delta, ires, 
     *     rpar, ipar, senpar, b, iwork)
      implicit none
      integer  n, ires, ipar(*), iwork(*)
      double precision v(*), cj, t, Y(*), yprime(*), delta(*), 
     *     rpar(*), senpar(*), b(*) 
c      
      integer i,j, jj, i_y, i_yp, i_delta
      integer info, outlen, indvec(200)
      double precision row(200)

      i_y = 0
      i_yp = i_y+n 
      i_delta = i_yp + n
         
      do i = 1, n
         call dspsd(iwork(i), i, 1.0d0, 1)
         call dspsd(iwork(n+i), i, 0.0d0, 1)
      end do
c
      call i_res (
     $     t, y, iwork(i_y+1), yprime, iwork(i_yp+1), cj, delta, 
     $     iwork(i_delta+1), ires, rpar, ipar, senpar)

      do j = 1, n
         b(j) = 0.0d0 
      end do
      do i = 1, n
         call dspxsq(indvec,row,n,iwork(i_delta+i),
     $        outlen,info)
         do j = 1, outlen
            jj = indvec(j)
            b(jj) = b(jj) + v(i)*row(j)
         end do
      end do
c 
      return
      end

      subroutine vmulFx3(n,v,t,y,yprime,cj, delta, ires, 
     *     rpar, ipar, senpar, b)

      implicit none
      integer  n, ires, ipar(*)
      double precision v(n), cj, t, Y(*), yprime(*), delta(*), 
     *     rpar(*), senpar(*), b(*), ady(200)

      integer i,j

      do i = 1, n
         b(i) = 0.0d0
         ady(i) = v(i)
      end do

      call res_ady(t, y, yprime, cj, ires, rpar, ipar, senpar,
     *     b, ady)

      return
      end

      

      
