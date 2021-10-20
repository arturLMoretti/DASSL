       program testad
*  
*      *** Demo program for DASPK ***
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, adresp, qres, adjac, jac
      external res_adyp, res_ady, g_res_ady
      DIMENSION Y(10),RWORK(3000),IWORK(1000),YPRIME(10), infob(10)
      DIMENSION INFO(30),RPAR(1),IPAR(20), senpar(1), ieopt(10)
      dimension qsen(30)
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-6
      ATOL = 1.0d-6
      LRW = 3000
      LIW = 1000
      T = 0.0d0
      TOUT = 1.0d0
      neq = 5
*     
*     Initialize y, yprime, senpar
*     
      call init(y,yprime,senpar)
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
*     Finite difference approx. (FDA) of Jacobian. (1 = no,0 = yes)
*     
      info(5) = 1
c
c...  exclude the algebraic variables from the error test
      info(16) = 1
      NEQ = 5
      do i = 1, 3
         iwork(40+i) = 1        ! differential variables
         iwork(40+5+i) = 0      ! differential equations
      end do

      iwork(44) = -2            !  index-2 variables
      iwork(45) = -2            !  index-2 variables
      iwork(49) = 1             !  index-2 constraint
      iwork(50) = 1             !  index-2 constraint
      iwork(38) = 0
*     
      info(19) = 1
      info(20) = 2
      info(14) = 1
      INFO(11) = 5
      neq = 10
      print *, ' do you want forward sensitivity with g(x):(0,1)'
      read(*,*) isen
 15   continue
      CALL DDASPK(
     *     RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *     LRW,IWORK,LIW,RPAR,IPAR,jac, psol,senpar, H_res, 
     *     K_RES, G_G_RES)
      if (info(11) .gt. 0) then
         info(11) = 0
         if (isen .eq. 0) goto 16
         goto 15
      end if
      gp = 0.0
      do i = 1, 5
         gp = gp + y(5+i)
      enddo
      print *, 'gp in forward mode =', gp, 4.*exp(2.*t)
      print *, (iwork(i),i=11,15)
      stop
 16   continue
      info(11) = 5
      neq = 5
c     
c...  IEOPT 
      do i = 1, 3
         ieopt (i) = 1
      end do
      ieopt(4) = -2
      ieopt(5) = -2
      do i = 1, neq
         ieopt(neq+i) = iwork(40+i) ! forward ID information
      end do
      infob(1) = 0              ! error tolerance
      infob(2) = 2              ! index-2 initialization for ODE
      infob(3) = 1              ! 0 -- ODE form
                                ! 1 -- DAE form with F_xt = const
                                ! 2 -- Input adres by user
                                ! 3 -- my adres
      nquad    = 1
c     print *, ' Input number of the derived functions(0,1,2):'
c     read(*,*) nquad
      infob(4) = 0              ! take it as variables or do the quadrature
                                ! independently 
      infob(5) = 1              ! index-2 case
      infob(6) = 1  
      infob(7) = 0
      info(18) = 0
      info(19) = 1
      info(20) = 2
      info(22) = 0              ! no parameters in the RES routine
      info(28) = 0
      neq = neq + info(28)
      nbuf = 200*neq
      t = 0.0d0

 20   IF(T .LT. TOUT) THEN 
*     
         if (infob(3) .eq. 3) then
            CALL DDASPKadjoint(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,jac, psol,senpar, 
     *           adresp, nquad, qres, qsen, infob, rtolb, atolb, nbuf,
     *           adjac, ieopt, res_adp, RES_ady,g_res_ady,
     *           adinit, k_res)
            if (idid.gt.0) GOTO 20
         else
            CALL DDASPKadjoint(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,jac, psol,senpar, 
     *           res_adyp, nquad, qres, qsen, infob, rtolb,atolb,nbuf,
     *           adjac, ieopt, res_adp, RES_ady, g_res_ady, 
     *           adinit, K_RES)
            if (idid.gt.0) GOTO 20
         end if            
      ENDIF 
*     
*     Output IWORK and RWORK vector values
      print *, 'g_p in backward mode   =', qsen(nquad*(neq+1)+1) 
      print *, 'The exact value of g_p =', 4.*exp(2.*t)
*     
      END
*------------------------------------------------------------------------*
      SUBROUTINE RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR, senpar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), DELTA(*), RPAR(*), IPAR(*),senpar(*) 
      cjl = 1.0d0
      if (ires .eq. 0) then
         DELTA(1) = YPRIME(1) - y(1) - y(4) - y(5)
         DELTA(2) = YPRIME(2) - 2.0d0*y(3) - y(4) + y(5)
         DELTA(3) = YPRIME(3) - y(1) - y(2)
         DELTA(4) = cjl*(y(1) - y(2))
         DELTA(5) = cjl*(y(1) + y(2) - 2.0d0*y(3))
      else if (ires .eq. 1) then
         DELTA(6) = YPRIME(6) - y(6) - y(9) - y(10)
         DELTA(7) = YPRIME(7) - 2.0d0*y(8) - y(9) + y(10)
         DELTA(8) = YPRIME(8) - y(6) - y(7)
         DELTA(9) =  cjl*(y(6) - y(7))
         DELTA(10) = CJl*(y(6) + y(7) - 2.0d0*y(8))
      else if (ires .eq. 2) then
         DELTA(4) = YPRIME(1) - yprime(2)
         delta(5) = yprime(1) + yprime(2) - 2.0d0*yprime(3)
         delta(9) = yprime(6) - yprime(7)
         delta(10) = yprime(6) + yprime(7) - 2.0d0*yprime(8)
      end if
      RETURN
      END 
*------------------------------------------------------------------------*

      subroutine jac(t, y, yp, pd, cj, rpar, ipar, senpar, ijac)
      implicit none
      double precision t, y(*), yp(*),pd(5,5),cj,rpar(*),senpar(*),cjl
      integer ipar(*), ijac, i, j

      do i = 1, 5
         do j = 1,5
            pd(i,j) = 0.0d0
         end do
      end do
      
      cjl = 1.0d0
      if (ijac .eq. 0 .or. ijac.eq.1) then
         pd(1,1) = cj - 1.0
         pd(1,4) = -1.0
         pd(1,5) = -1.0
         pd(2,2) = cj - 2.
         pd(2,4) = -1.
         pd(2,5) =  1.
         pd(3,1) = -1.
         pd(3,2) = -1.
         pd(3,3) =  cj
c$$$         pd(4,1) = cj
c$$$         pd(4,2) = -cj
c$$$         pd(5,1) = cj
c$$$         pd(5,2) = cj
c$$$         pd(5,3) = -2.*cj
         pd(4,1) = cjl
         pd(4,2) = -cjl
         pd(5,1) = cjl
         pd(5,2) = cjl
         pd(5,3) = -2.*cjl

      else if (ijac .eq. 2) then
         pd(1,1) = cj
         pd(1,4) = -1.0
         pd(1,5) = -1.0
         pd(2,2) = cj
         pd(2,4) = -1.
         pd(2,5) =  1.
         pd(3,3) =  cj
         pd(4,1) =  cj
         pd(4,2) = -cj
         pd(5,1) = cj
         pd(5,2) = cj
         pd(5,3) = -2.*cj         
      end if

      return
      end

      SUBROUTINE INIT(Y,YPRIME,senpar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), senPAR(*)
*     
*     Set y, yprime, IPAR, RPAR to 0, initially.
*     
      do i = 1,10
         y(i) = 1.0d0
         yprime(i) = 0.0d0
      enddo
      y(4) = 0.5d0
      y(5) = 0.5d0
      y(9) = 0.5d0
      y(10) = 0.5d0
*     
*     
*     Parameter values - placed somewhere in senPAR.
*     
      senpar(1) = 1.0d0
      return
      end
c
      subroutine adinit(t, neqad, nq, Y,YP,ady,adyp,
     *        Qsen,RPAR,IPAR,senpar)
      implicit none
      integer neqad, nq, ipar(*), i, j, iq, neqst, ires
      double precision t, y(*),yp(*),ady(*), adyp(*), rpar(*),
     *     qsen(*), senpar(*), delta(5)

      ady(1) = 0.0d0
      ady(2) = 0.0d0
      ady(3) = 4.0d0
      ady(4) = 0.0d0
      ady(5) = -4.0d0
      
      do i = 1, 5
         adyp(i) = 0.0d0
      end do
      call adresp(
     *     t, adY, adYP, 1.0d0, delta, IRES, RPAR,IPAR,SENPAR, y)
      do i = 1, 5
         adyp(i) = -delta(i)
      end do

      return
      end


      SUBROUTINE adresp(
     *     t, adY, adYP, CJ, delta, IRES, RPAR,IPAR,SENPAR, adi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ady(*), adyp(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)
      dimension adi(*)

      delta(1) = adyp(1) - ady(1) - ady(3) - ady(4) - ady(5)
      delta(2) = adyp(2) - (2.0d0*ady(2) + ady(3) - ady(4) + ady(5))
      delta(3) = adyp(3) + 2.0d0*ady(5) 
      delta(4) = cj*(ady(1) + ady(2))
      delta(5) = cj*(ady(1) - ady(2))
         
      return
      end
      subroutine adjac(t, y, yp, pd, cj, rpar, ipar, senpar, ijac)
      implicit none
      double precision t, y(*), yp(*),pd(5,5),cj,rpar(*),senpar(*)
      integer ipar(*), ijac, i, j

      do i = 1, 5
         do j = 1,5
            pd(i,j) = 0.0d0
         end do
      end do
      
      if (ijac .eq. 0) then
         pd(1,1) = cj -1.0
         pd(1,3) = -1.0d0
         pd(1,4) = -1.
         pd(1,5) = -1.
         pd(2,2) = cj - 2.
         pd(2,3) = -1.
         pd(2,4) = 1.0d0
         pd(2,5) = -1.0
         pd(3,3) = cj
         pd(3,5) = 2.
         pd(4,1) = cj
         pd(4,2) = cj
         pd(5,1) = cj
         pd(5,2) = -cj
      else if (ijac .eq. 2) then
         pd(1,1) = cj
         pd(1,4) = -1.
         pd(1,5) = -1.
         pd(2,2) = cj
         pd(2,4) =  1.0d0
         pd(2,5) = -1.0
         pd(3,3) = cj
         pd(3,5) = 2.
         pd(4,1) = cj
         pd(4,2) = cj
         pd(5,1) = cj
         pd(5,2) = -cj
      end if
      return
      end
      
      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i

      qsen(1) = y(1) + y(2) + y(3) + y(4) + y(5)
      if (ires .eq. 1) then
         qsen(2) = 1.0d0
         qsen(3) = 1.0d0
         qsen(4) = 1.0d0
         qsen(5) = 1.0d0
         qsen(6) = 1.0d0
      end if
      return
      end
      subroutine res_ady(t, y, yp, cj, ires, rpar, ipar, senpar, 
     *     delta, ady)
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     ady(*), delta(*)
      integer ires, ipar(*),i

      delta(1) = - ady(1) - ady(3) - cj*ady(4) - cj*ady(5)
      delta(2) = - (2.0d0*ady(2) + ady(3) - cj*ady(4) + cj*ady(5))
      delta(3) = 2.0d0*ady(5)*cj 
      delta(4) = (ady(1) + ady(2))
      delta(5) = (ady(1) - ady(2))
         
      return
      end

      subroutine res_adyp(t, y, yp, cj, ires, rpar, ipar, senpar, 
     *     adyp, addelta)
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     adyp(*), addelta(*)
      integer ires, ipar(*)

      adyp(1) = addelta(1)
      adyp(2) = addelta(2)
      adyp(3) = addelta(3)
      adyp(4) = 0.0d0
      adyp(5) = 0.0d0

      return
      end

      subroutine g_res_ady(t, g_t, y, g_y,yp, cj, ires, rpar,ipar,senpa
     *r, ady, g_ady, addelta, g_addelta)
        implicit none
        integer ires, ipar(*), i
        double precision t, y(*), yp(*), cj, rpar(*)
        double precision senpar(*), ady(*), addelta(*)
C
        double precision g_t, g_y(*), g_ady(*), g_addelta(*)

        g_ady(1) =-cj*g_addelta(5)-cj*g_addelta(4)-
     *       g_addelta(3)-g_addelta(1)
        g_ady(2) = (-cj*g_addelta(5))+cj*g_addelta(4)-g_addelta(3)-
     *       2*g_addelta(2)
        g_ady(3) = 2*g_addelta(5)*cj
        g_ady(4) = (g_addelta(2)+g_addelta(1))
        g_ady(5) = (-(g_addelta(2)*cj))+g_addelta(1)

        return
        end
        
      

