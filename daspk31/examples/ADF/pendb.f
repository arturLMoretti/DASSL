C Copyright 2000 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  FEKETE_ADJOINT
C***REFER TO  DDASPK, DASPKADJOINT
C***DATE WRITTEN   000725   (YYMMDD)
C
C***AUTHORS  Shengtai Li
C            University of California
C            Santa Barbara, CA  93106, USA
C***DESCRIPTION
C
C-----------------------------------------------------------------------
C Demo program for DDASPKADJOINT.
C Index-2 DAE system from single pendulum problem
C
C x1' = x3                      x1(0) = 0.5
C x2' = x4                      x2(0) = -sqrt(p1^2 - x1^2)
C x3' = -x1*x5                  x3(0) = 0.0
C x4' = -x2*x5 + g              x4(0) = 0.0
C 0   = x1*x3 + x2*x4
C
C The sensitivity parameter is p1. 
C The derived function is the kinetic energy 
C   g(1) = 0.5*x3**2 + 0.5*x4**2
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C***END PROLOGUE  FEKETE_ADJOINT
C
      program pendb 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, j_res, adresp, qres, adjac, jac,G_G_RES
      external res_adyp, res_ady, g_res_ady
      DIMENSION Y(10),RWORK(20000),IWORK(1000),YPRIME(10), infob(10)
      DIMENSION INFO(30),RPAR(1),IPAR(20), senpar(1), ieopt(10)
      dimension qsen(30)
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-6
      ATOL = 1.0d-6
      LRW = 20000
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
*     Adifor evaluation of the Jacobian. 
*     
      info(5) = 3
      info(5) = 3
*     
*     take out the algebraic variable from the error test.
*     
      info(16) = 1
      NEQ = 5
c
c     set the equation and variable property
      do i = 1, 4
         iwork(40+i) = 1
         iwork(40+5+i) = 0
      end do
c
c     fixed x1 and x2
      iwork(41) = 2
      iwork(42) = 2
c
c     x5 is an index-2 variable
      iwork(45) = -2
c
c     equation 5 is an index-2 constraint
      iwork(50) = 1
c
c     number of the sensitivity parameters
      info(19) = 1
c
c     evaluation method of the derivative g_x g_p 
      info(20) = 2              ! 2 -- user input
                                ! 0,1 -- finite difference
                                ! 3 -- ADIFOR
c
c     obtain the consisten initial condition first
      info(14) = 1
c
c     initialization for the index-2 system
      INFO(11) = 4
      neq = 10
      CALL DDASPK(
     *     RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *     LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, H_res, 
     *     K_RES, G_G_RES)
c
c     the index-2 constraint is already satisfied.
      info(11) = 5
      neq = 5
c     
c...  IEOPT 
      do i = 1, neq
         ieopt (i) = 1
      end do
c
c     index-2 contraints
      ieopt(5) = -2
      do i = 1, neq
         ieopt(i+neq) = iwork(40+i) ! forward ID information
      end do
c
c     options for the initialization of the adjoint system
c
      infob(2) = 2              ! 0 -- index-0 initialization
                                ! 1 -- index-1 initialization
                                ! 2 -- index-2 initialization for ODE
                                ! 3 -- user input initialization
c
c     options for the evaluation of the adjoint system
c
      infob(3) = 1              ! 0 -- input RES_ADY  in G_RES
                                ! 1 -- input RES_ADYP in ADRES
                                ! 2 -- Input G_RES_ADYP in ADRES
                                ! 3 -- my adres
c
c     Is the original system an index-2 system or not? 
c
      infob(5) = 1              ! index-2 
c
c     do you want to compute the sensitivities with the initial
c     conditions inside DASPKADJOINT?
c
      infob(6) = 1              ! 1 -- yes
                                ! 0 -- no
c
c     estimate of the number of time steps
      infob(7) = 0
c
c     number of the derived functions
c
      nq       = 1

C
C     Buffer size 
      nbuf = 65*neq
      t = 0.0d0
c
      if (infob(3) .eq. 3) then
         CALL DDASPKadjoint(
     *        RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *        LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, 
     *        adresp, nq, qres, qsen, infob, rtolb, atolb, nbuf,
     *        j_res, ieopt, res_adp, RES_ady,g_res_ady, g_res_adyp,
     *        adinit, k_res)
      else
         CALL DDASPKadjoint(
     *        RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *        LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, 
     *        res_adyp, nq, qres, qsen, infob, rtolb,atolb,nbuf,
     *        j_res, ieopt, res_adp, RES_ady, g_res_ady, g_res_adyp,
     *        adinit, K_RES)
      end if            
c
      print *, ' sensitivity g_p = ', qsen(nq*(1+neq)+1)
C------  End of main program for PENDB example program -----------    
      END
*------------------------------------------------------------------------*
      SUBROUTINE RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR, senpar)
C
C This is the user-supplied RES subroutine for this example.
C It computes the residuals for the index-2 pendulum equations
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), DELTA(*), RPAR(*), IPAR(*),senpar(*)     
      if (ires .eq. 0) then
c
c     state equations
         DELTA(1) = YPRIME(1) - Y(3)
         DELTA(2) = YPRIME(2) - Y(4)
         DELTA(3) = YPRIME(3) + Y(1)*Y(5)
         DELTA(4) = YPRIME(4) + Y(2)*Y(5) + 1.0D0
         DELTA(5) = cj*(Y(1)*Y(3) + Y(2)*Y(4))
      else if (ires .eq. 1) then
c
c     sensitivity equations
         DELTA(6) = YPRIME(6) - Y(8)
         DELTA(7) = YPRIME(7) - Y(9)
         DELTA(8) = YPRIME(8) + Y(6)*Y(5) + Y(1)*Y(10)
         DELTA(9) = YPRIME(9) + Y(7)*Y(5) + Y(2)*Y(10)
         DELTA(10) = cj*(Y(1)*Y(8) + Y(6)*Y(3) + Y(2)*Y(9) + Y(7)*Y(4))
      else if (ires .eq. 2) then
c
c     time derivatives of the index-2 constraints
         delta(5) = Y(1)*Yprime(3) + YPrime(1)*y(3) +
     *        Y(2)*yprime(4) + Yprime(2)*y(4)
         delta(10) = yprime(1)*y(8) + y(1)*yprime(8) + yprime(6)*y(3)
     *        + y(6)*yprime(3) + yprime(2)*y(9) + y(2)*yprime(9) +
     *          y(7)*yprime(4) + yprime(7)*y(4)
      end if
      RETURN
      END 
*------------------------------------------------------------------------*
      subroutine adinit(t, neqad, nq, Y,YP,ady,adyp,
     *        Qsen,RPAR,IPAR,senpar)
c
c This routine computes the consistent initial conditions for
c the adjoint index-2 system
c
      implicit none
      integer neqad, nq, ipar(*), i, j, iq, neqst, ires
      double precision t, y(*),yp(*),ady(*), adyp(*), rpar(*),
     *     qsen(*), senpar(*), delta(5)

      ady(1) = 0.0d0
      ady(2) = 0.0d0
      ady(3) = y(3)
      ady(4) = y(4)
      ady(4) = -y(1)*y(3)/y(2)
      ady(5) =(y(1)*ady(1)+y(2)*ady(2) - yp(1)*ady(3) - yp(2)*ady(4))/
     *     (y(1)**2 + y(2)**2)
      
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
c
c This is the ADRES routine for the adjoint system. 
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ady(*), adyp(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)
      dimension adi(*)

      delta(1) = adyp(1) + adi(5)*ady(3) + adi(3)*ady(5)
      delta(2) = adyp(2) + adi(5)*ady(4) + adi(4)*ady(5)
      delta(3) = adyp(3) - ady(1) + adi(1)*ady(5)
      delta(4) = adyp(4) - ady(2) + adi(2)*ady(5)
      delta(5) = cj*(adi(1)*ady(3) + adi(2)*ady(4))
         
      return
      end
      
      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
C
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
C
c
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i
c
c     derived function
      qsen(1) = 0.5d0*y(3)*y(3) + 0.5d0*y(4)*y(4)
      if (ires .eq. 1) then
c
c     derivative of the derived function with respect to y
         qsen(2) = 0.0d0
         qsen(3) = 0.0d0
         qsen(4) = y(3)
         qsen(5) = y(4)
         qsen(6) = 0.0d0
      end if
      return
      end

      subroutine adjac(t, y, yp, pd, cj, rpar, ipar, senpar, ijac, yt)
c
c This is the user-supplied ADJAC routine for this example.
C It computes the Jacobian matrix for the index-2 adjoint system
C
      implicit none
      double precision t, y(*), yp(*),pd(5,5),cj,rpar(*),senpar(*),yt(*)
      integer ipar(*), ijac, i, j

      do i = 1, 5
         do j = 1,5
            pd(i,j) = 0.0d0
         end do
      end do
      
      if (ijac .eq. 0) then
c
c     Jacobian for the integration
         pd(1,1) = cj
         pd(1,3) = yt(5)
         pd(1,5) = yt(3)
         pd(2,2) = cj
         pd(2,4) = yt(5)
         pd(2,5) = yt(4)
         pd(3,1) = -1.0d0
         pd(3,3) = cj
         pd(3,5) = yt(1)
         pd(4,2) = -1.0d0
         pd(4,4) = cj
         pd(4,5) = yt(2)
         pd(5,3) = cj*yt(1)
         pd(5,4) = cj*yt(2)
      else if (ijac .eq. 2) then
c
c     Jacobian for the initialization
         pd(1,1) = cj
         pd(1,5) = yt(3)
         pd(2,2) = cj
         pd(2,5) = yt(4)
         pd(3,3) = cj
         pd(3,5) = yt(1)
         pd(4,4) = cj
         pd(4,5) = yt(2)
         pd(5,3) = cj*yt(1)
         pd(5,4) = cj*yt(2)         
      end if

      return
      end
      subroutine jac(t, y, yp, pd, cj, rpar, ipar, senpar, ijac)
c
c This is the user-supplied JAC routine for this example.
C It computes the Jacobian matrix for the index-2 system
C

      implicit none
      double precision t, y(*), yp(*),pd(5,5),cj,rpar(*),senpar(*)
      integer ipar(*), ijac, i, j

      do i = 1, 5
         do j = 1,5
            pd(i,j) = 0.0d0
         end do
      end do
      
      if (ijac .eq. 0) then
c
c     Jacobian for the integration
         pd(1,1) = cj
         pd(1,3) = -1.0d0
         pd(2,2) = cj
         pd(2,4) = -1.0d0
         pd(3,1) = y(5)
         pd(3,3) = cj
         pd(3,5) = y(1)
         pd(4,2) = y(5)
         pd(4,4) = cj
         pd(4,5) = y(2)
         pd(5,1) = cj*y(3)
         pd(5,2) = cj*y(4)
         pd(5,3) = cj*y(1)
         pd(5,4) = cj*y(2)
      else if (ijac .eq. 1) then
c
c     Jacobian for the initialization of the first stage
         pd(1,1) = cj
         pd(1,3) = -1.0d0
         pd(2,2) = cj
         pd(2,4) = -1.0d0
         pd(3,3) = cj
         pd(3,5) = y(1)
         pd(4,4) = cj
         pd(4,5) = y(2)
         pd(5,3) = cj*y(1)
         pd(5,4) = cj*y(2)
      else if (ijac .eq. 2) then
c
c     Jacobian for the initialization of the second stage
         pd(1,1) = cj
         pd(2,2) = cj
         pd(3,3) = cj
         pd(3,5) = y(1)
         pd(4,4) = cj
         pd(4,5) = y(2)
         pd(5,1) = cj*y(3)
         pd(5,2) = cj*y(4)
         pd(5,3) = cj*y(1)
         pd(5,4) = cj*y(2)         
      end if
      return
      end

      SUBROUTINE INIT(Y,YPRIME,senpar)
C
C This routine computes and loads the vector of initial values.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), SENPAR(*)
*     
*     Set y, yprime, IPAR, RPAR to 0, initially.
*     
      do i = 1,10
         y(i) = 0.0d0
         yprime(i) = 0.0d0
      enddo
*          
*     Parameter values - placed somewhere in SENPAR.
*     
      senpar(1) = 1.0d0
*     
*     Specify the initial values for the SP problem.
*     
      pi = 4.0d0*datan(1.0d0)
      theta = pi/6.0d0
      Y(1) = senpar(1)*dsin(theta)
      Y(2) = -senpar(1)*dcos(theta)
      y(3) = 5.0d0
      y(4) = 5.0d0
C     
      y(5) = 0.0d0
      y(6) = dsin(theta)
      y(7) = -dcos(theta)
c     ... the state variables that satisfy the constraints.
c      y(3) = 11.829694d0
c      y(4) = -y(1)*y(3)/y(2)
*     y(5) = y(3)*y(3) + y(4)*y(4) - y(2) 
*     
*    
c      y(8) = 15.775812
c      y(9) = -27.3165081
      return
      end
c
      subroutine res_adyp(t, y, yp,adyp,cj, delta, addelta, 
     *     ires, rpar, ipar, senpar)
C-------------------------------------------------------------------
C This is the user-supplied RES_ADP subroutine for this example.
C It computes a vector-matrix product and return it in ADYP, i.e.
C     adyp = addelta * F_y'
C--------------------------------------------------------------------
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     adyp(*), addelta(*), delta(*)
      integer ires, ipar(*)

      adyp(1) = addelta(1)
      adyp(2) = addelta(2)
      adyp(3) = addelta(3)
      adyp(4) = addelta(4)
      adyp(5) = 0.0d0

      return
      end


      

