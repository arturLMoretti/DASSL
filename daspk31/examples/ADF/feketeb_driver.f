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
C Example program for DDASPKADJOINT.
C Index-2 DAE system derived from an optimization problem.
C 
C The detail about this problem can be seen from
C
c    http://www.cwi.nl/cwi/projects/IVPtestset/
C 
C This is the double precision version.
C
C-----------------------------------------------------------------------
C  sensitivity analysis
C-----------------------------------------------------------------------
C   For sensitivity analysis, we take eight parameters in the intial
C   conditions as our sensitivity parameters. The derived function is
C   the cost function of the optimization problem, which is
C      g(x) = \prod_{i<j} (|x_i - x_j|^2)
C
C   Because this is a global optimization problem, the sensitivity should
C   be zero or close to zero (within the tolerance).
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C***END PROLOGUE  FEKETE_ADJOINT
      program fekete_adjoint
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, res_adyp, j_res, qres, res_ady, g_res_ady, h_qres
      DIMENSION Y(160*9), YPRIME(160*9), RWORK(300000),IWORK(2000)
      DIMENSION INFO(30),RPAR(1),IPAR(1), senpar(8), YTRUE(160)
      dimension ysen(8*160), ypsen(160*8), tmp(8*8), qsen(200)
      dimension ieopt(160*2), infob(10), qsenb(8)

      logical consis
      real*4   t1, t2, t3, tt(2)
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-4
      ATOL = 1.0d-4
      LRW = 300000
      LIW = 2000
      T = 0.0d0
      TOUT = 1000.0d0
      neq = 160
*     
*     Initialize sensitivity parameters
*     
      senpar(1) = 1d0/13d0
      senpar(2) = 3d0/8d0
      senpar(3) = 1d0/29d0 
      senpar(4) = 1d0/8d0
      senpar(5) = 1d0/7d0
      senpar(6) = -2d0/15d0
      senpar(7) = 1d0/17d0
      senpar(8) = -0.3d0
c
c     initialize y and yprime
c
      call init(neq,t,y,yprime,consis,senpar)
      do i = 1, 64
         tmp(i) = 0.0d0
      end do
      do i = 1, 8
         tmp(i+(i-1)*8) = 1.0d0
      end do
c
c...  calculate the initial condition for the sensitivity variables
c     by calling ADIFOR generated routine
      call l_init(8,neq,t,y,ysen,8,yprime,ypsen,8,consis,senpar,
     *     tmp, 8)
      do i = 1, 8
         do j = 1, neq
            y(i*neq + j) = ysen(i+(j-1)*8)
            yprime(i*neq + j) = ypsen(i+(j-1)*8)
         end do
      end do
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
C     
C     Here set INFO(14) = 1 to get the consistent initial values.
      INFO(14) = 1
C--------------------------------------------------------------------
C     Set up the iwork(*) array for the information about
C     the variables and equations.
C--------------------------------------------------------------------
C
C     First taking the algebraic variables out of the error test
C
      info(16) = 1
C
C     Set all the variables as differential variables, set all the 
C     equations as differential equations.
      do i = 1, 120
         iwork(40+i) = 1
         iwork(40+neq+i) = 0
      end do
c
c...  indicate index-2 variables and constraints
      do i = 121, 160
         iwork(40+i) = -2
         iwork(40+i+neq) = 1
      end do
c
c...  initialization for the index-2 system
      INFO(11) = 5              ! 4 or 5 
      if (consis) info(11) = 0     
c     
c...  IEOPT 
      do i = 1, 120
         ieopt (i) = 1          ! differential equation
      end do
      do i = 121, 160
         ieopt(i) = -2          ! index-2 constraint
      end do
      do i = 1, neq
         ieopt(i+neq) = iwork(40+i) ! forward ID information
      end do
C
C     Here we set up the INFOB array, which describes the various options
C     in the way we want to solve the adjoint equations.
C
C     We first initialize the entire INFOB array to zero, then set select
C     entries to nonzero values for desired solution options.
C
      do i = 1, 10
         infob(i) = 0
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
c     number of the derived functions
c
      nq       = 1              ! number of derived functions
c
c     Is the original system an index-2 system or not? 
c
      infob(5) = 1              ! 1 -- index-2 system
                                ! 0 -- not an index-2 system
c
c     do you want to compute the sensitivities with the initial
c     conditions inside DASPKADJOINT?
c
      infob(6) = 1              ! 1 -- yes
                                ! 0 -- no
c
c...  estimation of the total number of time steps in forward problem
      infob(7) = 400            ! 0 -- not know
                                ! nsteps
c
c     total number of the sentivity parameters
      info(19) = 8              ! number of the parameters
      info(20) = 0
      nbuf = 120*neq
      t = 0.0d0
c
c...  call cpu time 
      t1 = dtime(tt)
 20   IF(T .LT. TOUT) THEN 
*     
         if (info(20) .gt. 2) then
            CALL DDASPKadjoint(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, 
     *           res_adyp, nq, h_qres, qsen, infob, rtolb,atolb,nbuf,
     *           j_res, ieopt, res_adp, RES_ady, 
     *           g_res_ady, g_res_adyp, adinit, K_RES)
         else
            CALL DDASPKadjoint(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, 
     *           res_adyp, nq, qres, qsen, infob, rtolb,atolb,nbuf,
     *           j_res, ieopt, res_adp, RES_ady, 
     *           g_res_ady, g_res_adyp, adinit, K_RES)
         end if
         if (idid.gt.0) GOTO 20
      ENDIF 

      t2 = dtime(tt)
      write(*,*)
      write(*,*) 'Performance statistics:'
      write(*,*) ' Time it took:', tt(1)+tt(2)
      print *, ' maximum value of log(V(x)) = ', value(20,y)
c
c...  sensitivity about the initial condition
      print *, ' The sensitivities with the parameters are:'
      print *, (qsen(nq*(1+neq)+i),i=1,8)
      STOP
C------  End of main program for FEKETE_ADJOINT example program ------------------
      END

      double precision function value(n,v)
C
C...  This routine calculates the cost function for the optimization
C     problem (n is the number of Fekete points, v is the location)
C
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
      SUBROUTINE QRES (T,Y,YPRIME,QSEN,ires,RPAR,IPAR,senpar)
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
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
c
      subroutine res_adyp(t, y, yp, adyp, cj, delta, addelta, 
     &     ires, rpar, ipar, senpar)
C-------------------------------------------------------------------
C This is the user-supplied RES_ADP subroutine for this example.
C It computes a vector-matrix product and return it in ADYP, i.e.
C     adyp = addelta * F_y'
C--------------------------------------------------------------------
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     adyp(*), addelta(*), delta(*)
      integer ires, ipar(*), i
      
      do i = 1, 120
         adyp(i) = addelta(i)
      end do
      do i = 121, 160
         adyp(i) = 0.0d0
      end do

      return
      end
