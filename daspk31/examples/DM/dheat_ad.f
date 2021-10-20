C Copyright 2000 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  DHEAT
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900122   
C***REVISION DATE  920825   BC via algebraic eqns, minor revs.
C***REVISION DATE  920929   CJ in RES call sequence
C***REVISION DATE  950914   New names, minor revisions throughout.
C***REVISION DATE  990114   revise for sensitivity analysis
C***REVISION DATE  000714   revise for sensitivity analysis by adjoint method
C***REVISION DATE  070309   include resh.f directly instead of having 2 copies of res
C
C
C***AUTHORS  A. C. Hindmarsh, P. N. Brown
C            Lawrence Livermore National Laboratory
C            Livermore, CA 94551, USA
C
C            L. R. Petzold
C            University of California
C            Santa Barbara, CA  93106, USA
C
C            Shengtai Li
C            University of California
C            Santa Barbara, CA  93106, USA
C
C***DESCRIPTON
C
C-----------------------------------------------------------------------
C Example program for DDASPKADJOINT.
C DAE system derived from the discretized heat equation on a square.
C
C This is the double precision version.
C-----------------------------------------------------------------------
C
C This program solves a DAE system that arises from the heat equation,
C   du/dt = p_1*u   + p_2*u
C                xx        yy
C posed on the 2-D unit square with zero Dirichlet boundary conditions.
C An M+2 by M+2 mesh is set on the square, with uniform spacing 1/(M+1).
C The spatial deriviatives are represented by standard central finite
C difference approximations.  At each interior point of the mesh,
C the discretized PDE becomes an ODE for the discrete value of u.
C At each point on the boundary, we pose the equation u = 0.  The
C discrete values of u form a vector U, ordered first by x, then by y.
C The result is a DAE system G(t,U,U') = 0 of size NEQ = (M+2)*(M+2).
C
C Initial conditions are posed as u = 16x(1-x)y(1-y) at t = 0.
C The problem is solved by DDASPK on the time interval t .le. 10.24.
C
C p_1 and p_2 are the sensitivity parameters. In this example, we also 
C take the initial conditions as sensitivity parameters. We have two
C derived functions in this example
C    g(1) = sum_1^N (u_i*u_i)
C    g(2) = \int sum(u_i)
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   UINIT, DDASPKADJOINT, RES_ADY, RES_ADP, QRES
C
C***END PROLOGUE  DHEAT_AD
C
C     Here are necessary declarations.  The dimension statements use a
C     maximum value for the mesh parameter M, and assume ML = MU = 1.
C
      IMPLICIT NONE
      EXTERNAL RESH, QRES, RES_ADY, RES_ADP, DUMMY
      INTEGER MAXM, MAXM2, MXNEQ
      PARAMETER (MAXM = 50, MAXM2 = MAXM+2, MXNEQ=MAXM2*MAXM2)
      INTEGER NBUF, LRWF, MXNQ
      PARAMETER (NBUF = 24*MXNEQ, 
     *     LRWF = 50 + 9*MXNEQ + (2*(MAXM+3)+1)*MXNEQ,
     *     MXNQ = 2)
      INTEGER LIWF, MXNCHK
      PARAMETER (LIWF = 40 + 3*MXNEQ, MXNCHK=100)
      INTEGER LENRW
      PARAMETER (LENRW = 2*LRWF + 10*MXNEQ + 3*MXNEQ*MXNQ + NBUF)
      INTEGER LENIW, LENSP
      PARAMETER (LENIW = 2*LIWF + 4*MXNCHK)
      PARAMETER (LENSP = 2)
      DOUBLE PRECISION U, UPRIME, RWORK, RPAR, SENPAR
      INTEGER IWORK, INFO, IPAR
      DIMENSION U(MXNEQ),UPRIME(MXNEQ),RWORK(LENRW),IWORK(LENIW)
      DIMENSION INFO(30), RPAR(10), IPAR(10), SENPAR(LENSP)
      INTEGER LQSEN, IEOPT, INFOB
      PARAMETER (LQSEN = MXNQ*(MXNEQ+1+MXNEQ))
      DOUBLE PRECISION QSEN, ADY, ADYP
      dimension Qsen(LQSEN), ieopt(mxneq), infob(10)
      INTEGER  MAXNEQAD, NEQAD
      PARAMETER (MAXNEQAD=MXNQ*(MXNEQ+LENSP))
      real*4   t1, t2, t3, tt(2), secnds
      INTEGER M, NEQ, LRW, LIW, I, MBAND, NQ, IDID, LSEN, K
      DOUBLE PRECISION DX, COEFF, RTOL, ATOL, T, TOUT, PSOL
      DOUBLE PRECISION RTOLB, ATOLB
      DOUBLE PRECISION UMAX, HU
      INTEGER NQU, NST, NNI, NLI
      INTEGER TOTAL, LRWORKB, LIWORKB, LICONT
C
C     Here set parameters for the problem being solved.  Use RPAR and IPAR
C     to communicate these to the other routines.
C
      M = 40
      DX = 1.0D0/(M+1)
      NEQ = (M+2)*(M+2)
      COEFF = 1.0D0/(DX*DX)
C
      IPAR(3) = NEQ
      IPAR(4) = M
      RPAR(1) = DX
      RPAR(2) = COEFF
      SENPAR(1) = 1.0d0
      SENPAR(2) = 1.0d0
      LRW = LENRW
      LIW = LENIW
C
C-----------------------------------------------------------------------
C     Here we set up the INFO array, which describes the various options
C     in the way we want DDASPK to solve the problem.
C
C     We first initialize the entire INFO array to zero, then set select
C     entries to nonzero values for desired solution options.
C-----------------------------------------------------------------------
C
      DO 10 I = 1,30
 10     INFO(I) = 0
      INFO(12) = 0        ! direct method
      info(5) = 0
C     Jacobian is banded
      info(6) = 1      
      mband = m+3
      iwork(1) = M+3
      iwork(2) = M+3
      iwork(38) = 0
      ipar(5) = 0
C
C     Call subroutine UINIT to initialize U and UPRIME.
C
      CALL UINIT (U, UPRIME, RPAR, IPAR)
C
C     Here we set tolerances for DDASPK to indicate how much accuracy 
C     we want in the solution, in the sense of local error control.
C     For this example, we ask for pure absolute error control with a
C     tolerance of 1.0D-3.
      RTOL = 1.0d-5
      ATOL = 1.0D-5
      RTOLB = RTOL
      ATOLB = ATOL
C
      t1 = dtime(tt)
      t1 = secnds(0.0e0)
c---------------------------------------------------------------
c  sensitivity analysis
c---------------------------------------------------------------
c
c...  IEOPT, used when INFOB(2) != 3.
c
      do i = 1, neq
         ieopt (i) = 1
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
C
c
c...  total number of derived functions
      nq    = 2
c
c...  total number of sensitivity parameters
      info(19) = 200            
c
c...  number of parameters that appear in the RES routine
      info(22) = 2
c
c...  number of the quadrature variables in the forward integration
      info(28) = 1
      neq = neq + info(28)
      t = 0.0d0
      tout = 0.16d0 
      
      IF (INFOB(4) .EQ. 0 .AND. INFO(22) .GT. 0) THEN
          NEQAD = NQ * (NEQ + INFO(22))
      ELSE
          NEQAD = NQ * NEQ
      END IF

      CALL DDASPKadjoint(
     1     RESH, NEQ, T, U, UPRIME, TOUT, TOUT, INFO, RTOL, ATOL,
     1     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, DUMMY,
     1     PSOL,SENPAR, 
     1     DUMMY, NQ, QRES, QSEN, INFOB, RTOLB, ATOLB, NBUF,
     1     DUMMY, IEOPT, RES_ADP, RES_ADY,
     1     NEQAD, ADY, ADYP)
      IF (IDID .LT. 0) THEN
         WRITE(*,65)T, idid
 65      FORMAT(//' Final time reached =',E12.4, ', idid =', I5 //)
         GO TO 80
      ENDIF
      
      UMAX = 0.0D0
      DO 50 I = 1,IPAR(3)
 50      UMAX = MAX (UMAX, ABS(U(I)) )
C
      HU = RWORK(7)
      NQU = IWORK(8)
      NST = IWORK(11)
      NNI = IWORK(19)
      NLI = IWORK(20)
      WRITE (6,60) T,UMAX,NQU,HU,NST,NNI,NLI
 60   FORMAT(E15.5,E12.4,I5,E14.3,I7,I9,I8)
         
C
C     Here we display some final statistics for the problem.
C
 80   CONTINUE
      t2 = dtime(tt)
      write(*,*)
      write(*,*) 'Performance statistics:'
      write(*,*) ' Time it took:', tt(1)+tt(2)
      TOTAL = IWORK(18) ! length of first chunk of RWORK (RWORKF)
      LICONT  = IWORK(39)
      LIWORKB = IWORK(LICONT+1)
      LRWORKB = IWORK(40)
      IF (LRWORKB.NE.1) THEN
        TOTAL = TOTAL + IWORK(LIWORKB+18) ! length of RWORKB
      END IF
      TOTAL = TOTAL + LRW - IWORK(LICONT+18) + 1 ! length of buffer
      print *, ' RWORK size = ', TOTAL
      TOTAL = IWORK(17) ! length of first section of IWORK (IWORKF)
      TOTAL = TOTAL + 90 ! length of ICONT
      IF (LIWORKB .NE. 1) THEN
        TOTAL = TOTAL + IWORK(LIWORKB+17) ! len of IWORKB
      END IF
      print *, ' IWORK size = ', TOTAL
c
c     Here we display the values for the derived functions and their
c     sensitivity with respect to the sensitivity parameters p
c
      print *, ' value of sum(u_i^2) and integral sum(u_i) at t=', tout
      print *, '        ', (qsen(i), i = 1, nq)
      print *, '        ', (qsen(i), i = nq+1, nq+2)
      write (14,*), (qsen(i), i=1,neqad)
      lsen = nq*(neq+1)+1
      print *, ((qsen(lsen+i*nq+k), i=0,1),k=0,nq-1)
c
c     From lsen+2*nquad, qsen(*) stores the sensitivity with
c     respect to the initial conditions.
c
c      do i = 13, 27
c         write(12,*) qsen(1+20*42 + i)
c      end do
C
C------  End of main program for DHEAT_AD example program --------------
      stop
      END

      SUBROUTINE UINIT (U, UPRIME, RPAR, IPAR)
C
C This routine computes and loads the vector of initial values.
C The initial U values are given by the polynomial u = 16x(1-x)y(1-y).
C The initial UPRIME values are set to zero.  (DDASPK corrects these
C during the first time step.)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), RPAR(*), IPAR(*)
C
      NEQ = IPAR(3)
      M = IPAR(4)
      M2 = M + 2
      DX = RPAR(1)
      COEFF = RPAR(2)
C
      do i = 1, neq+1
         u(i) = 0.0d0
         uprime(i) = 0.0d0
      end do
      
      DO 20 K = 1,M
        YK = K*DX
        IOFF = M2*K
        DO 10 J = 1,M
          XJ = J*DX
          I = IOFF + J + 1
          U(I) = 16.0D0*XJ*(1.0D0-XJ)*YK*(1.0D0-YK)
 10    CONTINUE
 20   CONTINUE
c ... unforce the consistent initial condition
      RETURN
C------------  End of Subroutine UINIT  --------------------------------
      END

C------------  Subroutine QRES  ---------------------------
      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
C
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
C
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i
C
C     Set problem constants using IPAR and RPAR.
      neq = ipar(3)
      M = IPAR(4)
      M2 = M + 2
C
C     define the value for the derived functions
C
      qsen(1) = 0.0
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
           I = IOFF + J + 1
           qsen(1) = qsen(1) + y(i)*y(i)
 20     CONTINUE
 30   CONTINUE
!       print *, 'qres called'
      qsen(2) = y(neq+1)
      return
      end
C -------------------------------------------------------------------
      subroutine test_res_ady(T,Y,YP,RPAR,IPAR,SENPAR)
      implicit none
      double precision T, Y(2704),YP(2704),RPAR(*),SENPAR(*)
      integer ipar(*), i, ires
      double precision CJ, DELTA(1766), ADDELTA(1766), ADY(2704)
      CJ = 1.0
      ires = 3
      do i = 1, 1766
        delta(i) = 1.0d0
        addelta(i) = 1.0d0
      end do
      do i = 1, 2704
        ady(i) = 1.0d0
      end do
      call res_ady(
     *     T, Y, ADY, YP, CJ, delta, ADDELTA,IRES, RPAR, IPAR, SENPAR)
      return
      end
C -------------------------------------------------------------------

C
C---- interface routine to the Tapenade-generated code
C
      SUBROUTINE RES_ADY(T, Y, ADY, YPRIME, CJ, DELTA, ADDELTA,
     *        IRES, RP, IP, SENPAR)
C
C This is the user-supplied RES_ADY subroutine for this example.
C It computes a vector-matrix product and return it in ady, i.e.
C     ady = addelta * F_y
C
      IMPLICIT NONE
      DOUBLE PRECISION T, Y, ADY, YPRIME, DELTA, RP, ADDELTA
      DOUBLE PRECISION SENPAR,CJ
      INTEGER IP,IRES,II
      DIMENSION ADY(*),YPRIME(*),DELTA(*),Y(*),RP(*),ADDELTA(*)
      DIMENSION IP(*)
      DOUBLE PRECISION DELTAB(1766), DELTA_COPY(1766), ADY_SAV(2704)
          
      DO II = 1, 1766
        DELTAB(II) = ADDELTA(II)
        DELTA_COPY(II) = DELTA(II)
      END DO
      DO ii = 1, 2704
        ADY_SAV(II) = ADY(II)
      end do
      
      CALL RESH_B(t, y, ady, yprime, cj, delta_copy, deltab,
     +           ires, rp, ip, senpar)     
     
      DO II = 1, 2704
        ADY(II) = ADY(II) + ADY_SAV(II)
      END DO

      return
      end
         
      SUBROUTINE RES_ADP(
     *             T, Y, YPRIME, CJ, DELTA, ADDELTA, IRES, RP, 
     *             IP, SENPAR, ADSENPAR)
C
C This is the user-supplied RES_ADP subroutine for this example.
C It computes a vector-matrix product and return it in ADSENPAR, i.e.
C     adsenpar = addelta * F_p
C
      IMPLICIT NONE
      DOUBLE PRECISION T, Y, YPRIME, DELTA, RP, ADDELTA
      DOUBLE PRECISION SENPAR,CJ,ADSENPAR(*)
      INTEGER IP,IRES,II
      DIMENSION YPRIME(*),DELTA(*),Y(*),RP(*),ADDELTA(*)
      DIMENSION IP(*)
      DOUBLE PRECISION DELTAB(1766), ADSENPAR_SAV(2)

      DO II = 1, 1766
        DELTAB(II) = ADDELTA(II)
      END DO
      DO II = 1, 2
        ADSENPAR_SAV(II) = ADSENPAR(II)
      END DO
      CALL RESH_OFP_B(t, y, yprime, cj, delta, deltab, ires, rp, 
     +               ip, senpar, adsenpar)
      DO II = 1, 2
        ADSENPAR(II) = ADSENPAR(II) + ADSENPAR_SAV(II)
      END DO
      return
      end

      subroutine DUMMY()
      return
      end
      