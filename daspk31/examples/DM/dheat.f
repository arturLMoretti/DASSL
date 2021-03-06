C Work performed under the auspices of the U.S. Department of Energy
C by Lawrence Livermore National Laboratory under contract number 
C W-7405-Eng-48.
C
C Copyright 1995 the Regents of the University of California.
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
C***REVISION DATE  070116   added sensitivity of derived function
C
C***DESCRIPTON
C
C-----------------------------------------------------------------------
C Example program for DDASPK.
C DAE system derived from the discretized heat equation on a square.
C
C This is the double precision version.
C-----------------------------------------------------------------------
C
C This program solves a DAE system that arises from the heat equation,
C   du/dt = u   + u
C            xx    yy
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
C The Krylov linear system solution method, with preconditioning, is
C selected.  The preconditioner is a band matrix with half-bandwidths
C equal to 1, i.e. a tridiagonal matrix.  (The true half-bandwidths
C are equal to M+2.)  This corresponds to ignoring the y-direction
C coupling in the ODEs, for purposes of preconditioning.  The extra
C iterations resulting from this approximation are offset by the lower
C storage and linear system solution costs for a tridiagonal matrix.  
C
C The routines DBANJA and DBANPS that generate and solve the banded
C preconditioner are provided in a separate file for general use.
C
C The output times are t = .01 * 2**n (n = 0,...,10).  The maximum of
C abs(u) over the mesh, and various performance statistics, are printed.
C
C For details and test results on this problem, see the reference:
C   Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
C   Using Krylov Methods in the Solution of Large-Scale Differential-
C   Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   UINIT, DDASPK
C
C***END PROLOGUE  DHEAT
C
C Here are necessary declarations.  The dimension statements use a
C maximum value for the mesh parameter M, and assume ML = MU = 1.
C
      IMPLICIT NONE
      EXTERNAL RESH, DBANJA, DBANPS
      EXTERNAL QRES
      INTEGER  Npara, NQ, LQSEN
      parameter (Npara=20,NQ=1,LQSEN=NQ*(Npara+1))
      INTEGER  MAXM, MAXM2, MXNEQ
      PARAMETER (MAXM = 50, MAXM2 = MAXM+2, MXNEQ=(npara+6)*MAXM2*MAXM2)
*      PARAMETER (LENRW = 91 + 18*MXNEQ, LENIW  = 40+MXNEQ,
      INTEGER LENRW, LENIW, LENWP, LENIWP
      PARAMETER (LENRW = 91 + 18*MXNEQ + 2*NQ, LENIW  = 40+MXNEQ,
     *          LENWP  = 5*MXNEQ + 2*((MXNEQ/3) + 1), LENIWP = MXNEQ)
      DOUBLE PRECISION U, UPRIME, RWORK, RPAR, SENPAR
      INTEGER IWORK, INFO, IPAR
      DIMENSION U(MXNEQ),UPRIME(MXNEQ),
     *          RWORK(LENRW+LENWP),IWORK(LENIW+LENIWP)
      DIMENSION INFO(30), RPAR(10), IPAR(10), SENPAR(2)
*      DIMENSION INFO(30), RPAR(2), IPAR(4)
      DOUBLE PRECISION QSEN, RTOL, ATOL, T, TOUT
      DIMENSION QSEN(LQSEN)
      INTEGER M, NEQ, ML, MU, MBAND, MSAVE, LENPD
      INTEGER LWP, LIWP, LRW, LIW, I, INSO, LOUT, NOUT
      INTEGER IOUT, IDID, NQU, NST, NNI, NLI, NPE, NRE, NPS, NSE
      INTEGER NETF, NCFN, NCFL, IRES
      DOUBLE PRECISION DX, COEFF, G_RES, UMAX, AVDIM, HU
      real*4   t1, t2, t3, tt(2), secnds
      common /seninfo/info, neq
      common /mband/mband
c
c === CPU time set-up
c
C
C Here set parameters for the problem being solved.  Use RPAR and IPAR
C to communicate these to the other routines.
C
c      M = MAXM
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
C
C Here set the half-bandwidths and load them into IPAR for use by the
C preconditioner routines.
      ML = 1
      MU = ML
      IPAR(1) = ML
      IPAR(2) = MU
C
C Here set the lengths of the preconditioner work arrays WP and IWP,
C load them into IWORK, and set the total lengths of WORK and IWORK.
      LENPD = (2*ML + MU + 1)*NEQ
      MBAND = ML + MU + 1
      MSAVE = (NEQ/MBAND) + 1
      LWP = LENPD + 2*MSAVE
      LIWP = NEQ
      IWORK(27) = LWP
      IWORK(28) = LIWP
      LRW = LENRW + LWP
      LIW = LENIW + LIWP
C
C-----------------------------------------------------------------------
C Here we set up the INFO array, which describes the various options
C in the way we want DDASPK to solve the problem.
C In this case, we select the iterative preconditioned Krylov method,
C and we supply the band preconditioner routines DBANJA/DBANPS.
C
C We first initialize the entire INFO array to zero, then set select
C entries to nonzero values for desired solution options.
C
C To select the Krylov iterative method for the linear systems,
C we set INFO(12) = 1.
C
C Since we are using a preconditioner that involves approximate
C Jacobian elements requiring preprocessing, we have a JAC routine,
C namely subroutine DBANJA, and we must set INFO(15) = 1 to indicate
C this to DDASPK.
C
C No other entries of INFO need to be changed for this example.
C-----------------------------------------------------------------------
C
      DO 10 I = 1,30
 10     INFO(I) = 0
      INFO(12) = 0        ! direct method
      if (info(12) .eq. 0) then
         info(6) = 1
         info(5) = 0
         mband = m+3
         iwork(1) = M+3
         iwork(2) = M+3
      end if
      iwork(38) = 0
      INFO(15) = 1
c---------------------------------------------------------------
c  sensitivity analysis
c---------------------------------------------------------------
      inso = 18
      info(inso+1) = npara
      print *, ' number of sensitivities? (0..20) '
      read(*,*) info(inso+1)
      neq = neq*(info(inso+1) + 1)
      ipar(5) = info(inso+1)
      if (info(inso+1) .eq. 0) go to 25
      info(inso+2) = 2          ! 0, central, 1, forward, 2 supply
      print *, ' method for computing sensitivity residues:(0,1,2) '
      print *, '    0: central differencing '
      print *, '    1: one-side forward differencing '
      print *, '    2: analytic input'
      read (*,*) info(inso+2)
      info(inso+3) = 0          ! default perturbed factor
      info(inso+4) = 2          ! two parameter appears in the RES
      info(inso+5) = 0          ! 0, full control 1. partial control
      print *, ' Error control: (0,1)'
      print *, '    0: full error control with sensitivity '
      print *, '    1: partial error control on state variables only '
      read(*, *) info(inso+5)
      if (info(20) .NE. 2) then
          info(inso+6) = NQ         ! two derived functions
      else
          info(inso+6) = 0          ! no derived information
      end if
      info(inso+7) = 1          ! simultaneous corrector(0), staggered(1)
      print *, ' Corrector method: (0, 1)'
      print *, '    0:  simultaneous corrector method '
      print *, '    1:  staggered corrector method'
      read(*,*) info(inso+7)
      info(inso+9) = 0
 25   continue
C
C Call subroutine UINIT to initialize U and UPRIME.
C
      t1 = dtime(tt)
      t1 = secnds(0.0e0)
      CALL UINIT (U, UPRIME, RPAR, IPAR)
C
C Here we set tolerances for DDASPK to indicate how much accuracy 
C we want in the solution, in the sense of local error control.
C For this example, we ask for pure absolute error control with a
C tolerance of 1.0D-3.
      RTOL = 1.0d-4
      ATOL = 1.0D-4
C
C Here we generate a heading with important parameter values.
C LOUT is the unit number of the output device.
      LOUT = 6
      WRITE (LOUT,30) M,NEQ,INFO(12),ML,MU,RTOL,ATOL
 30   FORMAT(' Heat Equation Example Program for DDASPK'//
     1       '    M+2 by M+2 mesh, M =',I3,',  System size NEQ =',I5//
     1       '    Linear solver method flag INFO(12) =',I3,
     1       '    (0 = direct, 1 = Krylov)'/
     1       '    Preconditioner is a banded approximation with ML =',
     1            I3,'  MU =',I3//
     1       '    Tolerances are RTOL =',E10.1,'   ATOL =',E10.1//)
      WRITE (LOUT,40)
 40   FORMAT(5X,'t',12X,'UMAX',8X,'NQ',8X,'H',8X,'STEPS',
     1       5X,'NNI',5X,'NLI', 2x, 'Time'/)
C
C-----------------------------------------------------------------------
C Now we solve the problem.
C
C DDASPK will be called to compute 11 intermediate solutions from
C tout = 0.01 to tout = 10.24 by powers of 2.
C
C We pass to DDASPK the names DBANJA and DBANPS for the JAC and PSOL
C routines to do the preconditioning.
C
C At each output time, we compute and print the max-norm of the
C solution (which should decay exponentially in t).  We also print
C some relevant statistics -- the current method order and step size,
C the number of time steps so far, and the numbers of nonlinear and
C linear iterations so far.
C
C If DDASPK failed in any way (IDID .lt. 0) we print a message and
C stop the integration.
C-----------------------------------------------------------------------
C
      NOUT = 11
      T = 0.0D0
      TOUT = 0.01D0
      DO 70 IOUT = 1,NOUT
         CALL DDASPK (
     1        RESH, NEQ, T, U, UPRIME, TOUT, INFO, RTOL, ATOL,
     1        IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, DBANJA,
     1        DBANPS,SENPAR, G_RES)
         UMAX = 0.0D0
         DO 50 I = 1,NEQ/(info(inso+1) + 1)
 50        UMAX = MAX (UMAX, ABS(U(I)) )
C
         HU = RWORK(7)
         NQU = IWORK(8)
         NST = IWORK(11)
         NNI = IWORK(19)
         NLI = IWORK(20)
         t2 = secnds(t1)
         WRITE (LOUT,60) T,UMAX,NQU,HU,NST,NNI,NLI, t2
 60      FORMAT(E15.5,E12.4,I5,E14.3,I7,I9,I8, f7.2)
C

         IF (TOUT .GE. 0.16) THEN
           PRINT *, 'Finished at time ', TOUT
           GO TO 80
         END IF 
         
         IF (IDID .LT. 0) THEN
           WRITE(LOUT,65)T
 65        FORMAT(//' Final time reached =',E12.4//)
           GO TO 80
           ENDIF
C
         TOUT = TOUT*2.0D0
 70      CONTINUE
 998  format(3i6, 2x, f13.8)
 999  format(3f15.9)
C
C Here we display some final statistics for the problem.
C The ratio of NLI to NNI is the average dimension of the Krylov
C subspace involved in the Krylov linear iterative method.
 80   CONTINUE
      NST = IWORK(11)
      NPE = IWORK(13)
      NRE = IWORK(12)
      LIW = IWORK(17)
      LRW = IWORK(18)
      NNI = IWORK(19)
      NLI = IWORK(20)
      NPS = IWORK(21)
      NSE = iwork(22)
      IF (NNI .NE. 0) AVDIM = REAL(NLI)/REAL(NNI)
      NETF = IWORK(14)
      NCFN = IWORK(15)
      NCFL = IWORK(16)
      WRITE (LOUT,90) LRW,LIW,NST,NRE,NSE,
     $     NPE,NPS,NNI,NLI,AVDIM,NETF,NCFN,NCFL
 90   FORMAT(//' Final statistics for this run..'/
     1   '   RWORK size =',I6,'   IWORK size =',I4/
     1   '   Number of time steps ................ =',I5/
     1   '   Number of residual evaluations ...... =',I5/
     3   '   Number of sensitivity evaluations ....=',I5/
     1   '   Number of preconditioner evaluations  =',I5/
     1   '   Number of preconditioner solves ..... =',I5/
     1   '   Number of nonlinear iterations ...... =',I5/
     1   '   Number of linear iterations ......... =',I5/
     1   '   Average Krylov subspace dimension =',F8.4/
     1   '   error test failure .............. =',I5/
     1   I5,' nonlinear conv. failures,',I5,' linear conv. failures')
      t2 = dtime(tt)
      t3 = tt(1)
c      write(*,*) t2, tt(1), tt(2)
c      t2 = secnds(t1)
c      t3 = t2
      write(6,*) '   Time it took:',t2, ', System time:', t3

      IF (INFO(19) .GT. 0 .AND. INFO(20) .NE. 2) THEN
        CALL DSENSD(
     *       QRES, NEQ, T, U, UPRIME, QSEN, INFO, 
     *       RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)
        print *, '\n \\sum y_i**2:', QSEN(1)
        print *, 'sensitivity of \\sum y_i**2 to p1:', QSEN(2)
      END IF
      
      STOP
      END
C
C------  End of main program for DHEAT example program -----------------
       subroutine qres(T, Y, YP, QSEN, IRES, RPAR, IPAR, SENPAR)
C
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
C
      implicit none
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i, neq, m, m2

      neq = ipar(3)
C
C     define the value for the derived functions
C
      qsen(1) = 0.0d0
c      qsen(2) = 0.0d0
      DO 30 I = 1, NEQ
        qsen(1) = qsen(1) + y(i)*y(i)
c        qsen(2) = qsen(2) + y(i) ! This needs to be a quadrature variable
 30   CONTINUE
      return
      end
c -------------------------------------------------------------


      SUBROUTINE UINIT (U, UPRIME, RPAR, IPAR)
C
C This routine computes and loads the vector of initial values.
C The initial U values are given by the polynomial u = 16x(1-x)y(1-y).
C The initial UPRIME values are set to zero.  (DDASPK corrects these
C during the first time step.)
C
      IMPLICIT NONE
      DOUBLE PRECISION U, UPRIME, RPAR, DX, COEFF, YK, XJ
      INTEGER IPAR, NEQ, NJP, M, M2, K, J, I, IOFF, NP
      DIMENSION U(*), UPRIME(*), RPAR(*), IPAR(*)
C
      NEQ = IPAR(3)
      np  = ipar(5)
      M = IPAR(4)
      M2 = M + 2
      DX = RPAR(1)
      COEFF = RPAR(2)
C
      do i = 1, (np+1)*neq
         u(i) = 0.0d0
         uprime(i) = 0.0d0
      end do
      
      do i = 4, np+1
         u((i-1)*neq + M2 + i) = 1.0d0
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

      SUBROUTINE RESH1 (T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, 
     *                  SENPAR)
C
C This is the user-supplied RES subroutine for this example.
C It computes the residuals for the 2-D discretized heat equation,
C with zero boundary values.
C
      IMPLICIT NONE
      DOUBLE PRECISION T, U, UPRIME, CJ, DELTA, RPAR, SENPAR, COEFF
      INTEGER IPAR, NEQ, M, M2, I, K, J, IOFF, IRES
      DOUBLE PRECISION TEMX, TEMY
      DIMENSION U(*), UPRIME(*), DELTA(*), RPAR(*), IPAR(*),
     *          SENPAR(*)
C
C Set problem constants using IPAR and RPAR.
      NEQ = IPAR(3)
      M = IPAR(4)
      COEFF = RPAR(2)
      M2 = M + 2
C
C Load U into DELTA, in order to set boundary values.
      DO 10 I = 1,NEQ
 10     DELTA(I) = U(I)
C
C Loop over interior points, and load residual values.
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
          I = IOFF + J + 1
          TEMX = SENPAR(1)*(U(I-1)  + U(I+1)  - 2.0d0*u(i))
          TEMY = SENPAR(2)*(U(I-M2) + U(I+M2) - 2.0d0*u(i))
          DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF
 20       CONTINUE
 30     CONTINUE
C
      RETURN
C------------  End of Subroutine RESH  ---------------------------------
      END
      subroutine ressa(t,u,uprime,cj,delta,ires,rpar,ipar,SENPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),UPRIME(*),DELTA(*),RPAR(*),IPAR(*),SENPAR(*)
C
C Set problem constants using IPAR and RPAR.
      NEQ = IPAR(3)
      M = IPAR(4)
      np  = ipar(5)
      COEFF = RPAR(2)
      M2 = M + 2
C
C Load U into DELTA, in order to set boundary values.
      DO 10 I = 1,(np+1)*NEQ
 10     DELTA(I) = U(I)
C
C Loop over interior points, and load residual values.
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
           I = IOFF + J + 1
           TEMX = (U(I-1)  + U(I+1)  - 2.0d0*u(i))
           TEMY = (U(I-M2) + U(I+M2) - 2.0d0*u(i))
           Delta(i) = UPRIME(I) - 
     &          (SENPAR(1)*TEMX + SENPAR(2)*TEMY)*COEFF
           do ip = 1, np
              i = i + neq
              TEMXp = (U(I-1)  + U(I+1)  - 2.0d0*u(i))
              TEMYp = (U(I-M2) + U(I+M2) - 2.0d0*u(i))
              Delta(i) = UPRIME(I) -  (SENPAR(1)*TEMXp + 
     &             SENPAR(2)*TEMYp)*COEFF
              if (ip .eq. 1) delta(i) = delta(i) - temx*coeff
              if (ip .eq. 2) delta(i) = delta(i) - temy*coeff
           end do
 20     CONTINUE
 30   CONTINUE
      
      return
      end
      SUBROUTINE RESH (T,U,UPRIME,CJ,DELTA,IRES,RPAR,IPAR,SENPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),UPRIME(*),DELTA(*),RPAR(*),IPAR(*),SENPAR(*)

      if (ires .eq. 0) then
         call RESH1 (
     *        T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
      else if (ires. eq. 1) then
         call RESSA (
     *        T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
      else 
         ires = -2
      end if
      return
      end
