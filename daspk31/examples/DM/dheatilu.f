C Work performed under the auspices of the U.S. Department of Energy
C by Lawrence Livermore National Laboratory under contract number 
C W-7405-Eng-48.
C
C Copyright 1995 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  DHEATILU
C***REFER TO  DDASPK
C***DATE WRITTEN   951005   (YYMMDD)
C***REVISION DATE  951005   Reworking of heat example using banded J
C***REVISION DATE  990117   Sensitivity analysis for parameters
C***REVISION DATE  070116   Moved QRES into this file from reshilu.f
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
C selected.  The preconditioner is a sparse matrix with half-bandwidths
C equal to 1, i.e. a tridiagonal matrix.  (The true half-bandwidths
C are equal to M+2.)  This corresponds to ignoring the y-direction
C coupling in the ODEs, for purposes of preconditioning.  The extra
C iterations resulting from this approximation are offset by the lower
C storage and linear system solution costs for a tridiagonal matrix.  
C
C The routines DJACILU and DPSOLILU that generate and solve the sparse
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
C***END PROLOGUE  DHEATILU
C
C Here are necessary declarations.  The dimension statements use a
C maximum value for the mesh parameter M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RESH, DJACILU, DPSOLILU, QRES
C Load parameters for sparse preconditioner.
      PARAMETER (LENPFAC = 5)   ! the average number of non-zeros
      PARAMETER (LENPLUFAC = 5) 
      PARAMETER (IPREMETH = 1)  ! =1 means ILUT preconditioner used
                                ! =2 means ILUTP preconditioner used
                                ! =3 ILUD
                                ! =4 ILUDP
      PARAMETER (LFILILUT = 5)
      PARAMETER (IREORDER = 1)
      PARAMETER (ISRNORM = 0)
      PARAMETER (NORMTYPE = 2)
      PARAMETER (JACOUT = 0)
      PARAMETER (JSCALCOL = 0)
      PARAMETER (TOLILUT = 0.001)
      PARAMETER (PERMTOL = 0.01)
C Load workspace lengths.
      parameter (Npara=20)
      PARAMETER (MAXM = 50, MAXM2 = MAXM+2, MXNEQ=(npara+6)*MAXM2*MAXM2)
      PARAMETER (MXNY = MAXM2*MAXM2)
      PARAMETER (LENWP = 2*LENPFAC*MXNY + LENPLUFAC*MXNY
     .                   + ISRNORM*MXNY + 2*(MXNY+1) )
C      PARAMETER (LENWP = 2*LENPFAC*MXNEQ + LENPLUFAC*MXNEQ
C     .                   + ISRNORM*MXNEQ + 2*(MXNEQ+1) )
C      PARAMETER (LENIWP = 4*(MXNEQ+1) + 3*LENPFAC*MXNEQ
C     .                    + 2*LENPLUFAC*MXNEQ + IREORDER*2*MXNEQ
C     .                    + (IPREMETH-1)*2*MXNEQ )
      PARAMETER (LENIWP = 4*(MXNY+1) + 3*LENPFAC*MXNY
     .                    + 2*LENPLUFAC*MXNY + IREORDER*2*MXNY
     .                    + (IPREMETH-1)*2*MXNY )
      PARAMETER (LENRW = 91 + 18*MXNEQ, LENIW  = 40 + mxny)
      DIMENSION U(MXNEQ),UPRIME(MXNEQ),
     *          RWORK(LENRW+LENWP),IWORK(LENIW+LENIWP)
*      DIMENSION INFO(20), RPAR(4), IPAR(34)
      DIMENSION INFO(30), RPAR(4), IPAR(35), SENPAR(2), QSEN(100)
      real*4   t1, t2, t3, tt(2),secnds
c
c === CPU time set-up
c
C Set LOUT, the unit number of the output device.
      LOUT = 6

C Open output file if JACOUT .EQ. 1.
      IF (JACOUT .EQ. 1) THEN
         IPAR(29) = 1
         OPEN(UNIT=1, FILE='Heat_Test_Matrix.dat', STATUS='unknown')
      ENDIF
C
C Here set parameters for the problem being solved.  Use RPAR and IPAR
C to communicate these to the other routines.
C
      M = 40
c      M = MAXM
      DX = 1.0D0/(M+1)
      NEQ = (M+2)*(M+2)
      COEFF = 1.0D0/(DX*DX)
      IPAR(33) = NEQ
      IPAR(34) = M
      RPAR(3) = DX
      RPAR(4) = COEFF
      SENPAR(1) = 1.0d0
      SENPAR(2) = 1.0d0
C
C Here set the lengths of the preconditioner work arrays WP and IWP,
C load them into IWORK, and set the total lengths of WORK and IWORK.
      IWORK(27) = LENWP
      IWORK(28) = LENIWP
      LRW = LENRW + LENWP
      LIW = LENIW + LENIWP
C Load values into IPAR and RPAR for sparse preconditioner.
      ML = M+3
      MU = ML
      IPAR(1) = ML
      IPAR(2) = MU
      MBAND = ML + MU + 1
      IPAR(3) = LENPFAC
      IPAR(4) = LENPLUFAC
      IPAR(5) = IPREMETH
      IPAR(6) = LFILILUT
      IPAR(7) = IREORDER
      IPAR(8) = ISRNORM
      IPAR(9) = NORMTYPE
      IPAR(10) = JACOUT
      IPAR(11) = JSCALCOL
      IPAR(12) = 0          ! Jacobian evaluation method
      IPAR(13) = 0          ! data dependence for rpar(*)
      IPAR(30) = 0
      RPAR(1) = TOLILUT
      RPAR(2) = PERMTOL
C Check IPAR, RPAR, LENWP and LENIWP for illegal entries and long
C enough work array lengths.
      CALL DSPSETUP (NEQ, lenwp, leniwp, RPAR, IPAR, IERR,
     .               LWPMIN, LIWPMIN)
      IWORK(27) = LWPMIN
      IWORK(28) = LIWPMIN
      LRW = LENRW + IWORK(27)
      LIW = LENIW + IWORK(28)
      IF (IERR .NE. 0) THEN
         WRITE(LOUT,15) IERR
 15      FORMAT(' Error return from DSPSETUP: IERR = ',i5)
         IF (LWPMIN .GT. LENWP) THEN
            WRITE(LOUT,*) ' More WP work array length needed'
         ENDIF
         IF (LIWPMIN .GT. LENIWP) THEN
            WRITE(LOUT,*) ' More IWP work array length needed'
         ENDIF
         STOP
      ENDIF
C
C-----------------------------------------------------------------------
C Here we set up the INFO array, which describes the various options
C in the way we want DDASPK to solve the problem.
C In this case, we select the iterative preconditioned Krylov method,
C and we supply the sparse preconditioner routines DJACILU/DPSOLILU.
C
C We first initialize the entire INFO array to zero, then set select
C entries to nonzero values for desired solution options.
C
C To select the Krylov iterative method for the linear systems,
C we set INFO(12) = 1.
C
C Since we are using a preconditioner that involves approximate
C Jacobian elements requiring preprocessing, we have a JAC routine,
C namely subroutine DJACILU, and we must set INFO(15) = 1 to indicate
C this to DDASPK.
C
C No other entries of INFO need to be changed for this example.
C-----------------------------------------------------------------------
C
      DO 20 I = 1,30
 20     INFO(I) = 0
C
      INFO(12) = 1
      INFO(15) = 1
c---------------------------------------------------------------
c  sensitivity analysis
c---------------------------------------------------------------
      inso = 18
      info(inso+1) = npara
      print *, ' number of sensitivities? (0..20) '
      read(*,*) info(inso+1)
      neq = neq*(info(inso+1) + 1)
      ipar(35) = info(inso+1)
      if (info(inso+1) .eq. 0) go to 25
      info(inso+2) = 2          ! 0, central, 1, forward, 2 supply
      print *, ' method for computing sensitivity residues:(0,1,2) '
      print *, '    0: central differencing '
      print *, '    1: one-side forward differencing '
      print *, '    2: user supply '
      read (*,*) info(inso+2)
      info(inso+3) = 0          ! default perturbed factor
      info(inso+4) = 2          ! two parameter appears in the RES
      info(inso+5) = 0          ! 0, full control 1. partial control
      print *, ' Error control: (0,1)'
      print *, '    0: full error control with sensitivity '
      print *, '    1: partial error control on state variables only '
      read(*, *) info(inso+5)
      info(inso+6) = 1          ! one derived information
      info(inso+7) = 1          ! simultaneous corrector(0), staggered(1)
      print *, ' Corrector method: (0, 1)'
      print *, '    0:  simultaneous corrector method '
      print *, '    1:  staggered corrector method'
      read(*,*) info(inso+7)
 25   continue
      iwork(38) = 0
C
C Call subroutine UINIT to initialize U and UPRIME.
C
      CALL UINIT (U, UPRIME, RPAR, IPAR)
C
C Here we set tolerances for DDASPK to indicate how much accuracy 
C we want in the solution, in the sense of local error control.
C For this example, we ask for pure absolute error control with a
C tolerance of 1.0D-3.
      RTOL = 1.0d-5
      ATOL = 1.0D-5
C
C Here we generate a heading with important parameter values.
C
      WRITE (LOUT,30) M,NEQ,INFO(12),ML,MU,IPREMETH,RTOL,ATOL
 30   FORMAT(' Heat Equation for DDASPK'//
     1       '    M+2 by M+2 mesh, M =',I3,',  System size NEQ =',I4//
     1       '    Linear solver method flag INFO(12) =',I3,
     1       '    (0 = direct, 1 = Krylov)'/
     1       '    Preconditioner is a sparse approximation with ML =',
     1            I3,'  MU =',I3/
     1       '    Incomplete factorization option =',I2,
     1       '    (1=ILUT,  2=ILUTP, 3=ILUD,'/
     1       '                                       ',
     1       '     4=ILUDP, 5=ILU0,  6=MILU0)'/,
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
C We pass to DDASPK the names DJACILU and DPSOLILU for the JAC and PSOL
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
      NOUT = 5
      T = 0.0D0
      TOUT = 0.01D0
      t1 = secnds(1.0)
      DO 70 IOUT = 1,NOUT
 55      continue
         CALL DDASPK (
     1        RESH, NEQ, T, U, UPRIME, TOUT, INFO, RTOL, ATOL,
     1        IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, 
     1        DJACILU, DPSOLILU, SENPAR, H_RES)
         CALL DSENSD(
     *        QRES, NEQ, T, U, UPRIME, QSEN, INFO, 
     *        RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)
C
c         write(20, *) qsen(1), qsen(2), qsen(3)
         IF (IDID .EQ. 1) GOTO 55
         UMAX = 0.0D0
         DO I = 1,NEQ/(info(inso+1)+1)
            UMAX = MAX (UMAX, ABS(U(I)) )
         END DO
         HU = RWORK(7)
         NQU = IWORK(8)
         NST = IWORK(11)
         NNI = IWORK(19)
         NLI = IWORK(20)
         t2 = secnds(t1)
         WRITE (LOUT,60) T,UMAX,NQU,HU,NST,NNI,NLI, t2
 60      FORMAT(E15.5,E12.4,I5,E14.3,I7,I9,I8, f7.2)
C
         IF (IDID .LT. 0) THEN
           WRITE(LOUT,65)T
 65        FORMAT(//' Final time reached =',E12.4//)
           GO TO 80
        ENDIF
C
        TOUT = TOUT*2.0D0
 70   CONTINUE
C
C Here we display some final statistics for the problem.
C The ratio of NLI to NNI is the average dimension of the Krylov
C subspace involved in the Krylov linear iterative method.
 80   CONTINUE
      NST = IWORK(11)
      NPE = IWORK(13)
c      NRE = IWORK(12) + NPE*MBAND
      NRE = IWORK(12)
      LIW = IWORK(17)
      LRW = IWORK(18)
      NNI = IWORK(19)
      NLI = IWORK(20)
      NPS = IWORK(21)
      NSE = IWORK(22)
      IF (NNI .NE. 0) AVDIM = REAL(NLI)/REAL(NNI)
      NETF = IWORK(14)
      NCFN = IWORK(15)
      NCFL = IWORK(16)
      WRITE (LOUT,90) LRW,LIW,NST,NRE+IPAR(30),NSE,NPE,NPS,NNI,NLI,
     1                AVDIM,NETF,NCFN,NCFL
 90   FORMAT(//' Final statistics for this run..'/
     1   ' RWORK size =',I8,'   IWORK size   =',I6/
     2   ' Number of time steps              =',I5/
     3   ' Number of residual evaluations    =',I5/
     3   ' Number of sensitivity evaluations =',I5/
     4   ' Number of Jac. or prec. evals.    =',I5/
     5   ' Number of preconditioner solves   =',I5/
     6   ' Number of nonlinear iterations    =',I5/
     7   ' Number of linear iterations       =',I5/
     8   ' Average Krylov subspace dimension =',F8.4/
     1   ' Error test failure .............. =',I5/
     1   I5,' nonlinear conv. failures,',I5,' linear conv. failures')
      WRITE(LOUT,100) LWPMIN, LIWPMIN
 100  FORMAT(' Minimum lengths for work arrays WP and IWP: ',i7,1x,i7)
C
C------  End of main program for DHEATILU example program --------------
      t2 = dtime(tt)
c      t3 = tt(1)
c      t2 = secnds(t1)
c      t3 = t2
      print *, (qsen(i), i=1,3)
      do i = 4, 18
         write(11,*) qsen(i)
      end do
      write(6,*) '   Time it took:', tt(1)+tt(2)
      STOP
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
      NEQ = IPAR(33)
      M = IPAR(34)
      m2 = m + 2
      NP = IPAR(35)
      DX = RPAR(3)

      do i = 1, (np+1)*neq
         u(i) = 0.0d0
         uprime(i) = 0.0d0
      end do

      ii = 20*M2 + 12
      do i = 4, np+1
         ii = ii + 1
         u((i-1)*neq + ii) = 1.0d0
      end do

C
      DO 20 K = 0,M+1
        YK = K*DX
        IOFF = (M+2)*K
        DO 10 J = 0,M+1
           XJ = J*DX
           I = IOFF + J + 1
           U(I) = 16.0D0*XJ*(1.0D0-XJ)*YK*(1.0D0-YK)
 10     CONTINUE
 20   CONTINUE

      RETURN
C------------  End of Subroutine UINIT  --------------------------------
      END

      SUBROUTINE RESH (T, U, UPRIME, CJ, DELTA, IRES, RPAR,IPAR,SENPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)

      if (ires .eq. 0) then
         call RESH1 (T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR,SENPAR)
      else if (ires. eq. 1) then
         call RESSA (T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR,SENPAR)
      else 
         ires = -2
      end if
      return
      end

      SUBROUTINE RESH1 (T, U, UPRIME, CJ, DELTA,IRES,RPAR,IPAR,SENPAR)
C
C This is the user-supplied RES subroutine for this example.
C It computes the residuals for the 2-D discretized heat equation,
C with zero boundary values.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)
C
C Set problem constants using IPAR and RPAR.
      NEQ = IPAR(33)
      M = IPAR(34)
      COEFF = RPAR(4)
      M2 = M + 2
C
C Load U into DELTA, in order to set boundary values.
      DO 10 I = 1,NEQ
 10     DELTA(I) = UPRIME(I)
C
C Loop over interior points, and load residual values.
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
           I = IOFF + J + 1
           TEMX = SENPAR(1)*(U(I-1)  + U(I+1)  - 2.0d0*u(i))
           TEMY = SENPAR(2)*(U(I-M2) + U(I+M2) - 2.0d0*u(i))
           DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF
 20     CONTINUE
 30   CONTINUE
C
      RETURN
C------------  End of Subroutine RESH  ---------------------------------
      END

      subroutine ressa(t,u,uprime,cj,delta,ires,rpar,ipar,SENPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)
C
C Set problem constants using IPAR and RPAR.
      NEQ = IPAR(33)
      M = IPAR(34)
      NP = ipar(35)
      COEFF = RPAR(4)
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
      subroutine jackdaspk
      return
      end
      
      SUBROUTINE QRES (T,U,UPRIME,QSEN,ires,RPAR,IPAR,senpar)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), RPAR(*), IPAR(*),senpar(*)
C
C Set problem constants using IPAR and RPAR.
      M = IPAR(34)
      M2 = M + 2
C
C Loop over interior points, and load residual values.
      qsen = 0.0
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
           I = IOFF + J + 1
           qsen = qsen + u(i)*u(i)
 20     CONTINUE
 30   CONTINUE
C
      RETURN
      END
