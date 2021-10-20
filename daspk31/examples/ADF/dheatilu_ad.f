C Copyright 2000 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  DHEATILU
C***REFER TO  DDASPK
C***DATE WRITTEN   951005   (YYMMDD)
C***REVISION DATE  951005   Reworking of heat example using banded J
C***REVISION DATE  990117   Sensitivity analysis for parameters
C***REVISION DATE  000714   revise for sensitivity analysis by adjoint method
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
C p_1 and p_2 are the sensitivity parameters. In this example, we also 
C take the initial conditions as sensitivity parameters. We have two
C derived functions in this example
C    g(1) = \int sum(u_i)
C    g(2) = sum_1^N (u_i*u_i)
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   UINIT, DDASPKADJOINT
C
C***END PROLOGUE  DHEATILU_AD
C
C Necessary declarations.  The dimension statements use a
C maximum value for the mesh parameter M.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RESH, DJACILU, DPSOLILU, QRES, AD_DJACILU,
     *         RES_ADY, RES_ADP, DPSOL
C Load parameters for sparse preconditioner.
      PARAMETER (LENPFAC = 5)   ! the average number of non-zeros
      PARAMETER (LENPLUFAC = 5) 
      PARAMETER (IPREMETH = 1)  ! =1 means ILUT preconditioner used
                                ! =2 means ILUTP preconditioner used
                                ! =3 ILUD
                                ! =4 ILUDP
      PARAMETER (LFILILUT = 5)
      PARAMETER (IREORDER = 1)
      PARAMETER (ISRNORM = 1)
      PARAMETER (NORMTYPE = 2)
      PARAMETER (JACOUT = 0)
      PARAMETER (JSCALCOL = 1)
      PARAMETER (TOLILUT = 0.001)
      PARAMETER (PERMTOL = 0.01)
C Load workspace lengths.
      PARAMETER (MAXM = 50, MAXM2 = MAXM+2, MXNEQ=MAXM2*MAXM2)
      PARAMETER (LENWP = 2*LENPFAC*MXNEQ + LENPLUFAC*MXNEQ
     .                   + ISRNORM*MXNEQ + 2*(MXNEQ+1) )
      PARAMETER (LENIWP = 4*(MXNEQ+1) + 3*LENPFAC*MXNEQ
     .                    + 2*LENPLUFAC*MXNEQ + IREORDER*2*MXNEQ
     .                    + (IPREMETH-1)*2*MXNEQ )
      PARAMETER (LRWF = 91 + 18*MXNEQ + LENWP)
      PARAMETER (LIWF = 40 + MXNEQ + LENIWP)
      PARAMETER (NBUF = 26*MXNEQ, MXNQ = 2, MXNCHK=4)
      PARAMETER (LENRW = 2*LRWF + 10*MXNEQ + 3*MXNEQ*MXNQ + NBUF)
      PARAMETER (LENIW = 2*LIWF + 4*MXNCHK)
      DIMENSION U(MXNEQ),UPRIME(MXNEQ)
      DIMENSION INFO(30), RPAR(4), SENPAR(2),  IPAR(50+3*MXNEQ)
      PARAMETER (LQSEN = MXNQ*(MXNEQ+1+MXNEQ))
      DIMENSION QSEN(LQSEN), IEOPT(2*MXNEQ), INFOB(10)
      DIMENSION RWORK(LENRW),IWORK(LENIW)
      real*4   t1, t2, t3, tt(2),secnds
C
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
      LRW = LENRW
      LIW = LENIW
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
      info(15) = 0
      if (info(15) .eq. 1) then
         IPAR(12) = 1           ! Jacobian evaluation method
         print *, ' Jacobian evaluation for ILU preconditioner:(0,1)'
         print *, '    0 --- finite differencing '
         print *, '    1 --- ADIFOR '
         read(*,*) i
         if (i .eq. 0) ipar(12) = 0
         IPAR(13) = 0           ! data dependence for rpar(*)
         IPAR(30) = 0
         IPAR(14) = 0           ! No saving for Krylov space
         IPAR(15) = 50          ! start of Adifor work space in PSOL
         IPAR(16) = 50 + 3*mxneq
         RPAR(1) = TOLILUT
         RPAR(2) = PERMTOL
C     Check IPAR, RPAR, LENWP and LENIWP for illegal entries and long
C     enough work array lengths.
         CALL DSPSETUP (NEQ, LENWP, LENIWP, RPAR, IPAR, IERR,
     .        LWPMIN, LIWPMIN)
         IF (IERR .NE. 0) THEN
            WRITE(LOUT,15) IERR
 15         FORMAT(' Error return from DSPSETUP: IERR = ',i5)
            IF (LWPMIN .GT. LENWP) THEN
               WRITE(LOUT,*) ' More WP work array length needed'
            ENDIF
            IF (LIWPMIN .GT. LENIWP) THEN
               WRITE(LOUT,*) ' More IWP work array length needed'
            ENDIF
            STOP
         ENDIF
      else
         IWORK(27) = 0
         IWORK(28) = 0
      end if          
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
C-----------------------------------------------------------------------
C
      DO 20 I = 1,30
 20     INFO(I) = 0
C
      INFO(6)  = 1
      INFO(12) = 1
C
C Call subroutine UINIT to initialize U and UPRIME.
C
      CALL UINIT (U, UPRIME, RPAR, IPAR)
C
C Here we set tolerances for DDASPK to indicate how much accuracy 
C we want in the solution, in the sense of local error control.
C For this example, we ask for pure absolute error control with a
C tolerance of 1.0D-3.
      RTOL = 1.0d-6
      ATOL = 1.0D-6
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
C
C-----------------------------------------------------------------------
C Now we solve the problem.
C
C We pass to DDASPK the names DJACILU and DPSOLILU for the JAC and PSOL
C routines to do the preconditioning.
C
C If DDASPK failed in any way (IDID .lt. 0) we print a message and
C stop the integration.
C-----------------------------------------------------------------------
C
      NOUT = 1
      T = 0.0D0
      TOUT = 0.16D0
      t1 = dtime(tt)
      t1 = secnds(0.0e0)
c---------------------------------------------------------------
c  sensitivity analysis
c---------------------------------------------------------------
c
c...  total number of sensitivity parameters
      info(19) = 116            
c
c...  number of parameters that appear in the RES routine
      info(22) = 2
c
c...  number of the quadrature variables in the forward integration
      info(28) = 1
      neq = neq + info(28)
c
c...  IEOPT 
      do i = 1, neq
         ieopt (i) = 1
         ieopt (i+neq) = 1
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
c...  required initial consistent values of ady and adyprime, set infob(2)=1:
      infob(2) = 0
      infob(7) = 120            ! estimate of the number of time steps
c
c...  total number of derived functions
      nq    = 2
c
c     run daspkadjoint program
c
      if (info(15) .eq. 1) then
         CALL DDASPKADJOINT(
     1        RESH, NEQ, T, U, UPRIME, TOUT, INFO, RTOL, ATOL,
     1        IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, 
     1        DJACILU, DPSOLILU, SENPAR, 
     1        ADRESH, NQ, QRES, QSEN, INFOB, RTOLB, ATOLB, NBUF,
     1        AD_DJACILU, IEOPT, RES_ADP, RES_ADY)
      else
         CALL DDASPKADJOINT(
     1        RESH, NEQ, T, U, UPRIME, TOUT, INFO, RTOL, ATOL,
     1        IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, 
     1        DJACILU, DPSOL, SENPAR, 
     1        ADRESH, NQ, QRES, QSEN, INFOB, RTOLB, ATOLB, NBUF,
     1        AD_DJACILU, IEOPT, RES_ADP, RES_ADY)
      end if
C
C Here we display some final statistics for the problem.
 80   CONTINUE
      t2 = dtime(tt)
      
      write(*,*)
      write(*,*) 'Performance statistics:'
      write(*,*) ' Time it took:', tt(1)+tt(2)
c
c     Here we display the values for the derived functions and their
c     sensitivity with respect to the sensitivity parameters p
c
      print *, ' value of integral sum(u_i) and sum(u_i^2) at t=', tout
      print *, '      ', (qsen(i), i = 1, nq)
      lsen = nq*(neq+1)+1
      print *, ' sensitivities are:'
      print *, '      ', ((qsen(lsen+i*nq+k), i=0,1),k=0,nq-1)
c
c     From lsen+2*nq, qsen(*) stores the sensitivity with
c     respect to the initial conditions.
C
C------  End of main program for DHEATILU_AD example program --------------
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
      DX = RPAR(3)

      do i = 1, neq+1
         u(i) = 0.0d0
         uprime(i) = 0.0d0
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
C
C This is the user-supplied RES subroutine for this example.
C It computes the residuals for the 2-D discretized heat equation,
C with zero boundary values. It also defines a quadrature equation
C for one of the derived functions. U(neq+1) is the quadrature 
C variables
C
      NEQ = IPAR(33)
      if (ires .eq. 0) then
         call RESH1 (T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR,SENPAR)
      else if (ires .eq. 3) then
         sum = 0.0d0
         do i = 1, neq
            sum = sum + u(i)
         enddo
         delta(neq+1) = uprime(neq+1) - sum
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
 10     DELTA(I) = uprime(I)
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

      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
C
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
C
      implicit none
      double precision t, y(*), yp(*), rpar(*), senpar(*),
     *     qsen(*)
      integer ipar(*)
      integer ires, i, neq

      neq = ipar(33)
      qsen(1) = y(neq+1)
      qsen(2) = 0.0d0
      do i = 1, neq
         qsen(2) = qsen(2) + y(i)*y(i)
      end do
      return
      end
C
C---- interface routine to the TAMC-generated code
C
      subroutine res_ady(t, y, ady, yp, cj, delta, addelta, 
     &     ires, rpar, ipar, senpar)
C
C This is the user-supplied RES_ADY subroutine for this example.
C It computes a vector-matrix product and return it in ady, i.e.
C     ady = addelta * F_y
C
      implicit none
      double precision t, y(*),  yp(*), cj, rpar(*), senpar(*),
     *     ady(*), addelta(*), delta(*)
      integer ires, ipar(*)
C
C     call TAMC generated routine
      call adyres(ires,  rpar, ipar, senpar, ady, addelta )
      
      return
      end
         
      subroutine res_adp(
     *     T, Y, YP, CJ, delta,ADDELTA, IRES, RPAR, IPAR, SENPAR,
     *     ADSENPAR)
C
C This is the user-supplied RES_ADP subroutine for this example.
C It computes a vector-matrix product and return it in ADSENPAR, i.e.
C     adsenpar = addelta * F_p
C
      implicit none
      double precision t, y(*), yp(*), cj, rpar(*), senpar(*),
     *     adsenpar(*), addelta(*), delta(*)
      integer ires, ipar(*)
C
C     call TAMC generated routine
      call adpres(y, ires, rpar, ipar, addelta, adsenpar)

      return
      end
