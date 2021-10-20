C Copyright 2000 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  DWEB_AD
C***REFER TO  DDASPK, DASPKADJOINT
C***DATE WRITTEN   950914   (YYMMDD)
C***MODIFIED TO INCLUDE SENSITIVITY ANALYSIS 990115
C***MODIFIED TO INCLUDE SENSITIVITY ANALYSIS by ADJOINT METHOD 000715
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
C***DESCRIPTION
C
C-----------------------------------------------------------------------
C Example program for DDASPKADJOINT.
C Index-1 DAE system derived from ns-species interaction PDE in 2 dimensions.
C Sensitivity analysis with respect to parameters alpha and beta.
C
C This is the double precision version.
C-----------------------------------------------------------------------
C
C This program solves a DAE system that arises from a system
C of partial differential equations.  The PDE system is a food web
C population model, with predator-prey interaction and diffusion on
C the unit square in two dimensions.  The dependent variable vector is
C
C         1   2        ns
C   c = (c , c , ..., c  )
C
C and the PDEs are as follows..
C
C     i               i      i
C   dc /dt  =  d(i)*(c    + c   )  +  R (x,y,c)  (i=1,...,ns/2)
C                     xx     yy        i
C
C                     i      i
C   0       =  d(i)*(c    + c   )  +  R (x,y,c)  (i=(ns/2)+1,...,ns)
C                     xx     yy        i
C
C where
C                  i          ns         j
C   R (x,y,c)  =  c *(b(i) + sum a(i,j)*c )
C    i                       j=1
C
C The number of species is ns = 2*np, with the first np being prey and
C the last np being predators.  The coefficients a(i,j), b(i), d(i) are
C
C   a(i,i) = -a  (all i)
C   a(i,j) = -g  (i .le. np, j .gt. np)
C   a(i,j) =  e  (i .gt. np, j .le. np)
C   b(i) =  b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .le. np)
C   b(i) = -b*(1 + alpha*x*y + beta*sin(4*pi*x)*sin(4*pi*y))  (i .gt. np)
C   d(i) = dprey  (i .le. np)
C   d(i) = dpred  (i .gt. np)
C
C The various scalar parameters are set in subroutine setpar.
C
C The boundary conditions are.. normal derivative = 0.
C A polynomial in x and y is used to set the initial conditions.
C
C The PDEs are discretized by central differencing on a MX by MY mesh.
C
C The DAE system is solved by DDASPK direct band method for the linear 
C systems (internal Jacobian).
C
C-------------------------------sensitivity analysis--------------------
C For sensitivity analysis, alpha and beta are the two sensitivity 
C parameters. We consider two same derived functions in this example,
C    g(1) = g(2) = \sum (u_i*u_i) 
C DDASPKADJOINT can also output the sensitivities with the initial 
C conditions with little effort.
C
C-----------------------------------------------------------------------
C Note.. in addition to the main program and subroutines given below,
C this program requires the BLAS routine DAXPY.
C-----------------------------------------------------------------------
C References
C [1] Peter N. Brown and Alan C. Hindmarsh,
C     Reduced Storage Matrix Methods in Stiff ODE Systems,
C     J. Appl. Math. & Comp., 31 (1989), pp. 40-91.
C [2] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
C     Using Krylov Methods in the Solution of Large-Scale Differential-
C     Algebraic Systems, SIAM J. Sci. Comput., 15 (1994), pp. 1467-1488.
C [3] Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
C     Consistent Initial Condition Calculation for Differential-
C     Algebraic Systems, LLNL Report UCRL-JC-122175, August 1995;
C     submitted to SIAM J. Sci. Comp.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C***END PROLOGUE  DWEB_AD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RESWEB, ADRESWEB, QRES, RES_ADP, I_RESW, 
     *         RES_ADY, RES_ADYP
C
C Dimension solution arrays and work arrays.
C
C  with INFO(5) = 0, INFO(6) = 1:
C    The length required for RWORK is
C              50 + (2*ML+MU+11)*NEQ + 2*(NEQ/(ML+MU+1) + 1) + 4*NY
C    For MX = MY = (even number) and ML = MU = NS*MX, this length is
C              50 + (3*NS*MX + 11)*NEQ + MY + 4*NY
C    The length required for IWORK is  40 + 2*NEQ .
C
C
C The dimensions for the various arrays are set below using parameters
C   MAXN    which must be .ge. NEQ = (NPARM+1)*NS*MX*MY,
C   MAXS    which must be .ge. NS,
C   MAXM    which must be .ge. MAX(MX,MY).
C
C Note that for sensitivity analysis, NY=NS*MX*MY, MAXP=INFO(19)
C is the number of problem parameters, NEQ is the total number of
C equations, including sensitivity equations.
C
      PARAMETER (MAXP=2,MAXS=2,MAXM=20, MAXN=MAXS*MAXM*MAXM,
     1     LRWF = 50 + (3*MAXS*MAXM + 11)*MAXN + MAXM + 4*MAXN,
     2     LIWF = 40 + 2*MAXN + 4*MAXN)
C
      DIMENSION CC(MAXN), CCPRIME(MAXN), 
     1     INFO(30), RPAR(MAXN), IPAR(4), SENPAR(2)
C
C The COMMON blocks /PPAR1/ and /PPAR2/ contain problem parameters.
C
      COMMON /PPAR1/ AA, EE, GG, BB, DPREY, DPRED
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
      PARAMETER (MXNCHK=4, NBUF = MXNCHK*12*MAXN, MXNQ=2) 
      PARAMETER (LRW = 2*LRWF + 10*MAXN + 3*MAXN*MXNQ + NBUF +
     1     MXNCHK*(40+MAXN) + 1070)
      PARAMETER (LIW = 2*LIWF + 4*MXNCHK + MAXN)
      DIMENSION QSEN((MAXP+1+MAXN)*MXNQ), IEOPT(2*MAXN), INFOB(10)
      DIMENSION RWORK(LRW), IWORK(LIW)
      real t1, t2, t3, tt(2)
C
C Call SETPAR to set basic problem parameters.
      CALL SETPAR(SENPAR)
C
      ALPH = SENPAR(1)
      BETA = SENPAR(2)
C Set remaining problem parameters.
      NPARM = 2
      NEQ = NS*MX*MY
      MXNS = MX*NS
      DX = AX/REAL(MX-1)
      DY = AY/REAL(MY-1)
      DO 10 I = 1,NS
        COX(I) = DIFF(I)/DX**2
 10     COY(I) = DIFF(I)/DY**2
C
      WRITE(*,20)NS
 20   FORMAT(//' Example program for DDASPK package'//
     1   ' Food web problem with NS species, NS =',I4/
     2   ' Predator-prey interaction and diffusion on a 2-D square'/)
      WRITE(*,30) AA,EE,GG,BB,DPREY,DPRED, ALPH,BETA
 30   FORMAT(' Matrix parameters..  a =',E12.4,'   e =',E12.4,
     1   '   g =',E12.4/21x,' b parameter =',E12.4//
     2   ' Diffusion coefficients.. dprey =',E12.4,'   dpred =',E12.4//
     3   ' Rate parameters alpha =',E12.4,' and beta =',E12.4/)
      WRITE(*,40) MX,MY,NEQ
 40   FORMAT(' Mesh dimensions (MX,MY) =',2I4,
     1   5x,' Total system size is NEQ =',I7/)
C
C Here set the flat initial guess for the predators.
c=================================================================
      PREDIC = 1.0D2
c=================================================================
C
C Set remaining method parameters for DDASPK.
C These include the INFO array and tolerances.
C 
      DO 50 I = 1,30
 50   INFO(I) = 0
C
C Here set INFO(11) = 1, indicating I.C. calculation requested.
      INFO(11) = 1
C
C Here set INFO(14) = 1 to get the computed initial values.
      INFO(14) = 1
C
C Here set INFO(16) = 1 to get alternative error test (on the
C differential variables only).
      INFO(16) = 1
C
C Here set INFO(19) = 2 to indicate that a sensitivity analysis
C is to be done with two parameters
      inso = 18
      INFO(inso+1) = 2          ! two sensitivity parameters
      INFO(inso+2) = 0          ! 0 finite difference
      INFO(inso+4) = 2          ! two parameters in the RES routine
      NY = NEQ
C
C Here set the tolerances.      
      RTOL = 1.0D-5
      ATOL = RTOL
C
      WRITE(*,70)RTOL,ATOL,INFO(11),PREDIC,INFO(16)
 70   FORMAT(' Tolerance parameters.. RTOL =',E10.2,'   ATOL =',E10.2//
     1   ' Internal I.C. calculation flag INFO(11) =',I2,
     2   '   (0 = off, 1 = on)'/
     3   ' Predator I.C. guess =',E10.2//
     4   ' Alternate error test flag INFO(16) =',I2,
     5        '  (0 = off, 1 = on)')
C
      CALL SETID (MX, MY, NS, NP, 40, IWORK)
C
C
C In the case of the direct method, set INFO(6) = 1 to signal a banded 
C Jacobian, set IWORK(1) = IWORK(2) = MX*NS, the half-bandwidth, and
C call SETID to set the IWORK segment ID indicating the differential
C and algebraic components.
      INFO(6) = 1
      info(5) = 0
      info(5)  = 0
      print *, ' Jacobian evaluation method(0,2):'
      print *, '   0 --- finite differencing'
      print *, '   2 --- ADIFOR with SparsLinC'
      read(*,*) info(5)
      IWORK(1) = MXNS
      IWORK(2) = MXNS
      WRITE(*,90)MXNS
 90   FORMAT(' Banded Jacobian,', ' half-bandwidths =',I4)
C
C Set the initial T and TOUT, and call CINIT to set initial values.
      T = 0.0D0
      CALL CINIT (CC, CCPRIME, PREDIC, RPAR, SENPAR)
      lcout =10
      NLI = 0
      NNI = 0
      t2 = dtime(tt)
      tout = 5.0d0
c---------------------------------------------------------------
c  sensitivity analysis by adjoint method
c---------------------------------------------------------------
c...  IEOPT 
      DO JY = 1,MY
         IYOFF = MXNS*(JY-1)
         DO JX = 1,MX
            IC0 = IYOFF + NS*(JX-1)
            DO I = 1,NS
               ICI = IC0 + I
               IF (I .GT. NP) THEN
                  ieopt(ici) = -1
               ELSE
                  ieopt(ici) = 1
               ENDIF
            end do
         end do
      end do
      do i = 1, neq
         ieopt(neq + i) = iwork(40+i)
      end do
C
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
c...  total number of derived functions
      nq    = 2
c
c...  initialization for an index-1 system
      infob(2) = 1
c
c...  evaluation methods for the index-1 adjoint equations
      infob(3) = 3              
      print *, ' Evaluation method for the adjoint equations:(1,3):'
      print *, '     1 --- input RES_ADY and RES_ADYP for index-1 DAE,'
      print *, '     3 --- input user-defined ADRES'
      read(*,*) infob(3)
C      
C=======================================================================
C
      if (infob(3) .eq. 3) then
         CALL DDASPKadjoint (
     1        RESWEB, NEQ, T, CC, CCPRIME, TOUT, INFO, RTOL, ATOL, 
     1        IDID, RWORK,LRW, IWORK,LIW, RPAR, IPAR, I_RESW, PSOLRS,
     1        SENPAR, 
     1        ADRESWEB, NQ, QRES, QSEN, INFOB, RTOLB, ATOLB, NBUF,
     1        I_RESW, IEOPT, RES_ADP, RES_ADY)
      else if (infob(3) .eq. 1) then
         CALL DDASPKadjoint (
     1        RESWEB, NEQ, T, CC, CCPRIME, TOUT, INFO, RTOL, ATOL, 
     1        IDID, RWORK,LRW, IWORK,LIW, RPAR, IPAR, I_RESW, PSOLRS,
     1        SENPAR, 
     1        RES_ADYP, NQ, QRES, QSEN, INFOB, RTOLB, ATOLB, NBUF,
     1        I_RESW, IEOPT, RES_ADP, RES_ADY)
      end if

c      CALL OUTWEB (T, CC, NS, MX, MY, LCOUT)
c      CALL OUTWEB (T, qsen(2), NS, MX, MY, LCOUT+1)

      t3 = dtime(tt)
      write(*,*)
      write(*,*) 'Performance statistics:'
      write(*,*) ' Time it took:', tt(1)+tt(2)
      print *, ' RWORK size = ', iwork(40)+iwork(iwork(39)+18)
      print *, ' IWORK size = ', iwork(39)+iwork(iwork(39)+17)
 300  continue
c
c     Here we display the values for the derived functions and their
c     sensitivity with respect to the sensitivity parameters p
c
      print *, ' The values for the two derived functions (q1 = q2):'
      print *, '        ', (qsen(i), i = 1, nq)
      lsen = nq*(neq+1)+1
      print *, ' The sensitivities of the derived functions with two',
     *     ' parameters:'
      print *, '        ', ((qsen(lsen+i*nq+k), i=0,1),k=0,nq-1)

      STOP
C------  End of main program for DWEB_AD example program ------------------
      END

      SUBROUTINE SETPAR(SENPAR)
C-----------------------------------------------------------------------
C This routine sets the basic problem parameters, namely
C AX, AY, NS, MX, MY,  problem coefficients ACOEF, BCOEF, DIFF,
C ALPH, BETA, using parameters NP, AA, EE, GG, BB, DPREY, DPRED.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXP=2,MAXS=2)
      DIMENSION SENPAR(*)
      COMMON /PPAR1/ AA, EE, GG, BB, DPREY, DPRED
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      AX = 1.0D0
      AY = 1.0D0
      NP = 1
      MX = 20
      MY = 20
      AA = 1.0D0
      EE = 1.0D4
      GG = 0.5D-6
      BB = 1.0D0
      DPREY = 1.0D0
      DPRED = 0.05D0
      ALPH = 50.0D0
      BETA = 100.0D0
      SENPAR(1) = ALPH
      SENPAR(2) = BETA
      NS = 2*NP
      DO 20 J = 1,NP
        DO 10 I = 1,NP
          ACOEF(NP+I,J) = EE
          ACOEF(I,NP+J) = -GG
 10       CONTINUE
        ACOEF(J,J) = -AA
        ACOEF(NP+J,NP+J) = -AA
        BCOEF(J) = BB
        BCOEF(NP+J) = -BB
        DIFF(J) = DPREY
        DIFF(NP+J) = DPRED
 20     CONTINUE
      PI = 3.141592653589793D0
      FPI = 4.0D0*PI
C
      RETURN
C------------  End of Subroutine SETPAR  -------------------------------
      END

      SUBROUTINE SETID (MX, MY, NS, NSD, LID, IWORK)
C-----------------------------------------------------------------------
C This routine sets the ID array in IWORK, indicating which components
C are differential and which are algebraic.
C-----------------------------------------------------------------------
      DIMENSION IWORK(*)
C
      NSDP1 = NSD + 1
      DO 40 JY = 1,MY
        I00 = MX*NS*(JY-1) + LID
        DO 30 JX = 1,MX
          I0 = I00 + NS*(JX-1)  
          DO 10 I = 1,NSD
 10         IWORK(I0+I) = 1
          DO 20 I = NSDP1,NS
 20         IWORK(I0+I) = -1
 30       CONTINUE
 40     CONTINUE

      RETURN
C------------  End of Subroutine SETID --------------------------------
      END

      SUBROUTINE CINIT (CC, CCPRIME, PREDIC, rpar, SENPAR)
C-----------------------------------------------------------------------
C This routine computes and loads the vectors of initial values.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(*), CCPRIME(*), RPAR(*), SENPAR(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C

      ALPH = SENPAR(1)
      BETA = SENPAR(2) 
C Load CC.
      NPP1 = NP + 1
      DO 30 JY = 1,MY
        Y = REAL(JY-1)*DY
        ARGY = 16.0D0*Y*Y*(AY-Y)*(AY-Y)
        IYOFF = MXNS*(JY-1)
        DO 20 JX = 1,MX
          X = REAL(JX-1)*DX
          ARGX = 16.0D0*X*X*(AX-X)*(AX-X)
          IOFF = IYOFF + NS*(JX-1)
          FAC = 1.0D0 + ALPH*X*Y + BETA*SIN(FPI*X)*SIN(FPI*Y)
          DO 10 I = 1,NP
 10         CC(IOFF + I) = 10.0D0 + REAL(I)*ARGX*ARGY
          DO 15 I = NPP1,NS
 15         CC(IOFF + I) = PREDIC
 20       CONTINUE
 30     CONTINUE
C
C Load CCPRIME.
      T = 0.0D0
      CALL FWEB (T, CC, CCPRIME, RPAR, SENPAR)
      DO 60 JY = 1,MY
        IYOFF = MXNS*(JY-1)
        DO 50 JX = 1,MX
          IOFF = IYOFF + NS*(JX-1)
          DO 40 I = NPP1,NS
 40         CCPRIME(IOFF+I) = 0.0D0
 50     CONTINUE
 60   CONTINUE
C
      RETURN
C------------  End of Subroutine CINIT  --------------------------------
      END

      SUBROUTINE OUTWEB (T, C, NS, MX, MY, LUN)
C-----------------------------------------------------------------------
C This routine prints the values of the individual species densities
C at the current time T, to logical unit LUN.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION C(NS,MX,MY)
      DO 30 JY = MY,1,-1
         do jx = 1, mx
            write(lun, *) jx, jy, (c(i,jx, jy), i=1,ns)
         end do
         write(lun,*)
 30   CONTINUE
      write(lun, *)
C     
      RETURN
C------------  End of Subroutine OUTWEB  -------------------------------
      END
      SUBROUTINE RESWEB (
     *     T, U, UPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C-----------------------------------------------------------------------
C This routine computes the residual vector, using Subroutine FWEB
C for the right-hand sides.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),UPRIME(*),DELTA(*),RPAR(*),IPAR(*),SENPAR(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      CALL FWEB (T, U, DELTA, RPAR,SENPAR)
C
      DO 30 JY = 1,MY
         IYOFF = MXNS*(JY-1)
         DO 20 JX = 1,MX
            IC0 = IYOFF + NS*(JX-1)
            DO 10 I = 1,NS
               ICI = IC0 + I
               IF (I .LE. NP) THEN
                  DELTA(ICI) = UPRIME(ICI) - DELTA(ICI)
               ELSE
                  DELTA(ICI) = -DELTA(ICI)
               ENDIF
 10         CONTINUE
 20      CONTINUE
 30   CONTINUE
C
      RETURN
C------------  End of Subroutine RESWEB  -------------------------------
      END
      SUBROUTINE FWEB (T, CC, CRATE, RPAR, SENPAR)
C-----------------------------------------------------------------------
C This routine computes the right-hand sides of all the equations
C and returns them in the array CRATE.
C The interaction rates are computed by calls to WEBR, and these are
C saved in RPAR(1),...,RPAR(NEQ) for use later in preconditioning.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CC(*), CRATE(*), RPAR(*), SENPAR(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      DO 60 JY = 1,MY
        IYOFF = MXNS*(JY-1)
        IDYU = MXNS
        IF (JY .EQ. MY) IDYU = -MXNS
        IDYL = MXNS
        IF (JY .EQ. 1) IDYL = -MXNS
        DO 40 JX = 1,MX
          IC = IYOFF + NS*(JX-1) + 1
C Get interaction rates at one point (X,Y).
          CALL WEBR (T, JX, JY, CC(IC), RPAR(IC), SENPAR)
          IDXU = NS
          IF (JX .EQ. MX) IDXU = -NS
          IDXL = NS
          IF (JX .EQ. 1) IDXL = -NS
          DO 20 I = 1,NS
            ICI = IC + I - 1
C Do differencing in Y.
            DCYLI = CC(ICI) - CC(ICI-IDYL)
            DCYUI = CC(ICI+IDYU) - CC(ICI)
C Do differencing in X.
            DCXLI = CC(ICI) - CC(ICI-IDXL)
            DCXUI = CC(ICI+IDXU) - CC(ICI)
C Collect terms and load CRATE elements.
            CRATE(ICI) = COY(I)*(DCYUI - DCYLI) + COX(I)*(DCXUI - DCXLI)
     1                  + RPAR(ICI)
 20         CONTINUE
 40       CONTINUE
 60    CONTINUE
      RETURN
C------------  End of Subroutine FWEB  ---------------------------------
      END

      SUBROUTINE WEBR (T, JX, JY, C, CRATE, SENPAR)
C-----------------------------------------------------------------------
C This routine computes one block of the interaction term R of the 
C system, namely block (JX,JY), for use in preconditioning.
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SENPAR(*)
      DIMENSION C(*), CRATE(*)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      ALPH = SENPAR(1)
      BETA = SENPAR(2)
c      ALPH = 50.0D0
c      BETA = 100.0D0
      Y = REAL(JY-1)*DY
      X = REAL(JX-1)*DX
      DO 10 I = 1,NS
 10     CRATE(I) = 0.0D0
      DO 15 J = 1,NS
        CALL DAXPY (NS, C(J), ACOEF(1,J), 1, CRATE, 1)
 15     CONTINUE
      FAC = 1.0D0 + ALPH*X*Y + BETA*SIN(FPI*X)*SIN(FPI*Y)
      DO 20 I = 1,NS
         CRATE(I) = C(I)*(BCOEF(I)*FAC + CRATE(I))
c         CRATE(I) = 100.d0*(BCOEF(I)*FAC + CRATE(I))
 20      continue
      RETURN
C------------  End of Subroutine WEBR  ---------------------------------
      END

C
C------------- Interface to the ADIFOR-generated routines --------------
C
      subroutine i_resw(t, y, i_y, yp, i_yp, cj, delta, i_delta, ires,
     *     rpar, ipar, senpar)
C-----------------------------------------------------------------------
C This routine is a wrapper to the ADIFOR-generated routine for the Jacobian
C evaluation
C-----------------------------------------------------------------------
      implicit none
      integer  ires, ipar(*), i_y(*), i_yp(*), i_delta(*)
      double precision t, y(*), yp(*), delta(*), cj, rpar(*), senpar(*)
c
c...  local variables 
      integer i_rpar(800)
      data i_rpar/800*0/
      save i_rpar
c
c...  This is the ADIFOR-generated routine with SparseLinC option
c
      call i_res(t, y, i_y, yp, i_yp, cj, delta, i_delta, ires,
     *     rpar, i_rpar, ipar, senpar)

      return
C------------  End of Subroutine I_RESW  ---------------------------------
      end
C
C------------- Interface to the TAMC-generated routines --------------
C
      SUBROUTINE adresweb (
     *           t, adY, adYP, CJ, delta, IRES, RPAR,IPAR,SENPAR, 
     *           adi)
C-----------------------------------------------------------------------
C This routine computes the residual vector for the adjoint DAE
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ady(*), adyp(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)
      dimension adi(*)
      dimension addelta(1000), adsenpar(2), tmp(1000)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
      
      NEQ = NS*MX*MY

      if (ires .eq. 0) then
         do i = 1, neq
            delta(i) = 0.0d0
            addelta(i) = ady(i)
         end do

         call res_ady(t, adi, delta, adi(ny+1), cj, tmp,addelta, ires,  
     *        rpar, ipar, senpar)

         DO 30 JY = 1,MY
            IYOFF = MXNS*(JY-1)
            DO 20 JX = 1,MX
               IC0 = IYOFF + NS*(JX-1)
               DO 10 I = 1,NS
                  IF (I .le. NP) THEN
                     ICI = IC0 + I
                     DELTA(ICI) = adyp(ici) + delta(ici)
                  ENDIF
 10            CONTINUE
 20         CONTINUE
 30      CONTINUE
      else if (ires .eq. 3) then
c
c..   compute \dot c = ady*F_p      
         do i = 1, neq
            addelta(i) = ady(i)
         enddo
         do i = 1, 2
            delta(neq+i) = 0.0d0
         end do
         
         call res_adp(t, adi, adi(ny+1), cj, tmp, addelta, ires, rpar,  
     *        ipar, senpar, delta(neq+1))

         do i = 1, 2
            delta(neq +i) = adyp(neq+i) - delta(neq+i)
         end do
      end if
C
      return
C------------  End of Subroutine ADRESWEB  ---------------------------------
      end
      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
C
C
C This is the user-supplied QRES subroutine for this example.
C It computes the values of the derived funtions. You can also
C define the derivatives of the derived functions with respect 
C to Y and SENPAR if you choose INFO(20) = 2.
C
      implicit double precision(a-h,o-z)
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
C
      qsen(1) = 0.0d0
      do i = 1, ns*my*mx
         qsen(1) = qsen(1) + y(i)*y(i)
c         if (ires .eq. 1) qsen(1+i) = 2.0d0*y(i)
      end do
      qsen(2) = qsen(1)
C
      return
      end
c
c...  interface to the TAMC generated routines
c
      subroutine res_ady(
     *     T, Y, ADY, YP, CJ, delta, ADDELTA,IRES, RPAR, IPAR, SENPAR)
C--------------------------------------------------------------
C This is the user-supplied RES_ADY subroutine for this example.
C It computes a vector-matrix product and return it in ady, i.e.
C     ady = addelta * F_y
C--------------------------------------------------------------
      implicit none
      double precision t, y(*), yp(*), cj, rpar(*), senpar(*),
     *     ady(*), addelta(*), adrpar(1000), delta(*)
      integer ires, ipar(*), i

      do i = 1, 800
         adrpar(i) = 0.0d0
      end do
C
C     call TAMC generated routine
      call adyres(y, rpar, senpar, ady, addelta, adrpar)

      return
      end

      subroutine res_adyp(
     *     T, Y, YP, ADYP, CJ, delta, ADDELTA,IRES, RPAR, IPAR, SENPAR)
C-------------------------------------------------------------------
C This is the user-supplied RES_ADP subroutine for this example.
C It computes a vector-matrix product and return it in ADYP, i.e.
C     adyp = addelta * F_y'
C--------------------------------------------------------------------
      implicit none
      double precision t, y(*), yp(*), cj, rpar(*), senpar(*),
     *     adyp(*), addelta(*), delta(*)
      integer ires, ipar(*)
C
C     call TAMC generated routine
      call adypres(adyp,addelta)

      return
      end

      subroutine res_adp(
     *     T, Y, YP, CJ, delta, ADDELTA, IRES, RPAR, IPAR, SENPAR,
     *     ADSENPAR)
C-------------------------------------------------------------------
C This is the user-supplied RES_ADP subroutine for this example.
C It computes a vector-matrix product and return it in ADSENPAR, i.e.
C     adsenpar = addelta * F_p
C--------------------------------------------------------------------
      implicit none
      double precision t, y(*), yp(*), cj, rpar(*), senpar(*),
     *     adsenpar(*), addelta(*), adprpar(1000), delta(*)
      integer ires, ipar(*), i

      do i = 1, 800
         adprpar(i) = 0.0d0
      end do
C
C     call TAMC generated routine
      call adpres( y, addelta, adprpar, adsenpar )

      return
      end



