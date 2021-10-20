C Work performed under the auspices of the U.S. Department of Energy
C by Lawrence Livermore National Laboratory under contract number 
C W-7405-Eng-48.
C
C Copyright 2000 the Regents of the University of California.
C All rights reserved.
C
C***BEGIN PROLOGUE  DWEBS
C***REFER TO  DDASPK
C***DATE WRITTEN   950914   (YYMMDD)
C***MODIFIED TO INCLUDE SENSITIVITY ANALYSIS 990115
C***MODIFIED TO INCLUDE SENSITIVITY ANALYSIS by ADJOINT METHOD 990115
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
C Example program for DDASPK.
C DAE system derived from ns-species interaction PDE in 2 dimensions.
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
C systems (internal Jacobian),
C
C-----------------------------------------------------------------------
C Note.. in addition to the main program and subroutines given below,
C this program requires the BLAS routine DAXPY. A copy of this routine
C resides in resweb.f as local_DAXPY so that the Tapenade engine works.
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
C***END PROLOGUE  DWEB
C
!       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE
      EXTERNAL RESWEB, QRES, ADRESWEB
      EXTERNAL RES_ADY, RES_ADP, RES_ADYP
      INTEGER  MX, MY, MXNS
      DOUBLE PRECISION AA, EE, GG, BB, DPREY, DPRED, STUB
      INTEGER NP, NS
      DOUBLE PRECISION AX, AY, ACOEF, BCOEF, DX, DY, FPI, DIFF
      DOUBLE PRECISION COX, COY
      DOUBLE PRECISION DT
      INTEGER          NUMSTEPS

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
C Note that NY=NS*MX*MY, NPARM=INFO(19) is the number of problem parameters, and NEQ=NY.
C
      INTEGER MAXP, MAXS, MAXM, MAXN, LRWF, LIWF, LENSP
      PARAMETER (MAXP=2,MAXS=2,MAXM=20, MAXN=MAXS*MAXM*MAXM,
     1     LRWF = 50 + (3*MAXS*MAXM + 11)*MAXN + MAXM + 4*MAXN,
     2     LIWF = 40 + 2*MAXN + 4*MAXN, LENSP=2)
C
      DOUBLE PRECISION CC, CCPRIME, RPAR, SENPAR
      INTEGER          INFO, IPAR
      DIMENSION CC((MAXP+1)*MAXN), CCPRIME(MAXN), 
     1     INFO(30), RPAR(MAXN), IPAR(4), SENPAR(LENSP)
C
C The COMMON blocks /PPAR1/ and /PPAR2/ contain problem parameters.
C
      COMMON /PPAR1/ AA, EE, GG, BB, DPREY, DPRED
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
      INTEGER MXNCHK, NBUF, MXNQ
      PARAMETER (MXNCHK=4, NBUF = MXNCHK*12*MAXN, MXNQ=2) 
      INTEGER LRW
      PARAMETER (LRW = 2*LRWF + 10*MAXN + 3*MAXN*MXNQ + NBUF +
     1     MXNCHK*(40+MAXN)+10)
      INTEGER LIW
      PARAMETER (LIW = 2*LIWF + 4*MXNCHK + MAXN)
      INTEGER          IEOPT, INFOB, IWORK
      DOUBLE PRECISION QSEN, RWORK, T, TOUT, TFINAL
      DIMENSION QSEN((MAXP+1+MAXN)*MXNQ), IEOPT(2*MAXN), INFOB(10)
      DIMENSION RWORK(LRW), IWORK(LIW)
      real t1, t2, t3, tt(2)
      DOUBLE PRECISION ALPH, BETA, PREDIC, RTOL, ATOL, RTOLB, ATOLB
      INTEGER          NPARM, NEQ, I, INSO, NY, LCOUT, NLI, NNI, JY
      INTEGER          IYOFF, JX, IC0, ICI, NQUAD, IDID, LSEN, K
      INTEGER          MXNEQAD, NEQAD
C     NEQAD is set such that
C                   IF INFOBI(4) .EQ. 0 .AND. INFO(22) .GT> 0 THEN
C                      NEQAD = NQ * (NEQ + INFO(22))
C                   ELSE
C                      NEQAD = NQ * NEQ
C                   END IF
C     where NQ is the number of derived functions, and NEQ is the number
C     of state variables.
      PARAMETER (MXNEQAD=MXNQ*(MAXN+LENSP))
      DOUBLE PRECISION ADY(MXNEQAD), ADYP(MXNEQAD)
      INTEGER TOTI, LRWORKB, LIWORKB, LICONT
!       DOUBLE PRECISION ADY, ADYP ! Since we aren't getting the ADY/ADYP values, 
C                                  we will just leave these variables as scalars
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
      PRINT *, 'MAXP=', MAXP, ' MAXN=', MAXN
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
C is to be done with two parameters.
      inso = 18
      INFO(inso+1) = 2          ! one sensitivity parameter
      INFO(inso+2) = 0          ! 0 finite difference
      INFO(inso+4) = 2          ! two parameters in the RES routine
      NY = NEQ
C
C Here we set IWORK(38) to indicate RPAR is used to pass the information
C of Y or YPRIME.
      IWORK(38) = MAXN

C
C Here set the tolerances.      
      RTOL = 1.0D-5
      ATOL = RTOL
      RTOLB = RTOL
      ATOLB = ATOL
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
      IF (INFO(12) .EQ. 0) THEN
        INFO(6)  = 1
        INFO(5)  = 0
        IWORK(1) = MXNS
        IWORK(2) = MXNS
        WRITE(*,90)MXNS
 90     FORMAT(' Banded Jacobian,', ' half-bandwidths =',I4)
      ENDIF
C
C Set the initial T and TOUT, and call CINIT to set initial values.
      T = 0.0D0
      CALL CINIT (CC, CCPRIME, PREDIC, RPAR, SENPAR)
      lcout =10
      NLI = 0
      NNI = 0
      t2 = dtime(tt)
      tout = 5.0d0
      tfinal = tout
      dt = 1.0d0
      numsteps = 5
c
c===================set-up for adjoint method===========================
c
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
      nquad    = 2              ! number of derived functions
      infob(1) = 0              ! default error tolerance
      infob(2) = 1              ! initialization for index-1
      infob(3) = 1              ! index-1 adjoint equation evaluation method,
      print *, ' Evaluation method for the adjoint equations:(1,3):'
      print *, '     1 --- input RES_ADY and RES_ADYP for index-1 DAE,'
      print *, '     3 --- input user-defined ADRES'
      read(*,*) infob(3)
      infob(4) = 1              ! take quadrature as variables
      infob(5) = 0              ! not index-2
      infob(6) = 0              ! not calculated inside DASPK
      infob(7) = 0             ! don't know how many steps we need
      
      IF (INFOB(4) .EQ. 0 .AND. INFO(22) .GT. 0) THEN
          NEQAD = NQUAD * (NEQ + INFO(22))
      ELSE
          NEQAD = NQUAD * NEQ
      END IF
C      
c=======================================================================
c
c Call DASPK_Adjoint
c
      if (infob(3) .eq. 3) then
             CALL DDASPKadjoint (
     1        RESWEB, NEQ, T, CC, CCPRIME, TOUT, TFINAL, INFO, 
     1        RTOL, ATOL, 
     1        IDID, RWORK,LRW, IWORK,LIW, RPAR, IPAR, STUB, 
     1        STUB, SENPAR, adresweb, nquad, Qres, Qsen, 
     1        infob, rtolb, atolb, nbuf,
     1        STUB, IEopt, RES_ADP, RES_ADY, 
     1        NEQAD, ADY, ADYP)
      else if (infob(3) .eq. 1) then
             CALL DDASPKadjoint (
     1        RESWEB, NEQ, T, CC, CCPRIME, TOUT, TFINAL, INFO, 
     1        RTOL, ATOL, 
     1        IDID, RWORK, LRW, IWORK,LIW, RPAR, IPAR, STUB,
     1        STUB, SENPAR, RES_ADYP, nquad, Qres, Qsen, 
     1        infob, rtolb, atolb, nbuf,
     1        STUB, IEopt, RES_ADP, RES_ADY, 
     1        NEQAD, ADY, ADYP)
      end if
        
      IF (idid .lt. 0) GOTO 300

      t3 = dtime(tt)
      write(*,*)
      write(*,*) 'Performance statistics:'
      write(*,*) ' Time it took:', tt(1)+tt(2)
      TOTI = IWORK(18) ! length of first chunk of RWORK (RWORKF)
      LICONT  = IWORK(39)
      LIWORKB = IWORK(LICONT+1)
      LRWORKB = IWORK(40)
      IF (LRWORKB.NE.1) THEN
        TOTI = TOTI + IWORK(LIWORKB+18) ! length of RWORKB
      END IF
      TOTI = TOTI + LRW - IWORK(LICONT+18) + 1 ! length of buffer
      print *, ' RWORK size = ', TOTI
      print *, 'length of RWORKF ', IWORK(18)
      print *, 'LRWORKB ', LRWORKB
      print *, 'length of RWORKB ', IWORK(LIWORKB+18)
      print *, ' Buffer size = ', (LRW - IWORK(LICONT+18) + 1), LRW
      TOTI = IWORK(17) ! length of first section of IWORK (IWORKF)
      TOTI = TOTI + 90 ! length of ICONT
      IF (LIWORKB .NE. 1) THEN
        TOTI = TOTI + IWORK(LIWORKB+17) ! len of IWORKB
      END IF
      print *, ' IWORK size = ', TOTI
 300  continue
      print *, ' The values for the two derived functions (q1 = q2):'
      print *, '        ', (qsen(i), i = 1, nquad)
      lsen = nquad*(neq+1)+1
      print *, ' The sensitivities of the derived functions with two',
     *     ' parameters:'
      print *, '        ', ((qsen(lsen+i*nquad+k), i=0,1),k=0,nquad-1)

      STOP
C------  End of main program for DWEB example program ------------------
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
      END
C------------  End of Subroutine OUTWEB  -------------------------------

C------------  Subroutine ADRESWEB  -------------------------------
      SUBROUTINE adresweb (
     *           t, adY, adYP, CJ, delta, IRES, RPAR,IPAR,SENPAR, 
     *           adi)
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

      return
      end
C------------  End of Subroutine ADRESWEB  ---------------------------------
      
C------------  Subroutine QRES  ---------------------------------
      subroutine QRES(T,Y,YP,Qsen,IRES,RPAR,IPAR,SENPAR)
      implicit double precision(a-h,o-z)
      double precision t, y(*), yp(*), rpar(*), senpar(*), qsen(*)
      integer ipar(*)
      integer ires, i
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
      

      qsen(1) = 0.0d0
      do i = 1, ns*my*mx
         qsen(1) = qsen(1) + y(i)*y(i)
c         if (ires .eq. 1) qsen(1+i) = 2.0d0*y(i)
      end do
      qsen(2) = qsen(1)
C
      return
      end
C------------  End of Subroutine QRES  ---------------------------------
      
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

!       IF (CALLS .EQ. 2) THEN
!          CALL PRINTIT('C:\\Documents and Settings\\STaylor\\My Documents
!      *\\DASPK3P1\\examples\\DM\\tamc_yadj1.mat',800,1,0,1,Y,800)
!          CALL PRINTIT('C:\\Documents and Settings\\STaylor\\My Documents
!      *\\DASPK3P1\\examples\\DM\\tamc_ady1.mat',800,1,0,1,ADY,800)
!       END IF
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
