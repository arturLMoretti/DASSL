      SUBROUTINE AD_DJACILU (ADRES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
     .     WK, H, CJ, WP, IWP, IERR, RPAR, IPAR, SENPAR, ADI, 
     .     ISENFO,SENWRK,ISENWK,G_RES, P_RES)
C ... Version of 10-6-95

C ... Calculate ILU decomposition of the Jacobian matrix
C     for use as a preconditioner by the DDASPK solver.

      IMPLICIT NONE

C ... Input arguments:
      INTEGER NEQ        ! total number of equations
      REAL*8 T           ! independent variable t
      REAL*8 Y(NEQ)      ! most recent iterate of solution vector y
      REAL*8 YPRIME(NEQ) ! most recent iterate of solution vector y'
      REAL*8 SAVR(NEQ)   ! current residual evaluated at (T,Y,YPRIME)
      REAL*8 REWT(NEQ)   ! scale factors for Y and YPRIME
      INTEGER IRES       ! error flag for RES routine
      REAL*8 WK(NEQ)     ! work space available to this subroutine
      REAL*8 H           ! current step size
      REAL*8 CJ          ! scalar proportional to 1/H
      REAL*8 RPAR(*)     ! user real workspace
      INTEGER IPAR(*)    ! user integer workspace
      real*8 SENPAR(*)  ! sensitivity parameter array

C ... Output arguments:
      REAL*8 WP(*)       ! matrix elements of ILU
      INTEGER IWP(*)     ! array indices for elements of ILU
      INTEGER NRE        ! number of RES calls needed to evaluate
                         ! Jacobian NRE is returned in IPAR(30)
      REAL*8 ADI(*)      
      INTEGER IERR       ! error flag (0 means success, else failure)
      REAL*8  SENWRK(*)
      INTEGER ISENFO(*), ISENWK(*)
      EXTERNAL ADRES, G_RES, P_RES

C ... Local variables:
      REAL*8 TOLILUT, PERMTOL, SQRTN
      INTEGER I, LBW, UBW, LENPLUMX, LJAC, LJACI, LJACJ, LIPERM,
     .        LROWNRMS, LRWK1, LIWK1, IFMT
      INTEGER LENPFAC, LENPLUFAC, LFILILUT, IPREMETH, NEQP1, NNZMX
      INTEGER LPLU, LJU, LJLU, LPERM, LQPERM, LLEVELS, LMASK
      INTEGER ISRNORM  !=1 causes row normalization of Jac.
      INTEGER NORMTYPE !=0,1,2 for max-norm, 1-norm, or
                       ! 2-norm row scaling
      INTEGER IREORDER !=1 causes row and column reordering of Jac.
      INTEGER JACOUT   !=1 causes the Jacobian matrix and SAVR to
                       !   be written to a file and then exit with
                       !   IRES = -2 to signal a stop to DDASPK.
      INTEGER IUNIT    ! logical unit number to use when JACOUT .EQ. 1
      INTEGER JSCALCOL !=1 causes the columns of the Jacobian matrix
                       !   to be scaled by REWT-inverse
      CHARACTER*8 PMETH(4), PREMETH
      CHARACTER*72 TITLE
      CHARACTER*80 MSG
      SAVE
      DATA PMETH/'ILUT','ILUTP','ILU0','MILU0'/

C ... Zero out NRE counter
      NRE = 0

C ... Load values from IPAR and RPAR.
      LBW = IPAR(1)
      UBW = IPAR(2)
      LENPFAC = IPAR(3)
      LENPLUFAC = IPAR(4)
      IPREMETH = IPAR(5)
      LFILILUT = IPAR(6)
      IREORDER = IPAR(7)
      ISRNORM = IPAR(8)
      NORMTYPE = IPAR(9)
      JACOUT = IPAR(10)
      JSCALCOL = IPAR(11)
      TOLILUT = RPAR(1)
      PERMTOL = RPAR(2)
      PREMETH = PMETH(IPREMETH)

C...  Set pointers into the WP and IWP arrays.
      NEQP1 = NEQ + 1
      NNZMX = LENPFAC*NEQ
      LENPLUMX = NNZMX + LENPLUFAC*NEQ
C ... Set up pointers into WP
      LJAC = 1
      LROWNRMS = NNZMX + LJAC
      IF (ISRNORM .EQ. 1) THEN
         LPLU = LROWNRMS + NEQ
      ELSE
         LPLU = LROWNRMS
      ENDIF
      LRWK1 = LPLU + LENPLUMX
C ... Set up pointers into IWP
      LJACI = 1
      LJACJ = LJACI + NEQP1
      LJU = LJACJ + NNZMX
      LJLU = LJU + LENPLUMX

C ... Calculate Jacobian matrix.
      IERR = 0
*      CALL DJCALC (NEQ, T, Y, YPRIME, SAVR, LBW, UBW, WK, REWT, RES,
*      CALL AD_JAC (NEQ, T, Y, YPRIME, SAVR, LBW, UBW, WK, REWT, RES,
      CALL AD_DJCALC (NEQ, T, Y, YPRIME, SAVR, LBW, UBW, WK, REWT,adres,
     .             H, CJ, NNZMX, WP(LJAC), IWP(LJACJ), IWP(LJACI),
     .             WP(LPLU), IWP(LJLU), IWP(LJU), IPAR, RPAR,
     .             SENPAR, IRES, NRE, ADI, IERR,
     .             ISENFO,SENWRK,ISENWK,G_RES, P_RES)
      IF (IRES .LT. 0) RETURN
      IF (IERR .NE. 0) RETURN

C ... Save NRE value for user output.
      IPAR(30) = IPAR(30) + NRE
      
C ... Modify pointers into IWP
      LJLU = LJU + NEQP1
      IF (IREORDER .NE. 0) THEN
         LPERM = LJLU + LENPLUMX
         LQPERM = LPERM + NEQ
         LIWK1 = LQPERM + NEQ
         LLEVELS = LJLU + NNZMX   ! assumes that LENPLUFAC >= 2.
         LMASK = LLEVELS + NEQ
      ELSE
         LPERM = 0
         LQPERM = 0
         LLEVELS = 0
         LMASK = 0
         LIWK1 = LJLU + LENPLUMX
      ENDIF
      IF (PREMETH .EQ. 'ILUTP') THEN
         LIPERM = LIWK1 + 2*NEQ
      ELSE
         LIPERM = LIWK1
      ENDIF
C ... Multiply Jacobian columns by inverse of scaling vector REWT.
C     In PSOLILU, the WGHT array equals REWT/SQRT(NEQ), so we must
C     be consistent here.
      IF (JSCALCOL .EQ. 1) THEN
         SQRTN = SQRT(REAL(NEQ))
         DO 10 I = 1, NEQ
            WK(I) = SQRTN / REWT(I)
 10      CONTINUE
         CALL AMUDIA (NEQ, 0, WP(LJAC), IWP(LJACJ), IWP(LJACI), WK,
     .                WP(LJAC), IWP(LJACJ), IWP(LJACI))
      ENDIF

C ... Normalize Jacobian rows, if desired.
      IF (ISRNORM .EQ. 1) THEN
         CALL ROSCAL (NEQ,0,NORMTYPE,WP(LJAC),IWP(LJACJ),IWP(LJACI),
     .        WP(LROWNRMS),WP(LJAC),IWP(LJACJ),IWP(LJACI), IERR)
         IF (IERR .NE. 0) RETURN
      ENDIF

C ... Reorder Jacobian rows and columns, if desired.
      IF (IREORDER .NE. 0) THEN
         CALL DJREORD (NEQ, NEQP1, NNZMX, PREMETH,
     .                 WP(LJAC), IWP(LJACJ), IWP(LJACI),
     .                 WP(LPLU), IWP(LJLU), IWP(LJU),
     .                 IWP(LPERM), IWP(LQPERM), IWP(LLEVELS), 
     .                 IWP(LMASK), IREORDER)
      ENDIF

C ... Write matrix JAC and scaled RES to file if JACOUT .eq. 1.
      IF (JACOUT .EQ. 1) THEN
         IUNIT = IPAR(29)
         IF (ISRNORM .EQ. 1) THEN
            DO 20 I = 1, NEQ
               SAVR(I) = SAVR(I) * WP(LROWNRMS+I-1)
 20         CONTINUE
         ENDIF
         IF (IREORDER .NE. 0) CALL DVPERM (NEQ, SAVR, IWP(LPERM))
         TITLE = ' DDASPK Test Matrix '
         IFMT = 15
         CALL PRTMT (NEQ,NEQ,WP(LJAC),IWP(LJACJ),IWP(LJACI),SAVR,
     .        'NN',TITLE,'SPARSKIT','RUA',IFMT,3,IUNIT)
         MSG = 'DJACILU -- Jacobian Matrix written to file.'
         CALL XERRWD(MSG,80,0,0,0,0,0,0,0.0,0.0)
         IERR = 1
         IRES = -2
         RETURN
      ENDIF

C ... Compute ILU decomposition.
      CALL DJILU (NEQ, NEQ+1, NNZMX, WP(LJAC), IWP(LJACJ), 
     .            IWP(LJACI), IWP(LJU), WP(LPLU), IWP(LJLU),
     .            WP(LRWK1), IWP(LIWK1), LENPLUMX, TOLILUT,
     .            LFILILUT, PERMTOL, PREMETH, IWP(LIPERM), IERR)
      IF ((IERR .EQ. -2) .OR. (IERR .EQ. -3)) THEN
         IRES = -2  ! Stop since more storage needed.
      ENDIF

C ... Save pointers for use in DPSOLILU into IPAR array.
      IPAR(21) = LPLU
      IPAR(22) = LJU
      IPAR(23) = LJLU
      IPAR(24) = LROWNRMS
      IPAR(25) = LPERM
      IPAR(26) = LQPERM

      RETURN
C------------  End of Subroutine DJACILU  ------------------------------
      END

      SUBROUTINE AD_DJCALC (NEQ, T, Y, YPRIME, R0, ML,MU,R1,REWT,adRES,
     .                   H, CJ,  NNZMX, JAC, JA, IA, RCOO, JCOO, ICOO,
     .                   IPAR, RPAR, SENPAR,IRES, NRE, ADI, IERR,
     .                   ISENFO,SENWRK,ISENWK,G_RES, P_RES)

c ... Version of 10-6-95

C ... Calculate Jacobian matrix (derivatives with respect to each
C     dependent variable of the right-hand side of each rate equation).
C     Lower and upper bandwidths are used to select for computation
C     only those Jacobian elements that may be nonzero.
C     Estimates of Jacobian elements are computed by finite differences.
C     The Jacobian is stored in compressed sparse row format.

      IMPLICIT NONE

C ... Input arguments:
      INTEGER NEQ        ! total number of equations
      REAL*8 T           ! independent variable t
      REAL*8 Y(NEQ)      ! most recent iterate of solution vector y
      REAL*8 YPRIME(NEQ) ! most recent iterate of solution vector y'
      REAL*8 R0(NEQ)     ! current residual evaluated at (T,Y,YPRIME)
      REAL*8 REWT(NEQ)   ! array of scaling factors for Y and YPRIME
      INTEGER IRES       ! error flag for RES routine
      INTEGER ML, MU     ! lower and upper bandwidths
      INTEGER NNZMX      ! maximum number of nonzeros in Jacobian
      REAL*8 H           ! current step size
      REAL*8 CJ          ! scalar proportional to 1/H
      REAL*8 RPAR(*)     ! user real workspace
      INTEGER IPAR(*)    ! user integer workspace
      real*8 SENPAR(*)  ! sensitivity parameter array

C ... Work-array argument:
      REAL*8 R1(NEQ)   ! work space available to this subroutine

C ... Output arguments:
      REAL*8 JAC(NNZMX)   ! nonzero Jacobian elements
      INTEGER JA(NNZMX)   ! col indices of nonzero Jacobian elements
      INTEGER IA(NEQ+1)   ! pointers to beginning of each row in JAC,JA
      INTEGER IERR

C ... Workspace for temporary storage of Jacobian elements:
      REAL*8 RCOO(NNZMX)   ! nonzero Jacobian elements
      INTEGER JCOO(NNZMX)  ! col indices of nonzero Jacobian elements
      INTEGER ICOO(NNZMX)  ! row indices of nonzero Jacobian elements
      REAL*8 ADI(*)
      REAL*8  SENWRK(*)
      INTEGER ISENFO(*), ISENWK(*)
      EXTERNAL ADRES, G_RES, P_RES

C ... Local variables:
      INTEGER NNZ, I, I1, I2, J, JJ, MBA, MEBAND, MEB1, MBAND, NRE,NQ
      integer NY
      REAL*8 JACELEM, UROUND, D1MACH, SQUR, DEL, DELINV
      CHARACTER*80 MSG

C ... Set band parameters.
      NNZ = 1
      MBAND = ML + MU + 1
      MBA = MIN0(MBAND,NEQ)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1

C ... Set the machine unit roundoff UROUND and SQRT(UROUND), used to 
C ... set increments in the difference quotient procedure. 
      UROUND = D1MACH(4)
      SQUR = SQRT(UROUND)

C ... Initial error flags.
      IERR = 0
      IRES = 0

C ... Make MBA calls to RES to approximate the Jacobian.
C ... Here, R0(1),...,R0(neq) contains the base RES value, and 
C     R1(1),...,R1(NEQ) contains the perturbed values of RES.
c
      NY = NEQ + ISENFO(10)
      DO 40 J = 1,MBA
        DO 10 JJ = J,NEQ,MBAND
          JAC(JJ) = Y(JJ)
          JAC(JJ+NEQ) = YPRIME(JJ)
          DEL = SQUR*MAX(ABS(Y(JJ)),ABS(H*YPRIME(JJ)),ABS(1.0/REWT(JJ)))
          DEL = SIGN(DEL, H*YPRIME(JJ))
          DEL = (Y(JJ) + DEL) - Y(JJ)
          Y(JJ) = Y(JJ) + DEL
          YPRIME(JJ) = YPRIME(JJ) + CJ*DEL
 10       CONTINUE
          CALL DDRESAD(
     *         T, Y, YPRIME, CJ, R1, IRES, RPAR, IPAR, SENPAR,
     *         adres, NY, isenfo, SENWRK, ISENWK, g_res, p_res, adi)
        IF (IRES .LT. 0) RETURN
        NRE = NRE + 1
        DO 30 JJ = J,NEQ,MBAND
          Y(JJ) = JAC(JJ)
          YPRIME(JJ) = JAC(JJ+NEQ)
          DEL = SQUR*MAX(ABS(Y(JJ)),ABS(H*YPRIME(JJ)),ABS(1.0/REWT(JJ)))
          DEL = SIGN(DEL, H*YPRIME(JJ))
          DEL = (Y(JJ) + DEL) - Y(JJ)
          DELINV=1.0/DEL
          I1 = MAX(1,(JJ-MU))
          I2 = MIN(NEQ,(JJ+ML))
          DO 20 I = I1,I2
C ... Calculate possibly nonzero Jacobian elements for this variable,
C     and store nonzero elements in coordinate format.
            JACELEM = (R1(I) - R0(I))*DELINV
            IF (JACELEM .NE. 0.) THEN
               IF (NNZ .GT. NNZMX) THEN
                  MSG = 'DJCALC -- More storage needed for Jacobian.'
                  CALL XERRWD(MSG,80,0,0,0,0,0,0,0.0,0.0)
                  MSG = 'DJCALC -- Increase LENPFAC.'
                  CALL XERRWD(MSG,80,0,0,0,0,0,0,0.0,0.0)
                  MSG = 'DJCALC -- Storage exceeded at (I,J) = (I1,I2)'
                  CALL XERRWD(MSG,80,0,0,2,I,JJ,0,0.0,0.0)
                  IERR = 1
                  IRES = -2
	          RETURN
               ENDIF
               RCOO(NNZ) = JACELEM
               JCOO(NNZ) = JJ
               ICOO(NNZ) = I
               NNZ = NNZ + 1
            ENDIF
 20         CONTINUE
 30       CONTINUE
 40    CONTINUE
      NNZ = NNZ - 1
C
C ... Convert Jacobian from coordinate to compressed sparse row format.
      CALL COOCSR (NEQ, NNZ, RCOO, ICOO, JCOO, JAC, JA, IA)

      RETURN
C------------  End of Subroutine DJCALC  -------------------------------
      END
