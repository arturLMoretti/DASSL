c-----------------------------------------------------------------------
c
c     This file is part of the CWI Test Set for IVP solvers
c     http://www.cwi.nl/cwi/projects/IVPtestset/
c
c        Slider Crank, DAE format
c        index 2 DAE of dimension 24
c
c     DISCLAIMER: see
c     http://www.cwi.nl/cwi/projects/IVPtestset/disclaimer.html
c
c     The most recent version of this source file can be found at
c     http://www.cwi.nl/cwi/projects/IVPtestset/src/problems/crank.f
c
c     This is revision
c     $Id: crank.f,v 1.1 1998/11/25 14:19:30 walter Exp $
c
c-----------------------------------------------------------------------
c
c     These are the subroutines for the
c
c        SLIDER CRANK in index-2 formulation,
c        index 2 semi-explicit DAE of dimension 24
c
c     Integrator: RADAU5
c
c     Reference: Simeon, B.: Modelling a flexible slider crank
c          mechanism by a mixed system of DAEs and PDEs,
c          Math. Modelling of Systems 2, 1-18 (1996);
c          see also the IVPTestset report.
c
c     B. Simeon, 18.2.98
c
c------------------------------------------------------------------------
      integer function pidate()
      pidate = 19970704
      return
      end
c------------------------------------------------------------------------
      subroutine prob(fullnm,problm,type,
     +                neqn,ndisc,t,
     +                numjac,mljac,mujac,
     +                nummas,mlmas,mumas,
     +                ind)
      character*(*) fullnm,problm,type
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(*)
      double precision t(0:*)
      logical numjac, nummas

      integer i

      fullnm = 'Slider Crank'
      problm = 'crank'
      type   = 'DAE'
      neqn   = 24
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 0.1d0
      numjac = .true.
      mujac  = neqn
      mljac  = neqn
      mlmas  = 0
      mumas  = 0
      do 10 i=1,14
         ind(i) = 1
   10 continue
      do 20 i=15,24
         ind(i) = 2
   20 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime,consis)
      integer neqn
      double precision t,y(neqn),yprime(neqn)
      logical consis

      integer i
C
C     Initial values: 'Close' to smooth motion,
C     accelerations and multipliers consistent
C     for linear stiffness term and no damping
C     (ipar(1) = 0, ipar(2) = 0).
C
C     Position variables
C     phi1, phi2, x3, q1, q2, q3, q4
      y(1) = 0.d0
      y(2) = 0.d0
      y(3) = .450016933d+00
      y(4) = 0.d0
      y(5) = 0.d0
      y(6) = .103339863d-04
      y(7) = .169327969d-04
C     Initial values velocity variables
      y(8) =  .150000000d+03
      y(9) = -.749957670d+02
      y(10)= -.268938672d-05
      y(11)=  .444896105d+00
      y(12)=  .463434311d-02
      y(13)= -.178591076d-05
      y(14)= -.268938672d-05
C     Initial values acceleration variables
      y(15)= 0.d0
      y(16)= -1.344541576008661d-03
      y(17)= -5.062194923138079d+03
      y(18)= -6.833142732779555d-05
      y(19)=  1.449382650173157d-08
      y(20)= -4.268463211410861d+00
      y(21)=  2.098334687947376d-01
C     Lagrange multipliers
      y(22)= -6.397251492537153d-08
      y(23)=  3.824589508329281d+02
      y(24)= -4.376060460948886d-09
c
      do 10 i=1,24
         yprime(i) = 0.d0
  10  continue
      do 20 i=1,14
         yprime(i) = y(7+i)
  20  continue
c
      consis = .true.
      return
      end
      subroutine init2(neqn,t,y,yprime,consis)
      integer neqn
      double precision t,y(neqn),yprime(neqn)
      logical consis

      integer i
C
C     Initial values: On rigid motion trajectory,
C     leading to strong oscillations.
C     accelerations and multipliers consistent
C     for linear stiffness term and no damping
C     (ipar(1) = 0, ipar(2) = 0).
C
C     Position variables
C     phi1, phi2, x3, q1, q2, q3, q4
      y(1) = 0.d0
      y(2) = 0.d0
      y(3) = .4500d+00
      y(4) = 0.d0
      y(5) = 0.d0
      y(6) = 0.d0
      y(7) = 0.d0
C     Initial values velocity variables
      y(8) =  .150d+03
      y(9) = -.750d+02
      y(10)= 0.d0
      y(11)= 0.d0
      y(12)= 0.d0
      y(13)= 0.d0
      y(14)= 0.d0
C     Initial values acceleration variables
      y(15)= 0.d0
      y(16)= 0.d0
      y(17)= -3.789473684210526d+03
      y(18)= 0.d0
      y(19)= 0.d0
      y(20)=  1.924342105263158d+02
      y(21)=  1.273026315789474d+03
C     Lagrange multipliers
      y(22)= 0.d0
      y(23)= 2.863023157894737d+02
      y(24)= 0.d0
c
      do 10 i=1,24
         yprime(i) = 0.d0
  10  continue
      do 20 i=1,14
         yprime(i) = y(7+i)
  20  continue
c
      consis=.true.
      return
      end
      subroutine res(t,y,yprime,cj,delta,ires,rpar,ipar,senpar)
      implicit none
      integer ires,ipar(*), i
      double precision t,cj,y(*),yprime(*),delta(*),rpar(*),senpar(*)
      integer neqn
      neqn = 24
      call feval(neqn,t,y,yprime,delta,ires,rpar,ipar,senpar)
c
c...  for the index-2 constraints
      do i = 1,3
         delta(21+i) = cj*delta(21+i)
      end do
c
      return
      end
      subroutine feval(neqn,t,y,yprime,f,ierr,rpar,ipar2, senpar)
      integer neqn,ierr,ipar2(*)
      double precision t,y(neqn),yprime(neqn),f(neqn),rpar(*),senpar(*)

      integer i,ires, ityp, iequa, icall, ipar(2)
c
c     index-2 formulation:            ityp  = 1
c     evaluate all equations:         iequa = 0
c     new call:                       icall = 0
c     default - linear stiffness term ipar(1) = 0
c     default - no damping            ipar(2) = 0
c
      ipar(1)=0
      ipar(2)=0
      ityp  = 1
      iequa = 0
      icall = 0
c
c     evaluate residual
c     with derivatives dy set to zero
c     on return, fy contains the right
c     hand side up to a final multipl. by -1
c
      call resmbs(ityp, iequa, icall,
     *            t,y,yprime,f,ires,rpar,ipar,senpar)
      do 20 i=1,14
         f(i) = -f(i)
   20 continue
c
      end
c-----------------------------------------------------------------------
      subroutine jeval(ldim,neqn,t,y,yprime,dfdy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t, y(neqn), yprime(*),dfdy(ldim,neqn),rpar(*)
C
C ... Dummy subroutine for the Jacobian
C
      return
      end
c-----------------------------------------------------------------------
      subroutine meval(ldim,neqn,t,y,yprime,dfddy,ierr,rpar,ipar)
      integer    ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfddy(ldim,neqn),rpar(*)

      integer i
C
C     mass matrix in band storage
c
      do 10 i=1,14
            dfddy(1,i) = 1.d0
 10   continue
      do 20 i=15,24
         dfddy(1,i) = 0.d0
 20   continue
      return
      end
c-----------------------------------------------------------------------
      subroutine solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
c
c     MDOP5 applied to slider crank problem, tend = 0.1
c
c     ATOL = 1d-10 for p, v
c     RTOL = 1d-8  for p, v
c
c     # steps =             8411
c     # reject =             743
c     # f-eval =           50467
c     # projections =       2337
c
      y(1) = 1.500000000000000d+01
      y(2) = -3.311734987910881d-01
      y(3) = 1.697373326718410d-01
      y(4) = 1.893192460247178d-04
      y(5) = 2.375751865617931d-05
      y(6) = -5.323907988763734d-06
      y(7) = -8.363283141616840d-06
      y(8) = 1.500000000000000d+02
      y(9) = 6.025346682645789d+01
      y(10)= -8.753116989887888d+00
      y(11)= -3.005536801092212d-02
      y(12)= -5.500488291932075d-03
      y(13)= 4.978243404809343d-04
      y(14)= 1.104933470696396d-03
      y(15)= 0.d0
      y(16)= 6.488722210234531d+03
      y(17)= 2.167924253080623d+03
      y(18)= 3.391435115267547d+01
      y(19)= 1.699107480197843d-01
      y(20)= -1.415799354959001d+00
      y(21)= 9.903251655235532d-01
      y(22)= -6.232893262533717d+01
      y(23)= -1.637910131687472d+02
      y(24)= 2.529853213732781d+01
c
      return
      end

      SUBROUTINE RESMBS(ITYP, IEQUA, ICALL,
     *                  T,X,XD,DELTA,IRES,RPAR,IPAR,SENPAR)
      IMPLICIT   CHARACTER (A-Z)
      INTEGER    ITYP, IEQUA, ICALL, IRES, IPAR(2)

      DOUBLE PRECISION
     *           T, X(*), XD(*), DELTA(*), RPAR(*), SENPAR(*)
C
C     Slider crank - flexible model with sliding block
C     ------------------------------------------------
C  ** written by Bernd Simeon, TH Darmstadt, 13/06/95 **
C  ** extended version with special beam model,
C  ** for IVPTestset                         11/28/97 **
C
C     The flexible slider crank mechanism is described in
C
C     Simeon, B.: Modelling a Flexible Slider Crank Mechanism
C     by a Mixed System of DAEs and PDEs.
C     Math. Modelling of Systems 2, 1-18, 1995
C
C     This version contains all coupling terms (also for 2D)
C     for discretizations of the connecting rod.
C     Particular grid used: 2 modal functions lateral;
C                           2 nodes (quadratic ansatz) longitudinal.
C
C     phi1(t) = omega*t prescribed by constraint.
C
C     PARAMETERS ON ENTRY:
C
C     ITYP     This integer flag determines in which formulation the equations
C              of motion have to be evaluated:
C              ITYP = 0: index 3 system;     ITYP = 3: index 1 system;
C                   = 1: index 2 system;
C
C     IEQUA    This integer flag determines whether the complete residual
C              or only parts of it have to be evaluated.
C              IEQUA = 0: evaluate complete residual;
C                    = 2: evaluate only position+velocity constraints
C                         in DELTA(1:6)
C
C     ICALL    This integer flag indicates whether RESMBS has already
C              been called with the actual parameter set T, X, XD.
C              ICALL = 0: new values T, X, XD;
C                    = 1: T, X, XD have not changed since the last call.
C              (unused)
C
C     T        This real variable contains the current value of the
C              independent variable (time).
C
C     X(*)     This array contains the current values of the dependent
C              variables: The multibody system variables are arranged as
C              X = [ p                   rigid motion (3 coordinates)
C                    q                   elastic motion (4 nodes)
C                    pd                  velocity variables p
C                    qd                  velocity variables q
C                    w                   acceleration variables p, q
C                    lambda              3 Lagrange multipliers
C
C     XD(*)    This array contains the derivatives of the solution
C              components at T.
C
C     RPAR,IPAR These are real and integer parameter arrays which
C              are used for communication between the calling program
C              and this subroutine.
C       IPAR(1)  0: only linear stiffness term K*q
C                1: nonlinear stiffness term included
C
C       IPAR(2)  0: no physical damping - purely imaginary EV
C                1: damping matrix (1 %) included
C
C     ON RETURN:
C
C     DELTA(NEQ) This array contains the residual of the equations of motion.
C
C     IRES     This integer flag indicates a stop condition (unused).
C
C     SYSTEM PARAMETERS:
C
C       IPAR(1)  0: only linear stiffness term K*q
C                1: nonlinear stiffness term included
C
C       IPAR(2)  0: no physical damping - purely imaginary EV
C                1: damping matrix (0.5 %) included
C
        INTEGER  KU, KV
        DOUBLE PRECISION
     *           GRAV, OMEGA, T1, J1, J2, L1, L2, M1, M2, M3, PI,
     *           EE, NUE, BB, HH, RHO
C
C     Data set
C
        PARAMETER( M1 = 0.36D0,     M2 = 0.151104D0,
     *             M3 = 0.075552D0,
     *             J1 = 0.002727D0, J2 = 0.0045339259D0,
     *             PI = 3.1415927D0,
     *             EE = .20D12,     NUE= 0.30D0,
     *             BB = 0.0080D0,   HH = 0.0080D0,
     *             RHO= 7870.0D0,
     *             GRAV= 0.0D0, OMEGA = 150.D0          )
c     *             L1 = 0.15D0,     L2 = 0.30D0,
C
C     LOCAL VARIABLES:
C
C       Q, QD for FE coefficients and time derivatives,
C       MQ, KQ, DQ, BQ, c1, c2, c12, c21 for FE matrices and vectors,
C       up to NQMAX = 20 variables.
C
        INTEGER  NQMAX, I, J, JJ, NQ, NP, NL, NX
        PARAMETER( NQMAX = 20 )
        DOUBLE PRECISION
     *           Q(NQMAX), QD(NQMAX), MQ(NQMAX,NQMAX),
     *           KQ(NQMAX,NQMAX), BQ(NQMAX,NQMAX), DQ(NQMAX,NQMAX),
     *           c1(NQMAX), c2(NQMAX), c12(NQMAX), c21(NQMAX),
     *           MQQ(NQMAX), KQQ(NQMAX), DQQD(NQMAX),
     *           QTBQ(NQMAX), BQQD(NQMAX),
     *           c1TQ, c1TQD, c2TQ, c2TQD, c12TQ, c12TQD,
     *           QDTBQQD, QTMQQ, QDTMQQ, DDOT, V(2),
     *           ALC(3), PLC(3), VLC(3),
     *           AM(NQMAX+3,NQMAX+3), GP(3,NQMAX+3), F(NQMAX+3),
     *           COSP1, COSP2, SINP1, SINP2, COSP12, SINP12,
     *           QKU, QKV, QDKU, QDKV, FACM, FACK, FACB
        SAVE     MQ, KQ, DQ, BQ, c1, c2, c12, c21
C
C       FIRST for first call - evaluation of FE matrices.
C
        LOGICAL  FIRST
        DATA     FIRST / .TRUE. /
C
C_______________End of declaration part RESMBS____________________________
C
C
C...  INPUT FROM senpar
c
      L1 = senpar(1)
      L2 = senpar(2)

      IRES  = 0
      NQ    = 4
      NP    = 7
      NL    = 3
      NX    = 3*NP + NL
      KU    = 4
      KV    = 0
      IF (FIRST) THEN
C
C       Initialize grid data.
C
        FACM = RHO*BB*HH*L2
        FACK = EE*BB*HH/L2
        FACB = BB*HH*L2
C
        DO 5 I=1,NQ
           DO 4 J=1,NQ
              MQ(J,I) = 0.D0
              KQ(J,I) = 0.D0
              BQ(J,I) = 0.D0
              DQ(J,I) = 0.D0
 4         CONTINUE
           c1(I) = 0.D0
           c2(I) = 0.D0
           c12(I)= 0.D0
           c21(I)= 0.D0
 5      CONTINUE
C
        MQ(1,1) = FACM*.5D0
        MQ(2,2) = FACM*.5D0
        MQ(3,3) = FACM*8.D0
        MQ(3,4) = FACM*1.D0
        MQ(4,3) = FACM*1.D0
        MQ(4,4) = FACM*2.D0
C
        KQ(1,1) = FACK*PI**4/24.D0*(HH/L2)**2
        KQ(2,2) = FACK*PI**4*2.D0/3.D0*(HH/L2)**2
        KQ(3,3) = FACK*16.D0/3.D0
        KQ(3,4) = -FACK*8.D0/3.D0
        KQ(4,3) = -FACK*8.D0/3.D0
        KQ(4,4) = FACK*7.D0/3.D0
C
        BQ(1,3) = -FACB*16.D0/PI**3
        BQ(1,4) = FACB*(8.D0/PI**3-1.D0/PI)
        BQ(2,4) = FACB*0.5D0/PI
        BQ(3,1) = FACB*16.D0/PI**3
        BQ(4,1) = -FACB*(8.D0/PI**3-1.D0/PI)
        BQ(4,2) = -FACB*0.5D0/PI
C
        c1(3)  = FACB*2.D0/3.D0
        c1(4)  = FACB*1.D0/6.D0
        c2(1)  = FACB*2.D0/PI
        c12(3) = L2*FACB*1.D0/3.D0
        c12(4) = L2*FACB*1.D0/6.D0
        c21(1) = L2*FACB*1.D0/PI
        c21(2) = -L2*FACB*0.5D0/PI
C
        IF (IPAR(2) .EQ. 1) THEN
C
C       0.5 per cent damping
C
           DQ(1,1) = 5.D0
           DQ(2,2) = 25.D0
           DQ(3,3) = 0.5D0*2.308375455264791D+02
           DQ(3,4) = -0.5D0*2.62688487992052D+02
           DQ(4,3) = -0.5D0*2.626884879920526D+02
           DQ(4,4) = 0.5D0*4.217421837156818D+02
        END IF
        FIRST = .FALSE.
      END IF
C
      COSP1  = COS(X(1))
      COSP2  = COS(X(2))
      SINP1  = SIN(X(1))
      SINP2  = SIN(X(2))
      COSP12 = COS(X(1)-X(2))
      SINP12 = SIN(X(1)-X(2))
      V(1)   = X(NP+1)
      V(2)   = X(NP+2)
C
      DO 6 I=1,NQ
         Q(I)  = X(3+I)
         QD(I) = X(NP+3+I)
   6  CONTINUE
C
C     Evaluate scalar products and quadratic forms.
C
      c1TQ  = DDOT(NQ,c1,1,Q,1)
      c1TQD = DDOT(NQ,c1,1,QD,1)
      c2TQ  = DDOT(NQ,c2,1,Q,1)
      c2TQD = DDOT(NQ,c2,1,QD,1)
      c12TQ = DDOT(NQ,c12,1,Q,1)
      c12TQD= DDOT(NQ,c12,1,QD,1)
      DO 10 I=1,NQ
         MQQ(I) = DDOT(NQ,MQ(1,I),1,Q,1)
         KQQ(I) = DDOT(NQ,KQ(1,I),1,Q,1)
         DQQD(I)= DDOT(NQ,DQ(1,I),1,QD,1)
         QTBQ(I)= DDOT(NQ,Q,1,BQ(1,I),1)
         BQQD(I)= DDOT(NQ,BQ(I,1),NQMAX,QD,1)
  10  CONTINUE
      QTMQQ   = DDOT(NQ,Q,1,MQQ,1)
      QDTMQQ  = DDOT(NQ,QD,1,MQQ,1)
      QDTBQQD = DDOT(NQ,QD,1,BQQD,1)
C
C     Kinematic and dynamic equations.
C
      DO 50 I=1,NP
         DELTA(I)    = XD(I)    - X(NP+I)
         DELTA(NP+I) = XD(NP+I) - X(2*NP+I)
 50   CONTINUE
C
C     Compute mass matrix.
C
          AM(1,1) = J1 + M2*L1*L1
          AM(1,2) = .5D0*L1*L2*M2*COSP12
          AM(2,2) = J2
          AM(1,3) = 0.D0
          AM(2,3) = 0.D0
          AM(3,1) = 0.D0
          AM(3,2) = 0.D0
          AM(3,3) = M3
          AM(1,2) = AM(1,2) + RHO*L1*(SINP12*c2TQ+COSP12*c1TQ)
          AM(2,2) = AM(2,2) + QTMQQ + 2.0D0*RHO*c12TQ
          DO 100 I=1,NQ
             AM(1,3+I) = RHO*L1*(-SINP12*c1(I) + COSP12*c2(I))
             AM(2,3+I) = RHO*c21(I) + RHO*QTBQ(I)
             AM(3,3+I) = 0.D0
  100     CONTINUE
          DO 120 I=1,NQ
             DO 110 J=1,I
                AM(3+J,3+I) = MQ(J,I)
  110        CONTINUE
  120     CONTINUE
          DO 140 I=1,NP
             DO 130 J=I+1,NP
                AM(J,I) = AM(I,J)
  130        CONTINUE
  140     CONTINUE
C
C     Compute constraint matrix.
C
          IF (KU .EQ. 0) THEN
             QKU = 0.D0
          ELSE
             QKU = Q(KU)
          END IF
          IF (KV .EQ. 0) THEN
             QKV = 0.D0
          ELSE
             QKV = Q(KV)
          END IF
          GP(1,1) = L1*COSP1
          GP(1,2) = L2*COSP2 + QKU*COSP2 - QKV*SINP2
          GP(1,3) = 0.D0
          GP(2,1) = L1*SINP1
          GP(2,2) = L2*SINP2 + QKU*SINP2 + QKV*COSP2
          GP(2,3) = 1.D0
          GP(3,1) = 1.D0
          GP(3,2) = 0.D0
          GP(3,3) = 0.D0
          DO 150 I=1,NQ
             GP(1,3+I) = 0.D0
             GP(2,3+I) = 0.D0
             GP(3,3+I) = 0.D0
  150     CONTINUE
          IF (KU .NE. 0) THEN
             GP(1,3+KU) = SINP2
             GP(2,3+KU) = -COSP2
          END IF
          IF (KV .NE. 0) THEN
             GP(1,3+KV) = COSP2
             GP(2,3+KV) = SINP2
          END IF
C
C     Forces - rigid motion entries.
C
          F(1) = -.5D0*L1*GRAV*(M1+2.0D0*M2)*COSP1
     &           -.5D0*L1*L2*M2*V(2)*V(2)*SINP12
c          IF (T .LE. T1) THEN
c             F(1) =  F(1)+OMEGA/T1*(1.D0-COS(2.D0*PI*T/T1))
c          END IF
          F(2) = -.5D0*L2*GRAV*M2*COSP2
     &           +.5D0*L1*L2*M2*V(1)*V(1)*SINP12
          F(3) = 0.d0
C
C     Superposition of flexible motion (term f^e).
C
          F(1) = F(1)
     &         + RHO*L1*V(2)*V(2)*(-SINP12*c1TQ+COSP12*c2TQ)
     &         - 2.0D0*RHO*L1*V(2)*(COSP12*c1TQD+SINP12*c2TQD)
          F(2) = F(2)
     &         + RHO*L1*V(1)*V(1)*(SINP12*c1TQ-COSP12*c2TQ)
     &         - 2.0D0*RHO*V(2)*c12TQD - 2.0D0*V(2)*QDTMQQ
     &         - RHO*QDTBQQD - RHO*GRAV*(COSP2*c1TQ-SINP2*c2TQ)
C
C     Coriolis and gravity terms flexible motion (Gamma).
C
          DO 200 I=1,NQ
             F(3+I) = V(2)*V(2)*MQQ(I)
     &         + RHO*(V(2)*V(2)*c12(I)
     &                + L1*V(1)*V(1)*(COSP12*c1(I)+SINP12*c2(I))
     &                + 2.0D0*V(2)*BQQD(I) )
     &         - RHO*GRAV*(SINP2*c1(I)+COSP2*c2(I))
  200     CONTINUE
C
C         Stiffness + damping terms - K q - D q'.
C
          DO 210 I=1,NQ
             F(3+I) = F(3+I) - KQQ(I) - DQQD(I)
  210     CONTINUE
          IF (IPAR(1) .EQ. 1) THEN
C
C            Nonlinear stiffness term
C
             FACK = 0.5D0*EE*BB*HH/L2**2*PI**2
             FACB = 80.D0/(PI**2*9.D0)
             F(4) = F(4) -
     &              FACK*(Q(1)*Q(4)-FACB*Q(2)*(-4*Q(3)+2*Q(4)))
             F(5) = F(5) -
     &              FACK*(4*Q(2)*Q(4)-FACB*Q(1)*(-4*Q(3)+2*Q(4)))
             F(6) = F(6) -
     &              FACK*4.D0*FACB*Q(1)*Q(2)
             F(7) = F(7) -
     &              FACK*(0.5D0*Q(1)**2+2*Q(2)**2-2*FACB*Q(1)*Q(2))
          END IF
C
C     Dynamics part II ( M*w - f + G(T)*lambda ).
C
      DO 250 I=1,NP
         DELTA(2*NP+I) = DDOT(NP,AM(1,I),1,X(2*NP+1),1)
     & - F(I) + GP(1,I)*X(NX-2)+GP(2,I)*X(NX-1)+GP(3,I)*X(NX)
  250 CONTINUE
C
C     Acceleration level constraints.
C
      IF (KU .EQ. 0) THEN
          QDKU = 0.D0
      ELSE
          QDKU = QD(KU)
      END IF
      IF (KV .EQ. 0) THEN
          QDKV = 0.D0
      ELSE
          QDKV = QD(KV)
      END IF
      ALC(1) = -L1*SINP1*V(1)*V(1) - (L2+QKU)*SINP2*V(2)*V(2)
     *   +2.0D0*V(2)*(COSP2*QDKU-SINP2*QDKV) - COSP2*V(2)*V(2)*QKV
      ALC(2) =  L1*COSP1*V(1)*V(1) + (L2+QKU)*COSP2*V(2)*V(2)
     *   +2.0D0*V(2)*(SINP2*QDKU+COSP2*QDKV) - SINP2*V(2)*V(2)*QKV
      ALC(3) = 0.0D0
      DO 300 I=1,NP
         ALC(1) = ALC(1) + GP(1,I)*X(2*NP+I)
         ALC(2) = ALC(2) + GP(2,I)*X(2*NP+I)
         ALC(3) = ALC(3) + GP(3,I)*X(2*NP+I)
  300 CONTINUE
C
C     Position level constraints.
C
      PLC(1) = L1*SINP1 + L2*SINP2 + QKU*SINP2 + QKV*COSP2
      PLC(2) = X(3) - L1*COSP1 - L2*COSP2
     &         -QKU*COSP2 + QKV*SINP2
      PLC(3) = X(1) - OMEGA*T
C
C     Velocity level constraints.
C
      VLC(1) = 0.0D0
      VLC(2) = 0.0D0
      VLC(3) = -OMEGA
      DO 400 I=1,NP
         VLC(1) = VLC(1) + GP(1,I)*X(NP+I)
         VLC(2) = VLC(2) + GP(2,I)*X(NP+I)
         VLC(3) = VLC(3) + GP(3,I)*X(NP+I)
  400 CONTINUE
C
      IF (IEQUA .EQ. 2) THEN
C
C         Evaluate only the constraints.
C
          DELTA(1) = PLC(1)
          DELTA(2) = PLC(2)
          DELTA(3) = PLC(3)
          DELTA(4) = VLC(1)
          DELTA(5) = VLC(2)
          DELTA(6) = VLC(3)
      ELSE
C
C         Select constraints defined by ITYP.
C
          IF (ITYP .EQ. 0) THEN
C
C             Index 3 system.
C
              DELTA(NX-2) = PLC(1)
              DELTA(NX-1) = PLC(2)
              DELTA(NX)   = PLC(3)
          ELSE IF (ITYP .EQ. 1) THEN
C
C             Index 2 system.
C
              DELTA(NX-2) = VLC(1)
              DELTA(NX-1) = VLC(2)
              DELTA(NX)   = VLC(3)
          ELSE IF (ITYP .EQ. 3) THEN
C
C             Index 1 system.
C
              DELTA(NX-2) = ALC(1)
              DELTA(NX-1) = ALC(2)
              DELTA(NX)   = ALC(3)
          END IF
      END IF
C
C_______________End of subroutine RESMBS____________________________
C
      RETURN
      END
c$$$      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
c$$$C
c$$$C     FORMS THE DOT PRODUCT OF TWO VECTORS.
c$$$C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
c$$$C     JACK DONGARRA, LINPACK, 3/11/78.
c$$$C
c$$$      DOUBLE PRECISION DX(1),DY(1),DTEMP
c$$$      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
c$$$C
c$$$      DDOT = 0.0D0
c$$$      DTEMP = 0.0D0
c$$$      IF(N.LE.0)RETURN
c$$$      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
c$$$C
c$$$C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
c$$$C          NOT EQUAL TO 1
c$$$C
c$$$      IX = 1
c$$$      IY = 1
c$$$      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
c$$$      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
c$$$      DO 10 I = 1,N
c$$$        DTEMP = DTEMP + DX(IX)*DY(IY)
c$$$        IX = IX + INCX
c$$$        IY = IY + INCY
c$$$   10 CONTINUE
c$$$      DDOT = DTEMP
c$$$      RETURN
c$$$C
c$$$C        CODE FOR BOTH INCREMENTS EQUAL TO 1
c$$$C
c$$$C
c$$$C        CLEAN-UP LOOP
c$$$C
c$$$   20 M = MOD(N,5)
c$$$      IF( M .EQ. 0 ) GO TO 40
c$$$      DO 30 I = 1,M
c$$$        DTEMP = DTEMP + DX(I)*DY(I)
c$$$   30 CONTINUE
c$$$      IF( N .LT. 5 ) GO TO 60
c$$$   40 MP1 = M + 1
c$$$      DO 50 I = MP1,N,5
c$$$        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
c$$$     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
c$$$   50 CONTINUE
c$$$   60 DDOT = DTEMP
c$$$      RETURN
c$$$      END
