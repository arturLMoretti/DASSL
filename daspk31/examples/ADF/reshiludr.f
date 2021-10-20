      SUBROUTINE QRES (T,U,UPRIME,QSEN,ires,RPAR,IPAR,senpar)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), RPAR(*), IPAR(*),senpar(*), Qsen(*)
C
C Set problem constants using IPAR and RPAR.
      NEQ = IPAR(33)
      M = IPAR(34)
      M2 = M + 2
C
C Loop over interior points, and load residual values.
      qsen(1) = 0.0d0
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
           I = IOFF + J + 1
c           if (u(i) .gt. qsen) qsen = u(i)
           qsen(1) = qsen(1) + u(i)*u(i)
 20     CONTINUE
 30   CONTINUE
C
      RETURN
      END
