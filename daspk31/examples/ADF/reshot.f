      SUBROUTINE RES (T, U, UPRIME, CJ, DELTA,IRES,RPAR,IPAR,SENPAR)
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
      damk = dexp(20.0d0)/(20.0d0)
      DO 30 K = 1,M
        IOFF = M2*K
        DO 20 J = 1,M
           I = IOFF + J + 1
           etmp = dexp(-20.0d0/u(I))
           ftmp = senpar(1)*damk*(2.0d0 - u(I))*etmp
           TEMX = 0.1d0*(U(I-1)  + U(I+1)  - 2.0d0*u(i))
           TEMY = 0.1d0*(U(I-M2) + U(I+M2) - 2.0d0*u(i))
           DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF - ftmp
 20     CONTINUE
 30   CONTINUE
C
C
c... bottom 
      k = 0
      ioff = M2*k
      do j = 1, M
         I = IOFF + J + 1
         etmp = dexp(-20.0d0/u(I))
         ftmp = senpar(1)*damk*(2.0d0 - u(I))*etmp
         TEMX = 0.1d0*(U(I-1)+U(I+1)  - 2.0d0*u(i))
         TEMY = 0.1d0*(2.0d0*U(I+M2) - 2.0d0*u(i))
         DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF - ftmp
      end do
C
C... left
      do k = 1,M
         ioff = M2*k
         j = 0
         I = IOFF + J + 1
         etmp = dexp(-20.0d0/u(I))
         ftmp = senpar(1)*damk*(2.0d0 - u(I))*etmp
         TEMX = 0.1d0*(2.0d0*U(I+1)  - 2.0d0*u(i))
         TEMY = 0.1d0*(U(I-M2) + U(I+M2) - 2.0d0*u(i))
         DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF - ftmp
      end do
c
c.. low-left corner
      k = 0
      j = 0
      I = 1
      etmp = dexp(-20.0d0/u(I))
      ftmp = senpar(1)*damk*(2.0d0 - u(I))*etmp
      TEMX = 0.1d0*(2.0d0*U(I+1)  - 2.0d0*u(i))
      TEMY = 0.1d0*(2.0d0*U(I+M2) - 2.0d0*u(i))
      DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF - ftmp
      
      RETURN
C------------  End of Subroutine RESH  ---------------------------------
      END

