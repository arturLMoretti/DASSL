      SUBROUTINE RES (T, U, UPRIME, CJ, DELTA,IRES,RPAR,IPAR,senpar)
C
C This is the user-supplied RES subroutine for this example.
C It computes the residuals for the 2-D discretized heat equation,
C with zero boundary values.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*), UPRIME(*), DELTA(*), RPAR(*), IPAR(*),senpar(*)
C
C Set problem constants using IPAR and RPAR.
      NEQ = IPAR(33)
      M = IPAR(34)
      COEFF = RPAR(4)
      M2 = M + 2
C
C Load U into DELTA, in order to set boundary values.
      if (ires .eq. 0) then
         DO 10 I = 1,NEQ
 10      DELTA(I) = UPRIME(I)
C     
C     Loop over interior points, and load residual values.
         DO 30 K = 1,M
            IOFF = M2*K
            DO 20 J = 1,M
               I = IOFF + J + 1
               TEMX = senpar(1)*(U(I-1)  + U(I+1)  - 2.0d0*u(i))
               TEMY = senpar(2)*(U(I-M2) + U(I+M2) - 2.0d0*u(i))
               DELTA(I) = UPRIME(I) - (TEMX + TEMY)*COEFF
 20         CONTINUE
 30      CONTINUE
      else if (ires .eq. 3) then
         sum = 0.0d0 
         do i = 1, neq
            sum = sum + u(i)
         enddo
         delta(neq+1) = uprime(neq+1) - sum
      end if
C      
      RETURN
C------------  End of Subroutine RESH  ---------------------------------
      END
