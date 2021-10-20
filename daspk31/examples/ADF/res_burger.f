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
      dx = rpar(3)
      dy = dx
      dx2 = dx*dx
      dy2 = dx2
      
      M2 = M + 2
      eps = senpar(1)
C
C Load U into DELTA, in order to set boundary values.
C
C Loop over interior points, and load residual values.
      DO K = 1,M
        IOFF = M2*K
        DO J = 1,M
           I = IOFF + J + 1
           u2x = (u(i-1)+u(i+1) - 2.0d0*u(i))/dx2
           u2y = (U(I-M2) + U(I+M2) - 2.0d0*u(i))/dy2
           u1x = (u(i+1)**2 - u(i-1)**2)/(2.d0*dx)
           u1y = (u(i+M2)**2 - u(i-M2)**2)/(2.d0*dy)            
           delta(i) = uprime(i) - eps*(u2x + u2y)
           delta(i) = delta(i) + 0.5d0*(u1x + u1y)
        enddo
      enddo

c
c...  boundary
      k = 0                     ! bottom boundary
      do j = 1, m
         i = j+1
         ubot = 1.d0/(1.d0 + dexp((j*dx+(k-1)*dy-t)/(2.0d0*senpar(1))))
         u2x = (u(i-1)+u(i+1) - 2.0d0*u(i))/dx2
         u2y = (u(i+m2) + ubot - 2.0d0*u(i))/dy2
         u1x = (u(i+1)**2 - u(i-1)**2)/(2.d0*dx)
         u1y = (u(i+M2)**2 - ubot**2)/(2.d0*dy)
         delta(i) = uprime(i) - eps*(u2x + u2y)
         delta(i) = delta(i) + 0.5d0*(u1x + u1y)
      end do
      k = m+1
      do j = 1, m
         i = k*m2 + j+1
         utmp = 1.d0/(1.d0 + dexp((j*dx+(k+1)*dy-t)/(2.0d0*senpar(1))))
         u2x = (u(i-1)+u(i+1) - 2.0d0*u(i))/dx2
         u2y = (u(i-m2) + utmp - 2.0d0*u(i))/dy2
         u1x = (u(i+1)**2 - u(i-1)**2)/(2.d0*dx)
         u1y = (-u(i-M2)**2 + utmp**2)/(2.d0*dy)
         delta(i) = uprime(i) - eps*(u2x + u2y)
         delta(i) = delta(i) + 0.5d0*(u1x + u1y)
      end do
      do k = 1, m
         j = 0                  ! left
         i = M2*k + j +1
         utmp = 1.d0/(1.d0 + dexp(((j-1)*dx+k*dy-t)/(2.0d0*senpar(1))))
         u2x = (utmp+u(i+1) - 2.0d0*u(i))/dx2
         u2y = (u(i-m2) + u(i+m2) - 2.0d0*u(i))/dy2
         u1x = (u(i+1)**2 - utmp**2)/(2.d0*dx)
         u1y = (u(i+M2)**2 - u(i-m2)**2)/(2.d0*dy)
         delta(i) = uprime(i) - eps*(u2x + u2y)
         delta(i) = delta(i) + 0.5d0*(u1x + u1y)
      end do
      do k = 1, m
         j = m+1                ! right
         i = M2*k + j +1
         utmp = 1.d0/(1.d0 + dexp(((j-1)*dx+k*dy-t)/(2.0d0*senpar(1))))
         u2x = (utmp+u(i-1) - 2.0d0*u(i))/dx2
         u2y = (u(i-m2) + u(i+m2) - 2.0d0*u(i))/dy2
         u1x = (-u(i-1)**2 + utmp**2)/(2.d0*dx)
         u1y = (u(i+M2)**2 - u(i-m2)**2)/(2.d0*dy)
         delta(i) = uprime(i) - eps*(u2x + u2y)
         delta(i) = delta(i) + 0.5d0*(u1x + u1y)
      end do

c      low-left corner
      k = 0
      j = 0
      i = 1
      uleft = 1.d0/(1.d0 + dexp(((j-1)*dx+k*dy-t)/(2.0d0*senpar(1))))
      ubot  = 1.d0/(1.d0 + dexp((j*dx+(k-1)*dy-t)/(2.0d0*senpar(1))))
      u2x = (uleft+u(i+1) - 2.0d0*u(i))/dx2
      u2y = (ubot + u(i+m2) - 2.0d0*u(i))/dy2
      u1x = (-uleft**2 + u(i+1)**2)/(2.d0*dx)
      u1y = (u(i+M2)**2 - ubot**2)/(2.d0*dy)
      delta(i) = uprime(i) - eps*(u2x + u2y)
      delta(i) = delta(i) + 0.5d0*(u1x + u1y)

c      low-right corner
      k = 0
      j = m+1
      i = m+2
      uright = 1.d0/(1.d0 + dexp(((j+1)*dx+k*dy-t)/(2.0d0*senpar(1))))
      ubot  = 1.d0/(1.d0 + dexp((j*dx+(k-1)*dy-t)/(2.0d0*senpar(1))))
      u2x = (u(i-1)+uright - 2.0d0*u(i))/dx2
      u2y = (ubot + u(i+m2) - 2.0d0*u(i))/dy2
      u1x = (uright**2 - u(i-1)**2)/(2.d0*dx)
      u1y = (u(i+M2)**2 - ubot**2)/(2.d0*dy)
      delta(i) = uprime(i) - eps*(u2x + u2y)
      delta(i) = delta(i) + 0.5d0*(u1x + u1y)

c      top-right corner
      k = m+1
      j = m+1
      i = m2*m2
      uright = 1.d0/(1.d0 + dexp(((j+1)*dx+k*dy-t)/(2.0d0*senpar(1))))
      utop  = 1.d0/(1.d0 + dexp((j*dx+(k+1)*dy-t)/(2.0d0*senpar(1))))
      u2x = (u(i-1)+uright - 2.0d0*u(i))/dx2
      u2y = (utop + u(i-m2) - 2.0d0*u(i))/dy2
      u1x = (uright**2 - u(i-1)**2)/(2.d0*dx)
      u1y = (-u(i-M2)**2 + utop**2)/(2.d0*dy)
      delta(i) = uprime(i) - eps*(u2x + u2y)
      delta(i) = delta(i) + 0.5d0*(u1x + u1y)

c      top-left corner
      k = m+1
      j = 0
      i = (m+1)*m2+1
      uleft = 1.d0/(1.d0 + dexp(((j-1)*dx+k*dy-t)/(2.0d0*senpar(1))))
      utop  = 1.d0/(1.d0 + dexp((j*dx+(k+1)*dy-t)/(2.0d0*senpar(1))))
      u2x = (u(i+1)+uleft - 2.0d0*u(i))/dx2
      u2y = (utop + u(i-m2) - 2.0d0*u(i))/dy2
      u1x = (-uleft**2 + u(i+1)**2)/(2.d0*dx)
      u1y = (-u(i-M2)**2 + utop**2)/(2.d0*dy)
      delta(i) = uprime(i) - eps*(u2x + u2y)
      delta(i) = delta(i) + 0.5d0*(u1x + u1y)

C
      RETURN
C------------  End of Subroutine RESH  ---------------------------------
      END
