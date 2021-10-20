       program testad2
*  
*      *** Demo program for DASPK ***
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, jac
      DIMENSION Y(10),RWORK(3000),IWORK(1000),YPRIME(10), infob(10)
      DIMENSION INFO(30),RPAR(1),IPAR(20), senpar(1), ieopt(10)
      dimension qsen(30), rtol(6), atol(6)
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      LRW = 3000
      LIW = 1000
      T = 0.0d0
      TOUT = 1.0d0
      neq = 3
*     
*     Initialize y, yprime, senpar
*     
      call init(y,yprime,senpar)
*     
*     Initialize the INFO vector to 0.
*     
      DO 10 I = 1,30
         INFO(I) = 0
 10   CONTINUE
*     
*     Get solution at intermediate steps.
*     
      info(2)  = 1
      do i = 1, neq*2
         RTOL(i) = 1.0d-6
         ATOL(i) = 1.0d-6
      end do
      rtol(6) = 1.0d0
      INFO(3)  = 1
*     
*     Finite difference approx. (FDA) of Jacobian. (1 = no,0 = yes)
*     
      info(5) = 1
c
c...  exclude the algebraic variables from the error test
      info(16) = 1
      do i = 1, 2
         iwork(40+i) = 1        ! differential variables
         iwork(40+neq+i) = 0      ! differential equations
      end do

      iwork(43) = -2            !  index-2 variables
      iwork(46) = 1             !  index-2 constraint
      iwork(38) = 0
*     
      info(19) = 1
      info(20) = 2
      info(23) = 1
      info(14) = 1
      INFO(11) = 5
      neq = (info(19)+1)*neq
 15   continue
      if (t .lt. tout) then
         CALL DDASPK(
     *        RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *        LRW,IWORK,LIW,RPAR,IPAR,jac, psol,senpar, H_res, 
     *        K_RES, G_G_RES)
         if (idid .lt. 0) then
            print *, ' idid =', idid
            stop
         end if
         if (info(11) .gt. 0) then
            info(11) = 0
            goto 15
         end if
         if (idid .eq. 1) goto 15
      end if
      print *, (y(i), i = 1, 6)
      stop
    
      END
*------------------------------------------------------------------------*
      SUBROUTINE RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR, senpar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), DELTA(*), RPAR(*), IPAR(*),senpar(*)     
      if (ires .eq. 0) then
         DELTA(1) = YPRIME(1) - y(1) - 2.d0*y(3)
         DELTA(2) = YPRIME(2) - y(2) - y(3)
         DELTA(3) = cj*(y(1)*y(1) + y(2)*y(2) - 2.0d0)
      else if (ires .eq. 1) then
         DELTA(4) = YPRIME(4) - y(4) - 2.d0*y(6)
         DELTA(5) = YPRIME(5) - y(5) - y(6)
         DELTA(6) = CJ*(y(1)*y(4) + y(2)*y(5))
      else if (ires .eq. 2) then
         DELTA(3) = y(1)*YPRIME(1) + y(2)*yprime(2)
         delta(6) = y(1)*yprime(4) + yprime(1)*y(4) + y(2)*yprime(5) 
     *        + yprime(2)*y(5)
      end if
      RETURN
      END 
*------------------------------------------------------------------------*

      subroutine jac(t, y, yp, pd, cj, rpar, ipar, senpar, ijac)
      implicit none
      double precision t, y(*), yp(*),pd(3,3),cj,rpar(*),senpar(*)
      integer ipar(*), ijac, i, j

      do i = 1, 3
         do j = 1,3
            pd(i,j) = 0.0d0
         end do
      end do
      
      if (ijac .eq. 0 .or. ijac.eq.1) then
         pd(1,1) = cj - 1.0d0
         pd(1,3) = -2.0d0
         pd(2,2) = cj - 1.0d0
         pd(2,3) = -1.d0
         pd(3,1) = 2.d0*cj*y(1)
         pd(3,2) = 2.d0*cj*y(2)
      else if (ijac .eq. 2) then
         pd(1,1) = cj 
         pd(1,3) = -2.0d0
         pd(2,2) = cj
         pd(2,3) = -1.d0
         pd(3,1) = cj*y(1)
         pd(3,2) = cj*y(2)
      end if

      return
      end

      SUBROUTINE INIT(Y,YPRIME,senpar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
      DIMENSION Y(*), YPRIME(*), senPAR(*)
*     
*     Set y, yprime, IPAR, RPAR to 0, initially.
*     
      do i = 1,2
         y(i) = 1.0d0
         yprime(i) = 0.0d0
      enddo
      y(3) = -1.0d0
      y(4) = 1.0d0
      y(5) = -1.0d0
      y(6) = 0.0d0
      do i = 4, 6
         yprime(i) = 0.0d0
      end do
*     
*     
*     Parameter values - placed somewhere in senPAR.
*     
      senpar(1) = 0.0d0
      return
      end
c
