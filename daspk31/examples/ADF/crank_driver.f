       program crank
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, h_res, j_res, G_RES, G_G_RES
      DIMENSION Y(24*3), YPRIME(24*3), RWORK(4000),IWORK(1000)
      DIMENSION INFO(30),RPAR(10),IPAR(2), senpar(2), YTRUE(24)
      logical consis
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-6
      ATOL = 1.0d-6
      LRW = 4000
      LIW = 1000
      T = 0.0d0
      TOUT = 0.1d0
      neq = 24
*     
*     Initialize y, yprime, senpar
*     
c     *             L1 = 0.15D0,     L2 = 0.30D0,
      senpar(1) = 0.15d0
      senpar(2) = 0.30d0
      call init(neq,t,y,yprime,consis)
c
c...  free index-1 variable, acceleration
      do i = 15, 21
         y(i) = 1.0d0
      end do
c
c..   free index-2 variables
      do i = 22, 24
         y(i) = 1.0d0
         yprime(i) = 0.0d0
      end do
c
c... free YPRIME
      do i = 1, 24
         yprime(i) = 0.0d0
      end do
c
c...  free V
c      do i = 1, 7
c         y(i+7) = 0.0d0
c      end do
      consis = .false.
*     
*     Initialize the INFO vector to 0.
*     
      DO 10 I = 1,30
         INFO(I) = 0
 10   CONTINUE
*     
*     Get solution at intermediate steps.
*     
      INFO(3)  = 1
*     
*     Adifor evaluation of the Jacobian. (1 = no,0 = yes)
*     
      info(5) = 3
C     
C     Here set INFO(14) = 1 to get the computed initial values.
      INFO(14) = 1
*     
*     Set #12 #13 as index-2 algebraic variables, does not do the error test on it.
*     
      info(16) = 1
      do i = 1, 14
         iwork(40+i) = 1        ! differential variables
         iwork(40+neq+i) = 0    ! differential equations
      end do
      do i = 1, 7
         iwork(40+i) = 2        ! fixed variables
      end do
      do i = 15, 24
         iwork(40+i) = -2       ! unknown algebraic values
      end do

      do i = 1, 3
         iwork(40+neq+21+i) = 1 ! constraints need differentiating
      end do
      do i = 15, 21
         iwork(40+neq+i) = 0    ! constraints not need differentiating
      end do

      Ny = neq
*     
*     
*     Make the IC's consistent and relay this info. to DASSLSO.
*     
      INFO(11) = 5
      if (consis) info(11) = 0
      inso = 18
      info(inso+1) = 2
      neq = neq*(info(inso+1) + 1)
      info(inso+2) = 3
      info(inso+3) = 0          ! default perturbed factor
      info(inso+4) = 2          ! no sensitivity parameter in the RES
      info(inso+5) = 1          ! 0, full control 1. partial control
      info(inso+6) = 0          ! no derived information
      
      do i = ny+1, neq
         y(i) = 0.0d0
         yprime(i) = 0.0d0
      end do

 20   IF(T .LT. TOUT) THEN     
         if (info(19) .gt. 0) then
            CALL DDASPK(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, H_res, 
     *           K_RES, G_G_RES)
         else 
            CALL DDASPK(
     *           RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     *           LRW,IWORK,LIW,RPAR,IPAR,j_res, psol,senpar, H_res, 
     *           K_RES, G_RES)
         end if

         print *, t, iwork(10)
         if (info(11).ne.0)  then
            info(11) = 0
            print *, 'Y(7) .. Y(14):'
            print *, (y(i),i=7,14)
            print *, 'Y(15) .. Y(24):'
            print *, (y(i),i=15,24)
c            stop
         end if
         if (idid.gt.0) GOTO 20
      ENDIF 

      print *
      call solut(ny,t,ytrue)
      do i = 1, 24
         print *, y(i) - ytrue(i)
      end do
      STOP
      END

