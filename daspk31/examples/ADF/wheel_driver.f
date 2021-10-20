       program wheel
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, h_res, j_res, G_RES, G_G_RES
      DIMENSION Y(51), YPRIME(51), RWORK(3000),IWORK(1000)
      DIMENSION INFO(30),RPAR(1),IPAR(1), senpar(2), YTRUE(17)
      logical consis
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-6
      ATOL = 1.0d-6
      LRW = 3000
      LIW = 1000
      T = 0.0d0
      TOUT = 10.0d0
      neq = 17
*     
*     Initialize y, yprime, senpar
*     
C     *           MU    = 0.120000d0 ,
C     *           XL    = 0.19d0,
      senpar(1) = 0.120000d0
      senpar(2) = 0.19d0
      call init(neq,t,y,yprime,consis)
c
c$$$c...  free index-1 variables
c$$$      y(12) = 8.0d-5
c$$$c$$$c      y(13) = 0.15213d0
c$$$      y(14) = 9.0d-5
c$$$      y(15) = 0.2d0
c
c...  free index-2 variables
      do i = 16,17
         y(i) = 1.0d0
         yprime(i) = 0.0d0
      end do
c
c... free YPRIME
      do i = 1, 11
         yprime(i) = 0.0d0
      end do
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
C     info(5) = 0 only works when consist = .true.
*     
      info(5) = 3
C     
C     Here set INFO(14) = 1 to get the computed initial values.
      INFO(14) = 1
*     
*     Set #12 #13 as index-2 algebraic variables, does not do the error test on it.
*     
      info(16) = 1
      do i = 1, 11
         iwork(40+i) = 1
         iwork(40+neq+i) = 0
      end do
c
c...  fixed differential variables
      do i = 1, 9
         iwork(40+i) = 2 
      end do
      iwork(40+11) = 2
c
c...  index-1 variables, fixed
      do i = 1, 4
         iwork(40+11+i) = -1
      end do
c
c...  index-2 variables, not fixed
      iwork(40+16) = -2
      iwork(40+17) = -2
c
c...  index-2 constraints
      iwork(40+neq+12) = 1
      iwork(40+neq+13) = 1
c
c...  index-1 constraints
      do i = 14, 17
         iwork(40+neq+i) = 1
      end do

      iwork(38) = 0
      Ny = neq
c
c...  initialization options
      INFO(11) = 5              ! 4 or 5 
      if (consis) info(11) = 0
      inso = 18
      info(inso+1) = 0
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
            print *, (y(i),i=12,17)
c            stop
         end if
         if (idid.gt.0) GOTO 20
      ENDIF 

      print *
      call solut(ny,t,ytrue)
      do i = 1, 13
         print *, y(i) - ytrue(i)
      end do
      print *, (y(i*ny+16), y(i*ny+17), i = 1, info(19))
      STOP
      END

