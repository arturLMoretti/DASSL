       program fekete
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
      EXTERNAL RES, h_res, j_res, G_RES, G_G_RES, qres
      DIMENSION Y(160*9), YPRIME(160*9), RWORK(200000),IWORK(10000)
      DIMENSION INFO(30),RPAR(1),IPAR(1), senpar(8), YTRUE(160)
      dimension ysen(8*160), ypsen(160*8), tmp(8*8), qsen(10)

      real*4   t1, t2, t3, tt(2)
      logical consis
*     
*     Specify Error tolerances, length of iwork and rwork, NEQ and the
*     interval of integeration.
*     
      RTOL = 1.0d-4
      ATOL = 1.0d-4
      LRW = 200000
      LIW = 10000
      T = 0.0d0
      TOUT = 1000.0d0
      neq = 160
*     
*     Initialize y, yprime, senpar
*     
      senpar(1) = 1d0/13d0
      senpar(2) = 3d0/8d0
      senpar(3) = 1d0/29d0
      senpar(4) = 1d0/8d0
      senpar(5) = 1d0/7d0
      senpar(6) = -2d0/15d0
      senpar(7) = 1d0/17d0
      senpar(8) = -0.3d0
      
      call init(neq,t,y,yprime,consis,senpar)
      do i = 1, 64
         tmp(i) = 0.0d0
      end do
      do i = 1, 8
         tmp(i+(i-1)*8) = 1.0d0
      end do
c
c...  calculate the initial condition for the sensitivity variables
c     by calling ADIFOR generated routine
      call l_init(8,neq,t,y,ysen,8,yprime,ypsen,8,consis,senpar,
     *     tmp, 8)
      do i = 1, 8
         do j = 1, neq
            y(i*neq + j) = ysen(i+(j-1)*8)
            yprime(i*neq + j) = ypsen(i+(j-1)*8)
         end do
      end do
c
c
c...  free index-2 variables
      do i = 121, 160
         y(i) = 1.0d0
      end do
c
c... free YPRIME
      do i = 1, 120
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
      INFO(3)  = 0
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
      do i = 1, 120
         iwork(40+i) = 1
         iwork(40+neq+i) = 0
      end do
c
c...  index-2 variables and constraints
      do i = 121, 160
         iwork(40+i) = -2
         iwork(40+i+neq) = 1
      end do

      iwork(38) = 0
      Ny = neq
c
c...  initialization options
      INFO(11) = 5              ! 4 or 5 
      if (consis) info(11) = 0
      inso = 18
      info(inso+1) = 8
      neq = neq*(info(inso+1) + 1)
      info(inso+2) = 0
      info(inso+3) = 0          ! default perturbed factor
      info(inso+4) = 0          ! no sensitivity parameter in the RES
      info(inso+5) = 1          ! 0, full control 1. partial control
      info(inso+6) = 1          ! one derived function
      
      t1 = dtime(tt)
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
            print *, (y(i),i=121,125)
c            stop
         end if
         if (idid.gt.0) GOTO 20
      ENDIF 

      print *
      call solut(ny,t,ytrue)
      do i = 1, 13
         print *, y(i) - ytrue(i)
      end do

      print *, ' maximum value of log(V(x)) = ', value(20,y)
      CALL DSENSD(
     *     QRES, NEQ, T, Y, YPRIME, QSEN, INFO, 
     *     RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)
      print *, (qsen(i),i=1,9)
      t2 = dtime(tt)
      write(6,*) '   Time it took:', tt(1)+tt(2)
      STOP
      END

      double precision function value(n,v)
      implicit none
      integer n, i, j, k
      double precision v(3,n), prod, dij, dsqrt
      
      prod = 1.0d0
      do i = 1, n
         do j = i+1, n
            dij = 0.0d0
            do k = 1, 3
               dij = dij + (v(k,i)-v(k,j))**2
            end do
            prod = prod*dsqrt(dij)
         end do
      end do
      value = dlog(prod)
      return
      end
      SUBROUTINE QRES (T,Y,YPRIME,QSEN,ires,RPAR,IPAR,senpar)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(*), YPRIME(*), RPAR(*), IPAR(*),senpar(*), Qsen(*)
C
C Set problem constants using IPAR and RPAR.
C
      qsen(1) = value(20, Y)

      RETURN
      END
      

