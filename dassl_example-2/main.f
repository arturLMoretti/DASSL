C234567
      SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      double precision y(90),yprime(90),delta(90),rpar(1000)   
      double precision vel(93),pos(93),am(1000,1000),g(20),gp(22,22)
      double precision f(100),rl1(100),rl2(100),gi(20), rl(100)
      integer ires,ipar(1000)
      n=13
      do ii = 1, n
      pos(ii)=y(ii)  
      vel(ii)=y(ii+n)
      end do
      do i = 1,2
       rl(i) = y(2*n+i) 
      end do
      call slider(n,vel,pos,f,am,g,gi,gp,t,rl)       
      do i=1,n
       rl1(i)=0.0
      do j=1,2
        rl1(i)=rl1(i)+gp(j,i)*y(2*n+j)
      end do
      end do
      do ii = 1, n
      delta(ii)=yprime(ii)-y(n+ii)
      end do
      do i=1,n
      delta(i+n)=0.0
      do j=1,n
      delta(i+n)=delta(i+n)+am(i,j)*yprime(n+j)
      end do
      delta(i+n)=delta(i+n)-f(i)+rl1(i)
      end do
      do i=1,2
       delta(2*n+i)=0.0
      do j=1,n
        delta(2*n+i)=delta(2*n+i)+gp(i,j)*yprime(j+n)
      end do
        delta(2*n+i)=delta(2*n+i)+gi(i)
      end do
      return
      end
C
      program Slider_benchmark
      external res
      external jac
      INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(100000), LIW, IPAR(1000)
      DOUBLE PRECISION x(93),pos(93),vel(93),rl(2),
     * T, Y(90), YPRIME(90), TOUT, RTOL, ATOL, RWORK(100000),
     * RPAR(1000),am(1000,1000),f(100),gp(22,22),g(22),gi(20)
      t=0.0
c     rwork(1) = 1.00
      tout=1D0
      info(1)=0
      info(2)=0
      info(3)=1
      info(4)=0
      info(5)=0
      info(6)=0
      info(7)=0
      info(8)=0
      info(9)=0
      info(10)=0
      info(11)=0
      info(12)=0
      info(13)=0
      info(14)=0
      info(15)=0
      lrw=100000
      liw=100000
      n = 13
      do ii = 1,n
       y(ii) = 0.0
       y(ii+n) = 0.0
      end do
      do ii = 1,2
      y(2*n + ii) = 0.0
      rl(ii) = 0.0
      end do
      y(3)=0.45D0
c     y(8) = 0.0001
      do ii = 1, n
      pos(ii)=y(ii)  
      vel(ii)=y(ii+n)
      end do
      call slider(n,vel,pos,f,am,g,gi,gp,t,rl)
      call gauss(n,am,f,x)
      write(*,*) (x(i),i=1,13)    
c     pause
c     do ii = 1, n
c      yprime(ii+n) = x(ii) 
c     end do
      neq=2*n+2
      rtol = 0.00000001D0 
      atol = 0.00000001D0
      open(unit=12,file='sl.dat')
      idid = 0
  30  call DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
     + IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      write(*,*) idid
      write(12,24) t,(y(i),i=1,n)    
      write(*,24) t,(y(i),i=1,n)    
      if((idid .lt. 0) .or. (idid .gt. 3)) then 
       write(*,*) idid
       stop
       endif
      if (t .lt. tout ) go to 30
  24  format(E14.6,14G13.5)
      stop
      end

      subroutine jac
      stop 
      end
