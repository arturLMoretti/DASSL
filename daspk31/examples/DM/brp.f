       program brp 
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), integer(i-n)
       EXTERNAL RES, QRES, JAC
       DIMENSION Y(110),RWORK(9400),IWORK(9400),YPRIME(110) 
       DIMENSION INFO(30),RPAR(8),IPAR(8), QSEN(9), 
     *      atol(110), rtol(110)
*
* 
* Initialize the INFO vector.
*
        DO 10 I = 1,30
          INFO(I) = 0
  10    CONTINUE
*
* Set RTOL values. 
*
        info(2) = 1
        do i = 1, 110
           rtol(i) = 1.0d-7
           atol(i) = 1.0d-7
        end do
*
* Scale ATOL values for each sensitivity segment using the rule:
*
* RTOL(i)/O(pi) = ATOL(i), for i = 1,...,Ny.
*
* Set scaled ATOL values for p1.
       ny = 10
       info(19) = 8
       neq = ny*(info(19) + 1)
*
*
* Set LRW, LIW to large values (see DASPK documentation).
*
       LRW = 9400
       LIW = 9400
*
* For this problem the total number of equations, NEQ, is 90.
*
       T = 0.0
       TOUT = 10.0
*
* Initialize vectors.
*
       call init(y,yprime,ipar,rpar,neq,senpar)
*
       info(5) = 1
*
*
* Make the  IC's consistent and tell DASSLSO that they are.
*
       call consist(y,yprime,ipar,rpar)
       INFO(11) = 0
       info(16) = 1
c
c...  set up the IWORK for inconsistent initial conditions
       do i = 1, 6
          iwork(40+i) = 1
       end do
       do i = 7, 10
          iwork(40+i) = -1
       end do
       info(3) = 1
       info(24) = 1             ! derived function
*
* Do that DASSLSO magic....
        icnt = 0
  20    IF(T .LT. TOUT) THEN 
*
* Call DASPK until integration is complete.
*
        CALL DDASPK(RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL,IDID,RWORK,
     c  LRW,IWORK,LIW,RPAR,IPAR,JAC,PSOL, SENPAR, G_RES)
        if (idid .lt. 0) stop
*
!         print *, T, Y(1)
        WRITE(11,*) T, (Y(i), i=1,10)

        icnt = icnt + 1
*
* Take every 10th data point.
* 
        if( icnt.eq.5 )then
!           write(*,*)'time = ',T,'  CSS = ',rwork(7)
          icnt = 0
        endif
*
        GOTO 20
        ENDIF 
        CALL DSENSD(
     *       QRES, NEQ, T, Y, YPRIME, QSEN, INFO, 
     *       RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)

        do i = 1, 9
           write(12,*) qsen(i)
!            print *, qsen(i)
        end do
        
*
* Output IWORK and RWORK vector values
*
        WRITE(*,*)'SOLUTION FOR THE Batch Reactor'
        write(*,*)'PROBLEM REQUIRED:'
        PRINT*
        WRITE(*,*)'TOTAL # OF STEPS:',IWORK(11)
        PRINT*
        WRITE(*,*)'TOTAL # OF ERROR TEST FAILURES:',IWORK(14)
        PRINT*
        WRITE(*,*)'TOTAL # OF CONV. TEST FAILURES:',IWORK(15)
        PRINT*
        WRITE(*,*)'TOTAL # OF RES function evals:',IWORK(12)
*
        END
* ------------------------------------------------        
       subroutine qres(T, Y, YPRIME, Q, IRES, RPAR, IPAR, SENPAR)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION Y(*),YPrime(*),q(*),RPar(*),IPar(*), senpar(*)
       
       q(1) = 0.0d0
       do i = 1, 10
          q(1) = q(1) + y(i)*y(i)
       end do
!        print *, 'QRES returning ', q(1)
       
       return
       end

      
*--------------------------------------------------------------------------*
* Problem subroutine. 
*--------------------------------------------------------------------------*
       SUBROUTINE RES(T,Y,YP,CJ,DELTA,IRES,SENPAR,IP,RP) !! flipped senpar and rp
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION Y(*),YP(*),DELTA(*),RP(*),IP(*),SENPAR(*)
*
* Equation set-up. 
*
       delta(1) = yp(1) + senpar(ip(3))*y(2)*y(8)
       delta(2) = yp(2) + senpar(ip(1))*y(2)*y(6)
     *            - senpar(ip(2))*y(10)
     *            + senpar(ip(3))*y(2)*y(8)
       delta(3) = yp(3) - senpar(ip(3))*y(2)*y(8)
     *            - senpar(ip(4))*y(4)*y(6)
     *            + senpar(ip(5))*y(9)
       delta(4) = yp(4) + senpar(ip(4))*y(4)*y(6) 
     *            - senpar(ip(5))*y(9)
       delta(5) = yp(5) - senpar(ip(1))*y(2)*y(6) 
     *            + senpar(ip(2))*y(10)
       delta(6) = yp(6) + senpar(ip(1))*y(2)*y(6) 
     *            + senpar(ip(4))*y(4)*y(6)
     *            - senpar(ip(2))*y(10) - senpar(ip(5))*y(9)
       delta(7) = -1.31d-2 + y(6) + y(8) + y(9) + y(10) - y(7)
       delta(8) = senpar(ip(7))*y(1)  
     *            - y(8)*( senpar(ip(7)) + y(7) )
       delta(9) = senpar(ip(8))*y(3) 
     *            - y(9)*( senpar(ip(8)) + y(7) )
       delta(10) = senpar(ip(6))*y(5) 
     *             - y(10)*( senpar(ip(6)) + y(7) )

       RETURN
       END
*--------------------------------------------------------------------------*
       SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IP, senpar, ijac)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION Y(*), YPRIME(*), PD(10,*), RPAR(*), IP(*),senpar(*)
       integer ijac
*
* Analytic Jacobian.
*
       PD(1,1) =  CJ
       PD(1,2) =  senpar(ip(3))*y(8)
       PD(1,8) =  senpar(ip(3))*y(2)
*
       PD(2,2) =  senpar(1)*y(6) + senpar(3)*y(8) + CJ
       PD(2,6) =  senpar(1)*y(2)
       PD(2,8) =  senpar(3)*y(2)
       PD(2,10) = -senpar(2)
*
       PD(3,2) =  -senpar(3)*y(8)
       PD(3,3) =  CJ
       PD(3,4) =  -senpar(4)*y(6)
       PD(3,6) =  -senpar(4)*y(4)
       PD(3,8) =  -senpar(3)*y(2)
       PD(3,9) =  senpar(5)
*     
       PD(4,4) =  senpar(4)*y(6) + CJ
       PD(4,6) =  senpar(4)*y(4)
       PD(4,9) =  -senpar(5)
*
       PD(5,2) =  -senpar(1)*y(6)
       PD(5,5) =  CJ
       PD(5,6) =  -senpar(1)*y(2)
       PD(5,10) = senpar(2)
*
       PD(6,2) =  senpar(1)*y(6)
       PD(6,4) =  senpar(4)*y(6)
       PD(6,6) =  senpar(1)*y(2) + senpar(4)*y(4) + CJ
       PD(6,9) =  -senpar(5)
       PD(6,10) = -senpar(2)
*
       PD(7,6) =  1.0d0
       PD(7,7) =  -1.0d0
       PD(7,8) =  1.0d0
       PD(7,9) =  1.0d0
       PD(7,10) = 1.0d0
*
       temp =  senpar(7) + y(7)
c       PD(8,1) =  rpar(7)/temp
c       PD(8,7) =  rpar(7)*y(1)/( temp*temp )
c       PD(8,8) =  -1.0d0
       PD(8,1) = senpar(7)
       PD(8,7) = -y(8)
       PD(8,8) = -temp
*
       temp = senpar(8) + y(7)
c       PD(9,3) =  rpar(8)/temp
c       PD(9,7) =  rpar(8)*y(3)/( temp*temp )
c       PD(9,9) =  -1.0d0
       PD(9,3) = senpar(8)
       PD(9,7) = -y(9)
       PD(9,9) = -temp
*
       temp = senpar(ip(6)) + y(7)
c       PD(10,5) =  rpar(ip(6))/temp
c       PD(10,7) =  rpar(ip(6))*y(5)/( temp*temp )
c       PD(10,10) =  -1.0d0
       PD(10,5) = senpar(6)
       PD(10,7) = -y(10)
       PD(10,10) = -temp
*
       if (ijac .eq. 1) then
          do i = 1, 10
             do j = 1, 6
                if (i .ne. j) then
                   PD(i,j) = 0.0d0
                end if
             end do
          end do
       end if
                
       return
       end
*---------------------------------------------------------------------------*
       SUBROUTINE consist(Y,YP,IP,RP)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION Y(*), YP(*), RP(*),IP(*)
*
* Consistency acheived by setting YPRIME to RHS, initially.
*
       t1 = rp(ip(3))*y(2)*y(8)
       t2 = rp(ip(1))*y(2)*y(6)
       t3 = rp(ip(4))*y(4)*y(6)
       t4 = rp(ip(5))*y(9)
       t5 = rp(ip(2))*y(10)
*
       yp(1) = -t1
       yp(2) = -t2 + t5 - t1
       yp(3) = t1 + t3 - t4
       yp(4) = -t3 + t4
       yp(5) = t2 - t5
       yp(6) = -t2 - t3 + t5 + t4
*
       return
       end
*--------------------------------------------------------------------------* 
       SUBROUTINE init(y,yp,ip,rp,neq,senpar)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z), INTEGER(I-N)
       DIMENSION y(*), rp(*), ip(*), yp(*), senpar(*)
*
       do i = 1,8
        rp(i) = 0.0d0
        ip(i) = 0
       enddo
*
* Place parameter values in RPAR.
*
        rp(1) = 21.893d0
        rp(2) = 2.14d+9
        rp(3) = 32.318d0
        rp(4) = 21.893d0
        rp(5) = 1.07d+9
        rp(6) = 7.65d-18
        rp(7) = 4.03d-11
        rp(8) = 5.32d-18
        
        do i = 1,8
          senpar(i) = rp(i)
        end do
*
* Set IPAR values to point to corresponding parameter vlues in RPAR.
*
        ip(1) = 1  
        ip(2) = 2 
        ip(3) = 3 
        ip(4) = 4
        ip(5) = 5
        ip(6) = 6
        ip(7) = 7
        ip(8) = 8
*
* Set Y, YPRIME to initial values.
*
       do i = 1,neq
          yp(i) = 0.0d0
          y(i) = 0.0d0
       enddo

       y(1) = 1.5776d0
       y(2) = 8.32d0
       y(3) = 0.0d0
       y(4) = 0.0d0
       y(5) = 0.0d0
       y(6) = 1.31d-2
       t2 = rp(ip(7))*(rp(ip(7)) + 4.0d0*y(1))
       y(7) = 0.5d0*( dsqrt(t2) - rp(ip(7)) )
       y(8) = y(7)
       y(9) = 0.0d0
       y(10) = 0.0d0
       do i = 1, 10
          y(10*i+i) = 1.0d0
       end do
*
       return
       end
*--------------------------------------------------------------------------*



