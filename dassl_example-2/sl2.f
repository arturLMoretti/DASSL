c234567
      subroutine slider(n,v,p,f,am,g,gi,gp,t,rlam)
      implicit none
      integer n
c
      integer i,j
      INTEGER    NP, NV, NL, NG, NU, LDG, IPAR(5), IFAIL
      LOGICAL    LFLAG(10)
      double precision T, P(N), V(N), U(100), RL(100), AM(1000,1000),
     &      GP(22,22), F(N), PDOT(100), UDOT(100), G(22),
     &      GI(20), FL(100,100), RPAR(9),EN, RLAM(100)
C
      np=n
      nv=n
      nl=2
      ng=2
      ldg=n+1
       do i=1,8
       lflag(i)=.true.
       end do
       lflag(10)=.true.
c
      ipar(1)=4
      rpar(1)=0.2
      rpar(2)=0.3
      rpar(3)=9.81

c     rpar(3)=0.0

      do i=4,ipar(1)+3
       rpar(i)=0.075
      end do
c      write(*,*) 'ipar',(ipar(i),i=1,5)
C      write(*,*) 'rpar',(rpar(i),i=1,ipar(1)+3)
c
      ldg = 22
      
      CALL FMBSE (NP, NV, NL, NG, NU, LDG, T, P, V, U,
     &         RLAM, AM, GP, F, PDOT, UDOT, G, GI,
     &         FL, LFLAG, RPAR, IPAR, IFAIL)

      return
      end

