c234567
C
      subroutine slider(n,v,p,f,am,g,gi,gp,t)
      integer n
c 
c     implicit none
      integer i,j
      INTEGER    NP, NV, NL, NG, NU, LDG, IPAR(5), IFAIL
      LOGICAL    LFLAG(10)
      double precision T, P(N), V(N), U(100), RL(100), AM(100,100),
     &           GP(22,22), F(N), PDOT(100), UDOT(100), G(100),
     &           GI(100), FL(100,100), RPAR(100)
C
      np=n
      nv=n
      nl=2
      ng=2
      ldg=100
       do i=1,10
       lflag(i)=.true.
       end do
       lflag(8)=.false.
c
      rpar(1)=0.2
      rpar(2)=0.3 
      rpar(3)=9.81
c
      call FMBSE (NP, NV, NL, NG, NU, LDG, T, P, V, U,
     &                 RLAM, AM, GP, F, PDOT, UDOT, G, GI,
     &                 FL, LFLAG, RPAR, IPAR, IFAIL)

      return
      end 
