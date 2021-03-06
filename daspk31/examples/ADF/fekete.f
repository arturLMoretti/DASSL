c-----------------------------------------------------------------------
c
c     This file is part of the CWI Test Set for IVP solvers
c     http://www.cwi.nl/cwi/projects/IVPtestset/
c
c        Fekete problem (in stabilized index 2 formulation)
c        index 2 DAE of dimension 160
c
c     DISCLAIMER: see
c     http://www.cwi.nl/cwi/projects/IVPtestset/disclaimer.html
c
c     The most recent version of this source file can be found at
c     http://www.cwi.nl/cwi/projects/IVPtestset/src/problems/fekete.f
c
c     This is revision
c     $Id: fekete.f,v 1.1 1998/11/25 14:19:30 walter Exp $
c
c-----------------------------------------------------------------------
      integer function pidate()
      pidate = 19970704
      return
      end
c-----------------------------------------------------------------------
      subroutine prob(fullnm,problm,type,
     +                neqn,ndisc,t,
     +                numjac,mljac,mujac,
     +                nummas,mlmas,mumas,
     +                ind)
      character*(*) fullnm, problm, type
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(*)
      double precision t(0:*)
      logical numjac, nummas

      integer i, nart
      parameter (nart = 20)

      fullnm = 'Fekete problem'
      problm = 'fekete'
      type   = 'DAE'
      neqn   = 8*nart
      ndisc  = 0
      t(0)   = 0d0
      t(1)   = 1d3
      numjac = .false.
      mljac  = neqn
      mujac  = neqn
      mlmas  = 0
      mumas  = 0
      do 10 i=1,6*nart
         ind(i) = 1
   10 continue
      do 20 i=6*nart+1,8*nart
         ind(i) = 2
   20 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine init(neqn,t,y,yprime,consis,senpar)
      integer neqn
      double precision t,y(neqn),yprime(neqn), senpar(*)
      logical consis

      integer i,j,nart,ierr,ipar
      double precision pi,rpar,alpha,beta
      parameter(pi=3.141592653589793238462643383d0)

      nart=neqn/8

      do 10 i=1,3
         alpha=2*pi*dble(i)/dble(3)+pi*senpar(1)
         beta=pi*senpar(2)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
   10 continue
      do 20 i=4,10
         alpha=2*pi*dble(i-3)/dble(7)+pi*senpar(3)
         beta=pi*senpar(4)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
   20 continue
      do 30 i=11,16
         alpha=2*pi*dble(i-10)/dble(6)+pi*senpar(5)
         beta=pi*senpar(6)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
   30 continue
      do 40 i=17,20
         alpha=2*pi*dble(i-17)/dble(4)+pi*senpar(7)
         beta=pi*senpar(8)
         y(3*(i-1)+1)=cos(alpha)*cos(beta)
         y(3*(i-1)+2)=sin(alpha)*cos(beta)
         y(3*(i-1)+3)=sin(beta)
   40 continue

      do 50 i=3*nart+1,6*nart
         y(i)=0d0
   50 continue
      do 60 i=6*nart+1,8*nart
         y(i)=0d0
   60 continue
      call feval(neqn,0d0,y,y,yprime,ierr,rpar,ipar)
      do 80 i=1,nart
         do 70 j=1,3
           y(6*nart+i)=y(6*nart+i)+y(3*(i-1)+j)*yprime(3*nart+3*(i-1)+j)
   70    continue
         y(6*nart+i)=-y(6*nart+i)/2d0
   80 continue
      call feval(neqn,0d0,y,y,yprime,ierr,rpar,ipar)
      consis = .true.
      return
      end
c-----------------------------------------------------------------------
      subroutine res(t,y,yprime,cj,delta,ires,rpar,ipar,senpar)
      implicit none
      integer ires,ipar(*), i
      double precision t,cj,y(*),yprime(*),delta(*),rpar(*),senpar(*)
      integer neqn, nart
      nart = 20
      neqn = 8*20
      call feval(neqn,t,y,yprime,delta,ires,rpar,ipar)
      do i = 1, 120
         delta(i) = yprime(i) - delta(i)
      end do
c
c...  index-2 equations
      do i = 121, 160
         delta(i) = cj*delta(i)
      end do
      end
c-----------------------------------------------------------------------
      subroutine feval(neqn,t,y,yprime,dy,ierr,rpar,ipar)
      integer neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dy(neqn),rpar(*)

      integer i,j,k,nart,maxn
      parameter(maxn=150)
      double precision p(maxn,3),q(maxn,3),lam(maxn),mu(maxn),
     +                 pp(maxn,3),qp(maxn,3),phi(maxn),gpq(maxn),
     +                 f(maxn,maxn,3),rn,alpha

      alpha=.5d0
      nart=neqn/8

      do 20 i=1,nart
         do 10 k=1,3
            p(i,k)=y(3*(i-1)+k)
            q(i,k)=y(3*nart+3*(i-1)+k)
   10    continue
         lam(i)=y(6*nart+i)
         mu(i)=y(7*nart+i)
   20 continue
      do 70 i=1,nart
         do 60 j=1,nart
         if(i.eq.j)then
            do 30 k=1,3
               f(i,j,k)=0d0
   30       continue
         else
            rn=0d0
            do 40 k=1,3
               rn=rn+(p(i,k)-p(j,k))**2
   40       continue
            do 50 k=1,3
               f(i,j,k)=(p(i,k)-p(j,k))/rn
   50       continue
         endif
   60    continue
   70 continue
      do 100 i=1,nart
         do 90 k=1,3
            pp(i,k)=q(i,k)+2*mu(i)*p(i,k)
            qp(i,k)=-alpha*q(i,k)+2*lam(i)*p(i,k)
            do 80 j=1,nart
               qp(i,k)=qp(i,k)+f(i,j,k)
   80       continue
   90    continue
  100 continue
      do 120 i=1,nart
         phi(i)=-1d0
         gpq(i)=0d0
         do 110 k=1,3
            phi(i)=phi(i)+p(i,k)**2
            gpq(i)=gpq(i)+2*p(i,k)*q(i,k)
  110    continue
  120 continue
      do 140 i=1,nart
         do 130 k=1,3
            dy(3*(i-1)+k)=pp(i,k)
            dy(3*nart+3*(i-1)+k)=qp(i,k)
  130    continue
         dy(6*nart+i)=phi(i)
         dy(7*nart+i)=gpq(i)
  140 continue

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine meval(ldim,neqn,t,y,yprime,dfddy,ierr,rpar,ipar)
      integer ldim,neqn,ierr,ipar(*)
      double precision t,y(neqn),yprime(neqn),dfddy(ldim,neqn),rpar(*)

      integer i,nart

      nart=neqn/8
      do 10 i=1,6*nart
         dfddy(1,i)=1d0
   10 continue
      do 20 i=6*nart+1,8*nart
         dfddy(1,i)=0d0
   20 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine solut(neqn,t,y)
      integer neqn
      double precision t,y(neqn)
c
c computed at medusa
c problem           fekete
c solver            RADAU5
c rtol              0.100E-11
c atol              0.100E-11
c h0                0.100E-11
c ymax
c # scd                  0.55
c # steps                1249
c # steps accepted       1249
c # f-eval               9076
c # Jac-eval              589
c # LU-decomp             770
c CPU-time              22.87
c
      y(  1) =  -0.4070263380333202d+00
      y(  2) =   0.3463758772791802d+00
      y(  3) =   0.8451942450030429d+00
      y(  4) =   0.7752934752521549d-01
      y(  5) =  -0.2628662719972299d+00
      y(  6) =   0.9617122871829146d+00
      y(  7) =   0.7100577833343567d+00
      y(  8) =   0.1212948055586120d+00
      y(  9) =   0.6936177005172217d+00
      y( 10) =   0.2348267744557627d+00
      y( 11) =   0.7449277976923311d+00
      y( 12) =   0.6244509285956391d+00
      y( 13) =  -0.4341114738782885d+00
      y( 14) =   0.8785430442262876d+00
      y( 15) =   0.1992720444237660d+00
      y( 16) =  -0.9515059600312596d+00
      y( 17) =   0.2203508762787005d+00
      y( 18) =   0.2146669498274008d+00
      y( 19) =  -0.6385191643609878d+00
      y( 20) =  -0.4310833259390688d+00
      y( 21) =   0.6375425027722121d+00
      y( 22) =  -0.1464175087914336d+00
      y( 23) =  -0.9380871635228862d+00
      y( 24) =   0.3139337298744690d+00
      y( 25) =   0.5666974065069942d+00
      y( 26) =  -0.6739221885076542d+00
      y( 27) =   0.4740073135462156d+00
      y( 28) =   0.9843259538440293d+00
      y( 29) =  -0.1696995357819996d+00
      y( 30) =  -0.4800504290609090d-01
      y( 31) =   0.1464175087914331d+00
      y( 32) =   0.9380871635228875d+00
      y( 33) =  -0.3139337298744656d+00
      y( 34) =  -0.7092757549979014d+00
      y( 35) =   0.5264062637139616d+00
      y( 36) =  -0.4688542938854929d+00
      y( 37) =  -0.8665731819284478d+00
      y( 38) =  -0.4813878059756024d+00
      y( 39) =  -0.1315929352982178d+00
      y( 40) =  -0.2347897778700538d+00
      y( 41) =  -0.8594340408013130d+00
      y( 42) =  -0.4541441287957579d+00
      y( 43) =   0.5530976940074118d+00
      y( 44) =  -0.7674370265615124d+00
      y( 45) =  -0.3242273140037833d+00
      y( 46) =   0.7711050969896927d+00
      y( 47) =   0.6357041816577034d+00
      y( 48) =   0.3573685519777001d-01
      y( 49) =   0.7103951209379591d+00
      y( 50) =   0.2403570431280519d+00
      y( 51) =  -0.6614886725910596d+00
      y( 52) =  -0.3038208738735660d-01
      y( 53) =   0.4501923293640461d+00
      y( 54) =  -0.8924145871442046d+00
      y( 55) =  -0.5772996158107093d+00
      y( 56) =  -0.1766763414971813d+00
      y( 57) =  -0.7971892020969544d+00
      y( 58) =   0.2414481766969039d+00
      y( 59) =  -0.3416456818373135d+00
      y( 60) =  -0.9082846503446250d+00
      y( 61) =   0.2409619682166627d-15
      y( 62) =  -0.1139818460497816d-15
      y( 63) =   0.1627536276556335d-15
      y( 64) =   0.1745651819597609d-15
      y( 65) =  -0.1914278710633076d-15
      y( 66) =  -0.6639600671806291d-16
      y( 67) =   0.1708576733899083d-15
      y( 68) =  -0.2277602521390053d-15
      y( 69) =  -0.1350782790950654d-15
      y( 70) =   0.2411941341109454d-15
      y( 71) =  -0.1438238671800488d-15
      y( 72) =   0.8087033550666644d-16
      y( 73) =   0.1618239105233347d-15
      y( 74) =   0.1837556152070701d-16
      y( 75) =   0.2715177369929503d-15
      y( 76) =   0.7930078658689191d-16
      y( 77) =   0.7482020588342764d-16
      y( 78) =   0.2746974939098084d-15
      y( 79) =   0.8849338913035911d-16
      y( 80) =  -0.5940734725324115d-16
      y( 81) =   0.4845984056889910d-16
      y( 82) =  -0.3728835248155620d-16
      y( 83) =  -0.4600332954062859d-16
      y( 84) =  -0.1548568884846698d-15
      y( 85) =   0.2507541692375411d-16
      y( 86) =  -0.1560155223230823d-15
      y( 87) =  -0.2517946296860555d-15
      y( 88) =  -0.3739779361502470d-16
      y( 89) =  -0.1381663620885020d-15
      y( 90) =  -0.2784051540342329d-15
      y( 91) =   0.6624397102887671d-16
      y( 92) =   0.4226207488883120d-16
      y( 93) =   0.1571821772296610d-15
      y( 94) =  -0.4112243677286995d-16
      y( 95) =   0.1939960344265876d-15
      y( 96) =   0.2800184977692136d-15
      y( 97) =  -0.9189023375328813d-16
      y( 98) =   0.1392943179389155d-15
      y( 99) =   0.9556003995587458d-16
      y(100) =  -0.2234188557495892d-15
      y(101) =   0.1276804778190781d-15
      y(102) =  -0.1261196211463950d-15
      y(103) =  -0.1887754149742397d-15
      y(104) =  -0.2140788698695373d-16
      y(105) =  -0.2713591291421657d-15
      y(106) =   0.1107887633060814d-15
      y(107) =  -0.1318443715631340d-15
      y(108) =  -0.4521275683078691d-16
      y(109) =  -0.1277688851278605d-15
      y(110) =   0.4850914012115388d-16
      y(111) =  -0.1195891666741192d-15
      y(112) =  -0.1569641653843750d-15
      y(113) =   0.1856239009452638d-15
      y(114) =   0.9898466095646496d-16
      y(115) =  -0.2068030800303723d-15
      y(116) =   0.2451470336752085d-15
      y(117) =   0.9542986459336358d-16
      y(118) =  -0.2456074075580993d-15
      y(119) =   0.1532475480661800d-15
      y(120) =  -0.1229326332276474d-15
      y(121) =  -0.4750000000000000d+01
      y(122) =  -0.4750000000000001d+01
      y(123) =  -0.4750000000000000d+01
      y(124) =  -0.4750000000000000d+01
      y(125) =  -0.4750000000000000d+01
      y(126) =  -0.4750000000000000d+01
      y(127) =  -0.4750000000000000d+01
      y(128) =  -0.4750000000000000d+01
      y(129) =  -0.4750000000000000d+01
      y(130) =  -0.4750000000000000d+01
      y(131) =  -0.4750000000000001d+01
      y(132) =  -0.4750000000000001d+01
      y(133) =  -0.4750000000000000d+01
      y(134) =  -0.4750000000000000d+01
      y(135) =  -0.4750000000000000d+01
      y(136) =  -0.4750000000000000d+01
      y(137) =  -0.4749999999999999d+01
      y(138) =  -0.4750000000000000d+01
      y(139) =  -0.4750000000000000d+01
      y(140) =  -0.4750000000000000d+01
      y(141) =  -0.3537526598492654d-19
      y(142) =   0.2338193888161182d-18
      y(143) =  -0.3267771993164953d-18
      y(144) =   0.2915679914072042d-18
      y(145) =   0.1965183195887647d-18
      y(146) =  -0.6224992924096233d-19
      y(147) =  -0.1715878416756298d-18
      y(148) =  -0.2704741705248803d-18
      y(149) =   0.3008700893194513d-18
      y(150) =  -0.2703121624910402d-18
      y(151) =   0.4243755291982164d-18
      y(152) =   0.2862063003949612d-18
      y(153) =   0.1222125408406218d-19
      y(154) =  -0.4958862706817728d-18
      y(155) =  -0.7070673036251212d-18
      y(156) =  -0.4454983024194383d-18
      y(157) =  -0.1125384872521777d-18
      y(158) =   0.1512898724592511d-18
      y(159) =  -0.6163704221424137d-19
      y(160) =   0.6255426995473074d-19
      return
      end
