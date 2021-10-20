C234567
      subroutine slider(n,vel,pos,f,am,g,gp,t)
      double precision vel(n),pos(n),am(100,100),gp(22,22),
     *  t,g(22),f(n),m1,m2,m3,l1,l2,j1,j2,pi,gv,omega,tm
      pi=4.0*atan(1.0)
      m1=0.36D0
      m2=0.151104D0
      m3=0.075552D0
      l1=0.15D0
      l2=0.3D0
      j1=0.002727D0
      j2=0.0045339259D0
      gv=9.81
      omega=0.2
      tm=0.3
      am(1,1)=j1+m2*l1*l1
      am(2,1)=0.5*l1*l2*m2*cos(pos(1)-pos(2))
      am(1,2)=am(2,1)
      am(2,2)=j2
      am(3,3)=m3
      am(1,3)=0.0
      am(3,1)=0.0
      am(2,3)=0.0
      am(3,2)=0.0
      g(1)=l1*sin(pos(1))+l2*sin(pos(2))
      g(2)=pos(3)-l1*cos(pos(1))-l2*cos(pos(2))
      gp(1,1)=l1*cos(pos(1))
      gp(1,2)=l2*cos(pos(2))
      gp(1,3)=0.0
      gp(2,1)=l1*sin(pos(1))
      gp(2,2)=l2*sin(pos(2))
      gp(2,3)=1.0
      f(1)=-0.5*l1*gv*(m1+2.0*m2)*cos(pos(1))
     *  -0.5*l1*l2*m2*vel(2)*vel(2)*sin(pos(1)-pos(2))
      if (t .le. tm) then
         f(1)=f(1)+omega/tm*(1.0-cos(2*pi*t/tm))
      endif
      f(2)=-0.5*l2*gv*m2*cos(pos(2))
     *     +0.5*l1*l2*m2*vel(1)*vel(1)*sin(pos(1)-pos(2))
      f(3)=0.0
      return 
      end 
