      subroutine gauss(n,a,b,x)
      implicit none
      integer n,i,j,k
      double precision a(1000,1000),s(100)
      double precision b(1000),x(1000)
      double precision smax,rmax,sum,pk,r,z
      integer p(1000)
      double precision max
c
C
      do 3 i=1,n
         p(i) = i
         smax = 0.0 
         do 2 j=1,n
            smax = max(smax,abs(a(i,j)))
 2       continue
         s(i) = smax
 3    continue
c
      do 7 k=1,n-1
         rmax = 0.0
         do 4 i=k,n
            r = abs(a(p(i),k))/s(p(i))
            if (r .gt. rmax) then
               j = i
               rmax = r
            endif
 4       continue
c
         pk = p(j)
         p(j) = p(k)
         p(k) = pk
c
         do 6 i=k+1,n      
            z = a(p(i),k)/a(p(k),k)       
            a(p(i),k) = z
            do 5 j=k+1,n    
               a(p(i),j) = a(p(i),j) - z*a(p(k),j)      
 5          continue
 6       continue  
 7    continue    
c
      do 9 k=1,n-1
         do 8 i=k+1,n
            b(p(i)) = b(p(i)) - a(p(i),k)*b(p(k))
 8       continue
 9    continue
      do 11 i=n,1,-1
         sum = b(p(i))
         do 10 j=i+1,n
            sum = sum - a(p(i),j)*x(j)
 10      continue
         x(i) = sum/a(p(i),i)
 11   continue
c
c
      return
      end 
