C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      subroutine adpfweb( cc, adpcrate, adprpar, adpsenpar )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
      implicit none

C==============================================
C define parameters                            
C==============================================
      integer maxs
      parameter ( maxs = 2 )

C==============================================
C define common blocks
C==============================================
      common /ppar2/ np, ns, ax, ay, acoef, bcoef, dx, dy, fpi, diff, 
     $cox, coy, mx, my, mxns
      double precision acoef(maxs,maxs)
      double precision ax
      double precision ay
      double precision bcoef(maxs)
      double precision cox(maxs)
      double precision coy(maxs)
      double precision diff(maxs)
      double precision dx
      double precision dy
      double precision fpi
      integer mx
      integer mxns
      integer my
      integer np
      integer ns

C==============================================
C define arguments
C==============================================
      double precision adpcrate(*)
      double precision adprpar(*)
      double precision adpsenpar(*)
      double precision cc(*)

C==============================================
C define local variables
C==============================================
      integer i
      integer ic
      integer ici
      integer iyoff
      integer jx
      integer jy

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      do jy = my, 1, -1
        iyoff = mxns*(jy-1)
        do jx = mx, 1, -1
          ic = iyoff+ns*(jx-1)+1
          do i = 1, ns
            ici = ic+i-1
            adprpar(ici) = adprpar(ici)+adpcrate(ici)
            adpcrate(ici) = 0.d0
          end do
          call adpwebr( jx,jy,cc(ic),adprpar(ic),adpsenpar )
        end do
      end do

      end


      subroutine adpres( u, adpdelta, adprpar, adpsenpar )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
      implicit none

C==============================================
C define parameters                            
C==============================================
      integer maxs
      parameter ( maxs = 2 )

C==============================================
C define common blocks
C==============================================
      common /ppar2/ np, ns, ax, ay, acoef, bcoef, dx, dy, fpi, diff, 
     $cox, coy, mx, my, mxns
      double precision acoef(maxs,maxs)
      double precision ax
      double precision ay
      double precision bcoef(maxs)
      double precision cox(maxs)
      double precision coy(maxs)
      double precision diff(maxs)
      double precision dx
      double precision dy
      double precision fpi
      integer mx
      integer mxns
      integer my
      integer np
      integer ns

C==============================================
C define arguments
C==============================================
      double precision adpdelta(*)
      double precision adprpar(*)
      double precision adpsenpar(*)
      double precision u(*)

C==============================================
C define local variables
C==============================================
      integer i
      integer ic0
      integer ici
      integer iyoff
      integer jx
      integer jy

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      do jy = my, 1, -1
        iyoff = mxns*(jy-1)
        do jx = mx, 1, -1
          ic0 = iyoff+ns*(jx-1)
          do i = ns, 1, -1
            ici = ic0+i
            if (i .gt. np) then
              adpdelta(ici) = -adpdelta(ici)
            else
              adpdelta(ici) = -adpdelta(ici)
            endif
          end do
        end do
      end do
      call adpfweb( u,adpdelta,adprpar,adpsenpar )

      end


      subroutine adpwebr( jx, jy, c, adpcrate, adpsenpar )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
      implicit none

C==============================================
C define parameters                            
C==============================================
      integer maxs
      parameter ( maxs = 2 )

C==============================================
C define common blocks
C==============================================
      common /ppar2/ np, ns, ax, ay, acoef, bcoef, dx, dy, fpi, diff, 
     $cox, coy, mx, my, mxns
      double precision acoef(maxs,maxs)
      double precision ax
      double precision ay
      double precision bcoef(maxs)
      double precision cox(maxs)
      double precision coy(maxs)
      double precision diff(maxs)
      double precision dx
      double precision dy
      double precision fpi
      integer mx
      integer mxns
      integer my
      integer np
      integer ns

C==============================================
C define arguments
C==============================================
      double precision adpcrate(*)
      double precision adpsenpar(*)
      double precision c(*)
      integer jx
      integer jy

C==============================================
C define local variables
C==============================================
      double precision adpalph
      double precision adpbeta
      double precision adpfac
      integer i
      double precision x
      double precision y

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adpalph = 0.d0
      adpbeta = 0.d0
      adpfac = 0.d0

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      y = real(jy-1)*dy
      x = real(jx-1)*dx
      do i = 1, ns
        adpfac = adpfac+adpcrate(i)*c(i)*bcoef(i)
        adpcrate(i) = adpcrate(i)*c(i)
      end do
      adpalph = adpalph+adpfac*x*y
      adpbeta = adpbeta+adpfac*sin(fpi*x)*sin(fpi*y)
      adpfac = 0.d0
      do i = 1, ns
        adpcrate(i) = 0.d0
      end do
      adpsenpar(2) = adpsenpar(2)+adpbeta
      adpbeta = 0.d0
      adpsenpar(1) = adpsenpar(1)+adpalph
      adpalph = 0.d0

      end


