      SUBROUTINE adresweb (
     *           t, adY, adYP, CJ, delta, IRES, RPAR,IPAR,SENPAR, adi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ady(*), adyp(*), DELTA(*), RPAR(*), IPAR(*),SENPAR(*)
      dimension adi(*)
      dimension adrpar(1000), addelta(1000), adsenpar(2)
      PARAMETER (MAXS = 2)
      COMMON /PPAR2/ NP, NS, AX, AY, ACOEF(MAXS,MAXS), BCOEF(MAXS),
     1               DX, DY, FPI, DIFF(MAXS),
     2               COX(MAXS), COY(MAXS), MX, MY, MXNS
      
      NEQ = NS*MX*MY

      if (ires .eq. 0) then
         do i = 1, neq
            adrpar(i) = 0.0d0
            delta(i) = 0.0d0
            addelta(i) = ady(i)
         end do

         call adres(adi, rpar, senpar, delta, addelta, adrpar)

         DO 30 JY = 1,MY
            IYOFF = MXNS*(JY-1)
            DO 20 JX = 1,MX
               IC0 = IYOFF + NS*(JX-1)
               DO 10 I = 1,NS
                  IF (I .le. NP) THEN
                     ICI = IC0 + I
                     DELTA(ICI) = adyp(ici) + delta(ici)
                  ENDIF
 10            CONTINUE
 20         CONTINUE
 30      CONTINUE
      else if (ires .eq. 3) then
c
c..   compute \dot c = ady*F_p      
         do i = 1, neq
            addelta(i) = ady(i)
         enddo
         do i = 1, 2
            delta(neq+i) = 0.0d0
         end do
         
         call p_res(t, adi, adi(ny+1), cj, ires, rpar, ipar, senpar, 
     *        delta(neq+1), addelta)

         do i = 1, 2
            delta(neq +i) = adyp(neq+i) - delta(neq+i)
         end do
      end if

      return
      end
