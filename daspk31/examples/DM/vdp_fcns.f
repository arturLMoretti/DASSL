*--------------------------------------------------------------------------*
* Function subroutine evaluating with change in Y(1) only
*--------------------------------------------------------------------------*
      SUBROUTINE VDP_FCN_Y1(T,Y,YP,YDEP,RESULT,RP)
      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846)
      DOUBLE PRECISION alpha
      DOUBLE PRECISION omega
      DOUBLE PRECISION y1
      DOUBLE PRECISION y2
      DOUBLE PRECISION T,Y,YP,YDEP,RESULT,RP
      DIMENSION Y(*),YP(*),RESULT(*),RP(*)
*
* Equation set-up. 
*


      alpha = RP(1)
      omega = RP(2)

      y1 = YDEP
      y2 = Y(2)


      RESULT(1) =y2
      RESULT(2) =-alpha*(y1**2.0000000000d+000-1.0000000000d+000)*y2-
     +   omega**2.0000000000d+000*y1
      RETURN
      END
*--------------------------------------------------------------------------*
* Function subroutine evaluating with change in Y(2) only
*--------------------------------------------------------------------------*
      SUBROUTINE VDP_FCN_Y2(T,Y,YP,YDEP,RESULT,RP)
      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846)
      DOUBLE PRECISION alpha
      DOUBLE PRECISION omega
      DOUBLE PRECISION y1
      DOUBLE PRECISION y2
      DOUBLE PRECISION T,Y,YP,YDEP,RESULT,RP
      DIMENSION Y(*),YP(*),RESULT(*),RP(*)
*
* Equation set-up. 
*


      alpha = RP(1)
      omega = RP(2)

      y1 = Y(1)
      y2 = YDEP


      RESULT(1) =y2
      RESULT(2) =-alpha*(y1**2.0000000000d+000-1.0000000000d+000)*y2-
     +   omega**2.0000000000d+000*y1
      RETURN
      END
*--------------------------------------------------------------------------*
* Function subroutine. 
*--------------------------------------------------------------------------*
      SUBROUTINE VDP_FCN(T,Y,YP,RESULT,RP)
      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979323846)
      INTEGER NEQ
      PARAMETER (NEQ=2)
      DOUBLE PRECISION T,Y,YP,RESULT,RP
      DOUBLE PRECISION alpha
      DOUBLE PRECISION omega
      DOUBLE PRECISION y1
      DOUBLE PRECISION y2
      DIMENSION Y(NEQ),YP(NEQ),RESULT(*),RP(*)
*
* Equation set-up. 
*


      alpha = RP(1)
      omega = RP(2)

      y1 = Y(1)
      y2 = Y(2)


      RESULT(1) =y2
      RESULT(2) =-alpha*(y1**2.0000000000d+000-1.0000000000d+000)*y2-
     +   omega**2.0000000000d+000*y1
      RETURN
      END

*
*--------------------------------------------------------------------------*
* Problem subroutine. 
*--------------------------------------------------------------------------*
      SUBROUTINE VDP_RES(T,Y,YP,CJ,DELTA,IRES,RP,IP,SENPAR)
      IMPLICIT NONE
      DOUBLE PRECISION RESULT
      DOUBLE PRECISION T, Y, YP, CJ, DELTA, RP, SENPAR
      INTEGER IRES, IP
      INTEGER NEQ
      PARAMETER (NEQ=2)
      DIMENSION Y(NEQ),YP(NEQ),DELTA(NEQ),RP(*),IP(*),SENPAR(*)
      DIMENSION RESULT(2)

      CALL VDP_FCN(T,Y,YP,RESULT,RP)
      DELTA(1) = YP(1) - RESULT(1)
      DELTA(2) = YP(2) - RESULT(2)
      RETURN
      END
