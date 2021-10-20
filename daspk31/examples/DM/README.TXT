DASPKROOT/examples/DM/README.TXT
March 9, 2007
Stephanie R. Taylor, Linda R. Petzold
University of California at Santa Barbara

--------
OVERVIEW
--------
The DM directory contains example code that uses DDASPK and 
DDASPK_ADJOINT without the Adifor library linked in. Some 
of the examples use finite difference approximations for
the requisite derivatives. Since it is more accurate to use
analytical derivatives, DASPK has been designed to work with
automatic differentiation tools - originally this meant Adifor
and TAMC. It has been updated to use a more recent package - 
Tapenade. Tapenade is a freely available package produced by
INRIA. 

--------
EXAMPLES
--------
The examples fall into three categories:
1. Forward Mode (uses no automatic differentiation engine)
2. Adjoint Mode with Sensitivity Analysis (uses Tapenade
   as the automatic differentiation engine)
3. Adjoint Mode for Retrieving the Adjoint Variable (uses Tapenade
   as the automatic differentiation engine)
   
---------------------------
RUNNING ON CYGWIN AND LINUX
---------------------------
To run on Cygwin, edit Makefile.cygwin to contain the proper pathnames 
and compiler names/flags. Then copy Makefile.cygwin to Makefile. To run 
on Linux, edit Makefile.linuz to contain the proper pathnames and compiler 
names/flags. Then copy Makefile.cygwin Makefile. Type "make" and you will 
be provided with the list of examples.
   
--------
TAPENADE
--------
Installing Tapenade: We recommend downloading a local version of Tapenade.
Go to http://www-sop.inria.fr/tropics/tapenade.html and familiarize
yourself with the tool and its operation. Download and expand the
tarball. For use on Cygwin, you may need to remove the file
TAPDIR/bin/windows/cygwin1.dll so that is doesn't interfere with
the .dll file associated with your installation of Cygwin. Additionally,
you may need to change the permissions of TAPDIR/bin/windows/fortranParser.exe
to make it executable. We have tested the code in this directory with 
version 2.2.1. Warning: Earlier versions have slightly different directory 
structures and our Makefile will not properly call Tapenade.

Routines output from Tapenade must have wrapper routines in order for them
to be used properly by DASPK_ADJOINT. The header comments for DASPK_ADJOINT
contain some examples. The general rules are
1) Tapenade will overwrite the independent variable and any dependent
   variables, so make sure you protect variables that you don't want
   over-written. For example, this means both DELTA and ADDELTA should be 
   protected in RES_ADY.
2) Tapenade produces code Xb = JT.Yb (where Xb and Yb are the backward
   differential variables and JT is the transpose ot the system Jacobian
   evaluated at X). However, DASPK_ADJOINT is expected Xb=Jt.Yb+Xb,
   so the wrapper must save Xb and then add it to the result of the
   call to the Tapenade-generated routine.
   
Warning: Be sure to properly size all variables. In the examples below, all 
variables happen to be of the same size (2), but this is not generally
the case. The sizes of the variables in the wrapper must match those
of the corresponding variable in the generated routine.

Example: 
Given a residual routine RES, 
C Use Tapenade to create the gradient of RES with respect to Y.
C                 Suppose the function RES (as described in the DASPK3.1 docs)
C                 resides in a file odes.f, then we will use tapenade in
C                 gradient (reverse) mode, to produce ADDELTA*F_y. We set
C                 "Y" as the independent variable (-vars "Y") and
C                 "DELTA" as the dependent variable (-outvars "DELTA"):
C
C                 java -mx256m -classpath \
C                   "C:\tapenade2.1.3\jars\tapenade.jar" \
C                   -Dali.tapenade_home="C:\TAPENA~1.3" \
C                   -Dali.browser="dummy" \
C                   topLevel.Differentiator \
C                   -backward -head RES -vars "Y" -outvars "DELTA" \
C                   "C:\DASPK3P1\examples\odes.f"
C                 
C                 The result will be a file res_b.f containing the function
C
C                 SUBROUTINE RES_B(t, y, yb, yp, cj, delta, deltab, ires, rp, ip, 
C     +                 senpar)
C
C Write the wrapper
C
C                 Even though DELTAB is the input, it is also treated as output.
C                 The implication is that if we pass ADDELTA into RES_B, it will
C                 be clobbered. Additionally, since RES presumably alters DELTA, 
C                 we must assume RES_B alters DELTA. This means RES_ADY must be a 
C                 wrapper for RES_B that protects both DELTA and ADDELTA.  For this 
C                 particular example, we have a two-dimensional system, so we dimension
C                 DELTA and DELTAB to be 2. Additionally, we need to save the value of 
C                 ADY, and add it to the result of RES_B.
C
C                 SUBROUTINE RES_ADY(T, Y, ADY, YPRIME, CJ, DELTA, ADDELTA,
C     *              IRES, RPAR, IPAR, SENPAR)
C                 IMPLICIT NONE
C                 DOUBLE PRECISION T, Y, ADY, YPRIME, CJ, DELTA, ADDELTA
C                 DOUBLE PRECISION RPAR, SENPAR
C                 INTEGER IRES,IPAR,II
C                 DIMENSION Y(*),YPRIME(*)ADY(*),DELTA(*),ADDELTA(*)
C                 DIMENSION RPAR(*),SENPAR(*),IPAR(*)
C                 DOUBLE PRECISION DELTAB(2), DELTA_COPY(2), ADY_SAV(2)
C          
C                 DO II = 1, 2
C                   DELTAB(II) = ADDELTA(II)
C                   DELTA_COPY(II) = DELTA(II)
C                   ADY_SAV(II) = ADY(II)
C                 END DO
C                 CALL RES_B(T, Y, ADY, YPRIME, CJ, DELTA_COPY, DELTAB,
C     +                 IRES, RP, IP, SENPAR)     
C                 DO II = 1, 2
C                   ADY(II) = ADY(II)+ADY_SAV(II)
C                 END DO
C                 RETURN
C                 END

--------
TAMC
--------
For those who may already have TAMC, there is a subdirectory in 
DM called TAMC which contains examples using TAMC.

Aside from cosmetic differences in the code, there is one fundamental
difference between Tapenade and TAMC code generated in reverse mode. 
For function Y = f(X), TAMC generates code
    Xb = (df/dX)^{T}.Yb + Xb (where Xb and Yb are the differentials in backward mode)
whereas Tapenade generates code
    Xb = (df/dX)^{T}.Yb
To make up for the difference, wrappers for Tapenade-generated
routines must save the values of the independent variable (Xb), 
call the Tapenade-generated code, and then add the saved
value onto the result.

   
