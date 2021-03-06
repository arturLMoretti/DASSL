C Copyright 2000 the Regents of the University of California.
C All rights reserved.
C
      SUBROUTINE DDASPK (
     *     RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL,ATOL,
     *     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL, 
     *     SENPAR, G_RES, K_RES, T_RES, A_RES)
C***VERSION 3.0
C***BEGIN PROLOGUE  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  910624   
C***REVISION DATE  920929   (CJ in RES call, RES counter fix.)
C***REVISION DATE  921215   (Warnings on poor iteration performance)
C***REVISION DATE  921216   (NRMAX as optional input)
C***REVISION DATE  930315   (Name change: DDINI to DDINIT)
C***REVISION DATE  940822   (Replaced initial condition calculation)
C***REVISION DATE  941101   (Added linesearch in I.C. calculations)
C***REVISION DATE  941220   (Misc. corrections throughout)
C***REVISION DATE  950125   (Added DINVWT routine)
C***REVISION DATE  950714   (Misc. corrections throughout)
C***REVISION DATE  950802   (Default NRMAX = 5, based on tests.)
C***REVISION DATE  950808   (Optional error test added.)
C***REVISION DATE  950814   (Added I.C. constraints and INFO(14))
C***REVISION DATE  950828   (Various minor corrections.)
C***REVISION DATE  951006   (Corrected WT scaling in DFNRMK.)
C***REVISION DATE  960129   (Corrected RL bug in DLINSD, DLINSK.)
C***REVISION DATE  960301   (Added NONNEG to SAVE statement.)
C***REVISION DATE  980813   (Added sensitivity analysis capability.)
C***REVISION DATE  981001   (Added staggered corrector method.)
C***REVISION DATE  981018   (Added I.C. calculation for index-2.)
C***REVISION DATE  981022   (Added ADIFOR processing routines.)
C***REVISION DATE  981103   (Added SENPAR and G_RES in arguments list.)
C***REVISION DATE  981105   (Parallel implementation for sensitivity.)
C***REVISION DATE  981112   (Staggered corrector for initial conditions.)
C***REVISION DATE  981117   (Two-step processing for index-two initialization.)
C***REVISION DATE  981208   (Modified I.C. calculation for index-2.)
C***REVISION DATE  981222   (Added matrix times vector method.)
C***REVISION DATE  990123   (Modified I.C. calculation for index-2.)
C***REVISION DATE  990125   (Reduced work space.)
C***REVISION DATE  990215   (Remove LNNIS index.)
C***REVISION DATE  990219   (Modified norm for initial conditions.)
C***REVISION DATE  990303   (Fixed the FLOAN=NY for Krylov method.)
C***REVISION DATE  990304   (Fixed a bug in simultaneous corrector method.)
C***REVISION DATE  990305   (Fixed a bug in Krylov work space.)
C***REVISION DATE  990306   (Fixed a bug in the IWORK array about IWP.)
C***REVISION DATE  990308   (Added staggered direct method.)
C***REVISION DATE  990309   (Added counter LNLIS=IWORK(38).)
C***REVISION DATE  990508   (Bugfixed for initialization of Krylov method.)
C***REVISION DATE  990714   (Bugfixed for for TN=TOUT.)
C***REVISION DATE  990723   (Remove some static variables.)
C***REVISION DATE  990723   (Remove bugs for 'warm restart'.)
C***REVISION DATE  990811   (Modification for ADIFOR with SparsLinC option.)
C***REVISION DATE  991105   (Fixed a bug for initialization of Krylov method.)
C***REVISION DATE  991122   (Fixed a bug for the staggered direct method.)
C***REVISION DATE  991207   (Added adjoint method.)
C***REVISION DATE  000207   (Various minor corrections for adjoint method.)
C***REVISION DATE  000212   (Added options for quadrature calculation.)
C***REVISION DATE  000217   (Combine two corrector methods into one.)
C***REVISION DATE  000223   (Various minor corrections for quadrature.)
C***REVISION DATE  000224   (Added new convergence rate for sensitivities.)
C***REVISION DATE  000225   (Added dummy routine A_RES for adjoint method.)
C***REVISION DATE  000228   (Modified ADIFOR with SparsLinc options.)
C***REVISION DATE  000301   (Adjoint sensitivity for \int g case.)
C***REVISION DATE  000302   (Modified error control for quadrature variables.)
C***REVISION DATE  000306   (Split error test for quadrature variables.)
C***REVISION DATE  000307   (Bugfixed for TN=TSTOP.)
C***REVISION DATE  000308   (Bugfixed for Krylov method when INFO(30)!=0.)
C***REVISION DATE  000310   (Bugfixed for counter NRE.)
C***REVISION DATE  000315   (Bugfixed for when LIW is exact IWORK(LNIW).)
C***REVISION DATE  000317   (Minor modification for Jacobian evaluation.)
C***REVISION DATE  000321   (New residual evaluation of adjoint equations.)
C***REVISION DATE  000331   (Added initialization for the ADIFOR space.)
C***REVISION DATE  000404   (Modified initialization for index-2 system.)
C***REVISION DATE  000406   (Added AD option K_RES for Krylov Jv product.)
C***REVISION DATE  000410   (Added AD option T_RES for index-2 initialization.)
C***REVISION DATE  000411   (Bugfixed for quadrature variable integration.)
C***REVISION DATE  000425   (Fixed output message when TSTOP < TN.)
C***REVISION DATE  000426   (Bugfixed for residuals of adjoint equations.)
C***REVISION DATE  000512   (Added initialization for index-2 adjoint system.)
C***REVISION DATE  000515   (Bugfixed in the cubic Hermit interpolation.)
C***REVISION DATE  000518   (Minor modification for index-2 adjoint system.)
C***REVISION DATE  000525   (Rewrite the documentation.)
C***REVISION DATE  000530   (Rewrite Jacobian evaluation by ADIFOR.)
C***REVISION DATE  000531   (Bugfixed for the adjoint integration.)
C***REVISION DATE  000620   (Bugfixed for the index-2 adjoint integration.)
C***REVISION DATE  000621   (Added LTFN=19 for adjoint integration.)
C***REVISION DATE  000718   (Minor modification on Newton iteration.)
C***REVISION DATE  000724   (Minor modification for index-0 initialization.)
C***REVISION DATE  000801   (Removed IWORK(38).)
C***REVISION DATE  000829   (Rewrote Krylov method for adjoint method.)
C***REVISION DATE  000831   (Bugfixed for quadrature variables in adjoint)
C***REVISION DATE  000908   (Added case F_y'(t,y) in adjoint method.)
C***REVISION DATE  000919   (Modified convergence rate for sensitivity.)
C***REVISION DATE  001113   (Added output variable (y_o = f(y)) options.)
C***REVISION DATE  001122   (Added user-input direct solver INFO(6)=2 option.)
C***REVISION DATE  001201   (Added ICOPT and ID in JAC arguments.)
C***REVISION DATE  001204   (Added matrix-vector product option INFO(5)=3.)
C***REVISION DATE  001213   (Added INFO(11) = 6 for initialization.)
C***REVISION DATE  001215   (Modified initialization in second pass.)
C***REVISION DATE  010302   (Fixed a bug for backward integration.)
C***REVISION DATE  020511   (Bug fixed for INFO(5)=2 and INFO(12)=1)
C***REVISION DATE  020722   (Bug fixed for INFO(24) < -1)
C***REVISION DATE  020730   (Bug fixed for INFO(20) = 1)
C***REVISION DATE  020730   (Documentation fixed for GRADIENT(s))
C***REVISION DATE  020801   (Changed argument list for seed matrix option)
C***REVISION DATE  020801   (Changed argument list for adjoint routines)
C***REVISION DATE  020801   (Added ADIFOR3.0 script file in documentation)  
C***REVISION DATE  020829   (Modification based on new DASPK2.0)  
C***REVISION DATE  030908   (Fixed a memory leak for ADIFOR SP option)  
C***REVISION DATE  070131   (Modification for ADIFOR SP option)  
C***REVISION DATE  070611   (moved init of ierj and iernew to beginning of DDASID)
C***REVISION DATE  070611   (Bug fixed for documentation of number of quad vars (should be INFO(28) but was described as INFO(30)).
C***CATEGORY NO.  I1A2        
C***KEYWORDS  DIFFERENTIAL/ALGEBRAIC, BACKWARD DIFFERENTIATION FORMULAS,
C             IMPLICIT DIFFERENTIAL SYSTEMS, KRYLOV ITERATION
C***AUTHORS   Linda R. Petzold, Peter N. Brown, Alan C. Hindmarsh and
C             Shengtai Li
C
C             Addresses:
C             Linda Petzold and Shengtai Li
C             Department of Computer Science
C             University of California
C             Santa Barbara, CA  93106
C             Email:  petzold@engineering.ucsb.edu
C                    shengtai@engineering.ucsb.edu
C
C             Peter Brown and Alan Hindmarsh
C             Center for Applied Scientific Computing
C             Lawrence Livermore National Laboratory
C             P.O. Box 808,
C             Livermore, CA 94551
C
C***PURPOSE  This code solves a system of differential/algebraic 
C            equations (DAEs) of the form 
C               G(t,y,y',p) = 0 , 
C            using a combination of Backward Differentiation Formula 
C            (BDF) methods and a choice of two linear system solution 
C            methods: direct (dense or banded) or Krylov (iterative).
C            If specified by the user, sensitivities of the solution
C            y with respect to perturbations in the parameters p are
C            also computed. The class of problems that DDASPK can solve 
C            includes index-0, index-1 and Hessenberg index-2 DAE system. 
C            This version is in double precision.
C-----------------------------------------------------------------------
C
C***This version DDASPK has been modified to perform an optional 
C   sensitivity analysis of problem parameters, to accommodate derivative  
C   evaluation via automatic differentiation (ADIFOR), to allow  
C   parallel computation of the sensitivities via MPI, to allow  
C   solution of Hessenberg index-2 DAE systems, to allow fast evaluation 
C   of quadrature, and to serve as the forward solver for DASPKADJOINT.
C
C   All modifications were made by:
C
C   Shengtai Li and Linda R. Petzold
C   Department of Computer Science
C   University of California
C   Santa Barbara, CA 93106
C-----------------------------------------------------------------------   
C***DESCRIPTION
C
C *Usage:
C
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR(*)
C      DOUBLE PRECISION T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*),
C         RWORK(LRW), RPAR(*), SENPAR(*)
C      EXTERNAL  RES, JAC, PSOL, G_RES
C
C      CALL DDASPK (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
C     *  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL, SENPAR, G_RES)
C
C  Quantities which may be altered by the code are:
C     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL, IDID, RWORK(*), IWORK(*)
C
C
C *Arguments:
C
C  RES:EXT          This is the name of a routine which you
C                   provide to define the residual function G(t,y,y')
C                   of the differential/algebraic system.
C
C  NEQ:IN           This is the number of equations in the system.
C                   If you are solving for the state variables only, NEQ
C                   is the number of state variables.  If you are
C                   also computing sensitivities, NEQ is the total
C                   number of variables, including state variables
C                   and sensitivity variables.
C
C  T:INOUT          This is the current value of the independent 
C                   variable.
C
C  Y(*):INOUT       This array contains the solution (and sensitivity) 
C                   components at T.
C
C  YPRIME(*):INOUT  This array contains the derivatives of the solution
C                   (and sensitivity) components at T.
C
C  TOUT:IN          This is a point at which a solution is desired.
C
C  INFO(N):IN       This is an integer array used to communicate details
C                   of how the solution is to be carried out, such as
C                   tolerance type, matrix structure, step size and
C                   order limits, and choice of nonlinear system method.
C                   N must be at least 30.
C
C  RTOL,ATOL:INOUT  These quantities represent absolute and relative
C                   error tolerances (on local error) which you provide
C                   to indicate how accurately you wish the solution to
C                   be computed.  You may choose them to be both scalars
C                   or else both arrays of length NEQ.
C
C  IDID:OUT         This integer scalar is an indicator reporting what
C                   the code did.  You must monitor this variable to
C                   decide what action to take next.
C
C  RWORK:WORK       A real work array of length LRW which provides the
C                   code with needed storage space.
C
C  LRW:IN           The length of RWORK.
C
C  IWORK:WORK       An integer work array of length LIW which provides
C                   the code with needed storage space.
C
C  LIW:IN           The length of IWORK.
C
C  RPAR,IPAR:IN     These are real and integer parameter arrays which
C                   you can use for communication between your calling
C                   program and the RES, JAC, and PSOL routines.
C
C  JAC:EXT          This is the name of a routine which you may
C                   provide (optionally) for calculating Jacobian 
C                   (partial derivative) data involved in solving linear
C                   systems within DDASPK. 
C                   When you use the ADIFOR option to calculate the Jacobian,
C                   this is the name of the ADIFOR-generated routine.
C
C  PSOL:EXT         This is the name of a routine which you must
C                   provide for solving linear systems if you selected
C                   a Krylov method.  The purpose of PSOL is to solve
C                   linear systems involving a left preconditioner P.
C  
C  SENPAR:IN        This is a real array for sensitivity parameters that 
C                   appear in the RES routine. If you wish to compute 
C                   sensitivities and RES depends on the problem parameters, 
C                   SENPAR must be used to store and pass the problem 
C                   parameters to RES.
C  
C  G_RES:EXT        This is the name of the ADIFOR-generated routine which 
C                   you may provide (optionally) for the evluation of the
C                   sensitivity equations.
C 
C  K_RES:EXT        This is the name of the ADIFOR-generated routine which 
C                   you may provide (optionally) for the matrix-vector product 
C                   in the Krylov iteration.
C
C  T_RES:EXT        This is the name of the ADIFOR-generated routine which 
C                   you may provide (optionally) for the initialization of 
C                   an index-2 system.
C 
C  A_RES:EXT        This is the name of a routine which is reserved only for
C                   DASPKADJOINT. The users can disregard.
C                   
C                   
C *Overview
C
C  The DDASPK solver uses the backward differentiation formulas of
C  orders one through five to solve a system of the form G(t,y,y') = 0
C  for y = Y and y' = YPRIME.  Values for Y and YPRIME at the initial 
C  time must be given as input.  These values should be consistent; 
C  that is, if T, Y, YPRIME are the given initial values, they should 
C  satisfy G(T,Y,YPRIME) = 0.  However, if consistent initial values are not
C  known, in many cases you can have DDASPK solve for them -- see INFO(11).
C  Note that if sensitivities are to be computed, the initial values
C  for the sensitivities must also be consistent.
C  (This and other options are described in more detail below.)
C
C  Normally, DDASPK solves the system from T to TOUT.  It is easy to
C  continue the solution to get results at additional TOUT.  This is
C  the interval mode of operation.  Intermediate results can also be
C  obtained easily by specifying INFO(3).
C
C  DDASPK includes an option to perform a sensitivity analysis of the
C  problem to be solved. Given a DAE depending on a vector of parameters
C  p, G(t,y,y',p) = 0, where y, y' and G are vectors of length Ny and p
C  a vector of length Np, DDASPK will compute s_i=dy/dp_i for i = 1....Np.
C  There are two types of sensitivity parameters: one set of the parameters
C  may appear only in the initial conditions; the other may also appear 
C  directly in the evaluation of the DAE (RES). During the sensitivity  
C  analysis, only the second type of parameters need to be stored. The 
C  sensitivities satisfy the equations -
C
C  dG/dy * s_i + dG/dy' * s'_i + dG/dp_i = 0, i = 1.....Np
C
C  or a finite difference approximation to this system. The sen -
C  sitivity values are stored in the Y-vector following the solution
C  of the DAE, so that -
C
C           Y = |   y    |    ,    Y' = |    y'   |
C               |  s_1   |              |   s'_1  |
C               |   .    |              |    .    |
C               |   .    |              |    .    |
C               |  s_Np  |              |   s'_Np |
C
C  DDASPK requires the first type of parameters to be stored following the 
C  second type of parameters in Y and Y'.
C
C  DDASPK includes an option to integrate quadrature variables 
C  separately from the other state variables. Quadrature variables are
C  defined as a special group of state variables that represent
C  the quadratures of functions of the other state variables, e.g.
C          y_q' = f(y),
C  where y_q is called quadrature variable. The quadrature
C  variables are always differential variables. The derivative (YPRIME) of 
C  every state variable depends only on non-quadrature state variables. 
C  Usually, the quadrature equations depend on many state variables, and 
C  they can destroy the banded-structure of the Jacobian matrix in the Newton 
C  iteration if the quadrature equations are not treated differently. 
C  DDASPK uses a staggered approach to integrate the quadrature
C  variables and exclude them from Jacobian evaluation and Newton iteration. 
C  DDASPK requires that the quadrature variables be placed in the last 
C  part of components of the state variables.
C
C  The quadrature variable option can also be used to integrate output
C  variables, e.g. 
C        y_out = f(y),
C  where y_out is called an output variable. Output variables are always
C  index-1 algebraic variables. Output variables can be treated as
C  special quadrature variables if the residual equations are properly scaled.
C  (see the description of RES in this documentation)
C
C  On each step taken by DDASPK, a sequence of nonlinear algebraic  
C  systems arises.  These are solved by one of two types of
C  methods:
C    * a Newton iteration with a direct method for the linear
C      systems involved (INFO(12) = 0), or
C    * a Newton iteration with a preconditioned Krylov iterative 
C      method for the linear systems involved (INFO(12) = 1).
C
C  The direct method choices are dense and banded matrix solvers, 
C  with either a user-supplied or an internal difference quotient 
C  Jacobian matrix, or an ADIFOR-generated Jacobian, as specified by 
C  INFO(5) and INFO(6). In the banded case, INFO(6) = 1, you must supply 
C  half-bandwidths in IWORK(1) and IWORK(2). 
C
C  The Krylov method is the Generalized Minimum Residual (GMRES) 
C  method, in either complete or incomplete form, and with 
C  scaling and preconditioning.  The method is implemented
C  in an algorithm called SPIGMR.  Certain options in the Krylov 
C  method case are specified by INFO(13) and INFO(15).
C
C  If the Krylov method is chosen, you may supply a pair of routines,
C  JAC and PSOL, to apply preconditioning to the linear system.
C  If the system is A*x = b, the matrix is A = dG/dY + CJ*dG/dYPRIME
C  (of order NEQ).  This system can then be preconditioned in the form
C  (P-inverse)*A*x = (P-inverse)*b, with left preconditioner P.
C  (DDASPK does not allow right preconditioning.)
C  The Krylov method is applied to this altered, but equivalent,
C  linear system, hopefully with much better performance than without
C  preconditioning.  (In addition, a diagonal scaling matrix based on
C  the tolerances is also introduced into the altered system.)
C
C  The JAC routine evaluates any data needed for solving systems
C  with coefficient matrix P, and PSOL carries out that solution.
C  In any case, in order to improve convergence, you should try to
C  make P approximate the matrix A as much as possible, while keeping
C  the system P*x = b reasonably easy and inexpensive to solve for x,
C  given a vector b.
C
C-------------------Dependence with other packages---------------------
C
C  *NOTE* Depending on the options selected, this package uses a few 
C         routines or functions generated by the ADIFOR automatic 
C         differentiation package and from the MPI parallel package. 
C         
C         If you do not have access to the ADIFOR package, you 
C         can still use DDASPK by setting the appropriate options and 
C         linking an
C            "adf_dummy.f" 
C         file when compiling the code; or you can use other AD-tool
C         to generate a code and write a wrapper in the required format
C         to call that code. 
C
C         If you do not have access to the MPI package, you can still 
C         use DDASPK for sequential programming by linking an
C            "mpi_dummy.f" 
C         file when compiling the code. Both dummy files are avaiable with
C         the DDASPK package. 
C
C         If you use ADIFOR option and RPAR is used to pass the 
C         information of Y or YPRIME in your RES routine, the ADIFOR-
C         generated code will contain the seed matrix of RPAR. You have to
C         write a wrapper in the required format to call that code.
C
C----------------------------------------------------------------------
C
C *Description
C
C------INPUT - WHAT TO DO ON THE FIRST CALL TO DDASPK-------------------
C
C
C  The first call to the code is defined to be the start of each new
C  problem.  Read through the descriptions of all the following items,
C  provide sufficient storage space for the designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  provide the requested information about how you want the problem to 
C  be solved.
C
C
C  RES -- Provide a routine of the form
C
C       SUBROUTINE RES(T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C
C         to define the system of differential/algebraic
C         equations which is to be solved. For the given values
C         of T, Y and YPRIME, the routine should return
C         the residual of the differential/algebraic system
C             DELTA = G(T,Y,YPRIME)
C         DELTA is a vector of length NEQ which is output from RES.
C
C         Subroutine RES must not alter T, Y, YPRIME, or CJ.
C         You must declare the name RES in an EXTERNAL
C         statement in your program that calls DDASPK.
C         You must dimension Y, YPRIME, and DELTA in RES.
C
C         The input argument CJ can be ignored, or used to rescale
C         constraint equations in the system (see Ref. 2, p. 145).
C         Note: In this respect, DDASPK is not downward-compatible
C         with DDASSL, which does not have the RES argument CJ.
C
C         IRES is an integer flag which can be equal to 0, 1, 2 or 3
C         on input.  Subroutine RES should alter IRES only if it
C         encounters an illegal value of Y or a stop condition.
C         Set IRES = -1 if an input value is illegal, and DDASPK
C         will try to solve the problem without getting IRES = -1.
C         If IRES = -2, DDASPK will return control to the calling
C         program with IDID = -11.
C
C         RPAR and IPAR are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and routine RES. They are not altered by DDASPK. If you
C         do not need RPAR or IPAR, ignore these parameters by treat-
C         ing them as dummy arguments. If you do choose to use them,
C         dimension them in your calling program and in RES as arrays
C         of appropriate length.  
C
C         SENPAR is a real parameter array for sensitivity computation. It
C         is not altered by DDASPK. If you do not need SENPAR, ignore it 
C         by treating it as dummy argument. If the finite difference or 
C         ADIFOR methods are selected to evaluate sensitivity equations and 
C         some sensitivity parameters appear in RES, SENPAR must contain 
C         the problem parameters related to RES. You must dimension it 
C         in your calling program and in RES as an array of appropriate 
C         length.  
C
C         Depending on the value of INFO(20), RES may also be used
C         to define the residuals of the sensitivity equations.
C         If INFO(20) = 2 (sensitivity equations defined analytically)
C         and if IRES = 0, RES should evaluate only the Ny residuals of
C         the original DAE. If INFO(20) = 2 and IRES = 1, RES should 
C         evaluate the residuals of of the sensitivity equations, storing 
C         them in DELTA immediately following the residuals of the state 
C         variables. If the parallel method is also selected (INFO(27)>1) 
C         for sensitivity computation, RES should evaluate the residuals 
C         of the sensitivity equations with respect to parameters starting 
C         from INFO(26)+1, and the ones every other INFO(27) as  
C         INFO(26)+1, INFO(26)+1+INFO(27),...,INFO(26)+1+j*INFO(27),... 
C         where the values of INFO(26) and INFO(27) are generated by a 
C         set-up routine of the message passing interface (MPI). 
C
C         For the initialization of a Hessenberg index-2 DAE problem, 
C         if INFO(20) = 2 and IRES = 2, RES is used to evaluate
C         the residuals of the time derivatives for index-2 constraints
C         which should include the index-2 sensitivity constraints if
C         sensitivity is considered (see the description of item INFO(11)
C         in the documentation).
C
C         Depending on the value of INFO(28), RES may also be used to 
C         define the residuals of the equations for the quadrature variables.
C         If INFO(28)--the number of quadrature variables--is greater than 0, 
C         the RES routine should evaluate only the residuals of the quadrature 
C         equations if IRES = 3. The residuals for the quadrature equations 
C         should be placed in the last part of the Ny components of the DELTA 
C         array for the state variables. If INFO(20) = 2 is also chosen, the 
C         corresponding sensitivity equations to the quadrature equations 
C         should also be evaluated. 
C
C         The INFO(28) should also include the number of output variables 
C         (i.e., it should be the sum of the number of quadrature variables 
C         and the number of output variables). The residual equations for 
C         the output variables MUST be scaled by CJ, i.e. transform
C             delta = y_out - f(y) 
C         into 
C             delta = CJ*(y_o - f(y))
C         The initial values of the output variables are specified by
C         the user, not calculated inside DDASPK even if the initialization 
C         is selected by the user. The quadrature and output variables can
C         be mixed together without any specific order.
C
C         Overall, an efficient RES routine has the following multiple-branch 
C         selection statement:
C
C            IF (IRES .EQ. 0) THEN
C               define the residuals for the state equations 
C               (excluding quadrature/output equations)
C            ELSE IF (IRES .EQ. 1) THEN
C               define the residuals for the sensitivity equations
C               (excluding quadrature/output sensitivity equations)
C            ELSE IF (IRES .EQ. 2) THEN
C               define the residuals of the time derivatives of the  
C               index-2 constraints and/or index-2 sensitivity constraints
C            ELSE IF (IRES .EQ. 3) THEN
C               define the residuals for the quadrature equations and/or
C               the corresponding sensitivity equations
C            END IF
C         
C  NEQ -- Set it to the number of equations in the system (NEQ .GE. 1).
C         This is the total number of equations, including both state
C         equations (including the quadrature and output equations) and 
C         sensitivity equations, if applicable.
C
C  T -- Set it to the initial point of the integration. (T must be
C       a variable.)
C
C  Y(*) -- Set this array to the initial values of the NEQ solution and 
C          sensitivity components at the initial point.  You must 
C          dimension Y of length at least NEQ in your calling program.
C
C  YPRIME(*) -- Set this array to the initial values of the NEQ first
C               derivatives of the solution and sensitivity components
C               at the initial point. You must dimension YPRIME at least 
C               NEQ in your calling program. 
C
C  TOUT - Set it to the first point at which a solution is desired.
C         You cannot take TOUT = T.  Integration either forward in T
C         (TOUT .GT. T) or backward in T (TOUT .LT. T) is permitted.
C
C         The code advances the solution from T to TOUT using step
C         sizes which are automatically selected so as to achieve the
C         desired accuracy.  If you wish, the code will return with the
C         solution and its derivative at intermediate steps (the
C         intermediate-output mode) so that you can monitor them,
C         but you still must provide TOUT in accord with the basic
C         aim of the code.
C
C         The first step taken by the code is a critical one because
C         it must reflect how fast the solution changes near the
C         initial point.  The code automatically selects an initial
C         step size which is practically always suitable for the
C         problem.  By using the fact that the code will not step past
C         TOUT in the first step, you could, if necessary, restrict the
C         length of the initial step.
C
C         For some problems it may not be permissible to integrate
C         past a point TSTOP, because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         TSTOP.  When you have declared a TSTOP point (see INFO(4)
C         and RWORK(1)), you have told the code not to integrate past
C         TSTOP.  In this case any TOUT beyond TSTOP is invalid input.
C
C  INFO(*) - Use the INFO array to give the code more details about
C            how you want your problem solved.  This array should be
C            dimensioned of length at least 30. DDASPK uses only the 
C            first 30 entries.  You must respond to all of the following
C            items, which are arranged as questions.  The simplest use
C            of DDASPK corresponds to setting all entries of INFO to 0.
C
C       INFO(1) - This parameter enables the code to initialize itself.
C              You must set it to indicate the start of every new 
C              problem.
C
C          **** Is this the first call for this problem ...
C                yes - set INFO(1) = 0
C                 no - not applicable here.
C                      See below for continuation calls.  ****
C
C       INFO(2) - How much accuracy you want of your solution
C              is specified by the error tolerances RTOL and ATOL.
C              The simplest use is to take them both to be scalars.
C              To obtain more flexibility, they can both be arrays.
C              The code must be told your choice.
C
C          **** Are both error tolerances RTOL, ATOL scalars ...
C                yes - set INFO(2) = 0
C                      and input scalars for both RTOL and ATOL
C                 no - set INFO(2) = 1
C                      and input arrays of length NEQ for both 
C                      RTOL and ATOL ****
C
C       INFO(3) - The code integrates from T in the direction of TOUT
C              by steps.  If you wish, it will return the computed
C              solution and derivative at the next intermediate step
C              (the intermediate-output mode) or TOUT, whichever comes
C              first.  This is a good way to proceed if you want to
C              see the behavior of the solution.  If you must have
C              solutions at a great many specific TOUT points, this
C              code will compute them efficiently.
C
C          **** Do you want the solution only at
C               TOUT (and not at the next intermediate step) ...
C                yes - set INFO(3) = 0
C                 no - set INFO(3) = 1 ****
C
C       INFO(4) - To handle solutions at a great many specific
C              values TOUT efficiently, this code may integrate past
C              TOUT and interpolate to obtain the result at TOUT.
C              Sometimes it is not possible to integrate beyond some
C              point TSTOP because the equation changes there or it is
C              not defined past TSTOP.  Then you must tell the code
C              this stop condition.
C
C         **** Can the integration be carried out without any
C              restrictions on the independent variable T ...
C                 yes - set INFO(4) = 0
C                  no - set INFO(4) = 1
C                       and define the stopping point TSTOP by
C                       setting RWORK(1) = TSTOP ****
C
C       INFO(5) - Options for partial derivative evaluation methods.
C              The numerical methods used in DDASPK make use of a matrix 
C              of partial derivatives of the system of differential- 
C              algebraic equations.  If you do not provide a routine 
C              to evaluate this matrix analytically (see description of 
C              the item JAC in the call list), it will be approximated 
C              by numerical differencing. Although it is less trouble for 
C              you to have DDASPK compute partial derivatives by 
C              numerical differencing, the solution will be more reliable 
C              if you provide the derivatives via JAC.  Usually numerical 
C              differencing is more costly than evaluating derivatives in 
C              JAC, but sometimes it is not - this depends on your problem.
C
C              This option is also used by the Krylov iterative method
C              in evaluating the matrix-vector product by K_RES (see the
C              documentation for K_RES for details).
C
C              It should be noted that the Jacobian matrix here is only
C              with respect to the state variables (excluding the 
C              quadrature/output and sensitivity variables).
C
C           **** Do you want the code to evaluate the partial deriv-
C                atives automatically by numerical differences ...
C                 yes - set INFO(5) = 0
C                  no - There are three options
C                    1. INFO(5) = 1. If the direct method (INFO(12)=0) is 
C                       chosen, provide routine JAC for evaluating the
C                       matrix of partial derivatives of G with respect 
C                       to Y. 
C
C                       If the Krylov method (INFO(12)=1) is chosen, 
C                       provide routine K_RES for evaluating the
C                       matrix-vector product.
C
C                    2. INFO(5) = 2. If the direct method (INFO(12)=0) is 
C                       chosen, provide the ADIFOR-generated routine in JAC,
C                       where ADIFOR with the SparsLinC option should be used.
C
C                       If the Krylov method (INFO(12)=1) is chosen, provide 
C                       the ADIFOR-generated routine in K_RES, where ADIFOR 
C                       with matrix-vector product option should be used
C.
C                    3. INFO(5) = 3. If the direct method (INFO(12)=0) is 
C                       chosen, provide the ADIFOR-generated routine in JAC,
C                       where ADIFOR with the seed matrix option should be 
C                       used.
C
C                       If the Krylov method (INFO(12)=1) is chosen, provide 
C                       routine K_RES for evaluating the matrix-vector product.
C                       (See documentation of K_RES for more information)
C                   
C
C       INFO(6) - used only when INFO(12) = 0 (direct methods).
C              DDASPK will perform much better if the matrix of
C              partial derivatives, dG/dY + CJ*dG/dYPRIME (here CJ is
C              a scalar determined by DDASPK), is banded and the code
C              is told this.  In this case, the storage needed will be
C              greatly reduced, numerical differencing will be performed
C              much cheaper, and a number of important algorithms will
C              execute much faster.  The differential equation is said 
C              to have half-bandwidths ML (lower) and MU (upper) if 
C              equation i involves only unknowns Y(j) with
C                             i-ML .le. j .le. i+MU ,
C              for all i=1,2,...,NEQ.  Thus, ML and MU are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded.  If you do not
C              indicate that the equation has a banded matrix of partial
C              derivatives, the code works with a full matrix of Ny**2
C              elements (stored in the conventional way).  Computations
C              with banded matrices cost less time and storage than with
C              full matrices if  2*ML+MU .lt. Ny.  If you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide routine JAC to
C              compute the partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of JAC.
C
C              INFO(6) is also used for the user to input an external
C              linear solver, for example a sparse direct solver. 
C
C          **** Do you want to use your own linear solver?
C                 no - Do you want to solve the problem using a full (dense)
C                      matrix (and not a special banded structure) ...
C                       yes - set INFO(6) = 0
C                        no - set INFO(6) = 1
C                             and provide the lower (ML) and upper (MU)
C                             bandwidths by setting
C                             IWORK(1)=ML
C                             IWORK(2)=MU 
C                yes - set INFO(6) = 2 and provide routines JAC and PSOL.
C                      Jac is used to evaluate and factorize the Jacobian,
C                      and PSOL is used to solve the linear system with
C                      the factorized Jacobian.   ****
C
C       INFO(7) - You can specify a maximum (absolute value of)
C              stepsize, so that the code will avoid passing over very
C              large regions.
C
C          ****  Do you want the code to decide on its own the maximum
C                stepsize ...
C                 yes - set INFO(7) = 0
C                  no - set INFO(7) = 1
C                       and define HMAX by setting
C                       RWORK(2) = HMAX ****
C
C       INFO(8) -  Differential/algebraic systems may occasionally
C              suffer from severe scaling difficulties on the first
C              step.  If you know a great deal about the scaling of 
C              your problem, you can help to alleviate this problem 
C              by specifying an initial stepsize H0.
C
C          ****  Do you want the code to define its own initial
C                stepsize ...
C                 yes - set INFO(8) = 0
C                  no - set INFO(8) = 1
C                       and define H0 by setting
C                       RWORK(3) = H0 ****
C
C       INFO(9) -  If storage is a severe problem, you can save some
C              storage by restricting the maximum method order MAXORD.
C              The default value is 5.  For each order decrease below 5,
C              the code requires NEQ fewer locations, but it is likely 
C              to be slower.  In any case, you must have 
C              1 .le. MAXORD .le. 5.
C          ****  Do you want the maximum order to default to 5 ...
C                 yes - set INFO(9) = 0
C                  no - set INFO(9) = 1
C                       and define MAXORD by setting
C                       IWORK(3) = MAXORD ****
C
C       INFO(10) - If you know that certain components of the
C              solutions to your equations are always nonnegative
C              (or nonpositive), it may help to set this
C              parameter.  Three options are available:
C              1.  To have constraint checking in Y only in the initial
C                  condition calculation.
C              2.  To enforce nonnegativity in Y during the integration.
C              3.  To enforce both options 1 and 2.
C
C              When selecting option 2 or 3, it is probably best to try the
C              code without using this option first, and only use
C              this option if that does not work well.
C
C          ****  Do you want the code to solve the problem without
C                invoking any special inequality constraints ...
C                 yes - set INFO(10) = 0
C                  no - set INFO(10) = 1 to have option 1 enforced 
C                  no - set INFO(10) = 2 to have option 2 enforced
C                  no - set INFO(10) = 3 to have option 3 enforced ****
C
C                  If you have specified INFO(10) = 1 or 3, then you
C                  will also need to identify how each component of Y
C                  in the initial condition calculation is constrained.
C                  You must set:
C                  IWORK(40+I) = +1 if Y(I) must be .GE. 0,
C                  IWORK(40+I) = +2 if Y(I) must be .GT. 0,
C                  IWORK(40+I) = -1 if Y(I) must be .LE. 0, while
C                  IWORK(40+I) = -2 if Y(I) must be .LT. 0, while
C                  IWORK(40+I) =  0 if Y(I) is not constrained.
C
C        ***Notes: The constraint control acts only on the state variables, 
C                  excluding the quadrature/output and sensitivity variables.
C                  We recommend to set this parameter only when absolutely
C                  necessary (i.e., if the solution fails to satisfy those 
C                  constraints without setting the parameters).
C
C       INFO(11) - DDASPK normally requires the initial T, Y, and
C              YPRIME to be consistent.  That is, you must have
C              G(T,Y,YPRIME) = 0 at the initial T.  Note that this
C              includes the sensitivity equations, if applicable.
C              If you do not know the initial conditions precisely,
C              in some cases DDASPK may be able to compute them.
C
C              We will refer to those variables where their time 
C              derivatives appear explicitly in G as differential 
C              variables. The rest of the variables in Y 
C              are called algebraic variables. Note that algebraic 
C              variables can be further divided into index-1, index-2 
C              and index-3, etc.(see the references at the end of this
C              documentation)
C              
C              Differential equations are those equations that explicitly 
C              include time derivatives YPRIME. The rest of the equations 
C              are called algebraic constraints, which can be further 
C              divided into index-1, index-2 and index-3, etc. 
C
C              Denoting the differential variables in Y by Y_d, and the 
C              algebraic variables by Y_a, which includes Y_a1 (index-1 
C              algebraic variables) and Y_a2 (index-2 algebraic 
C              variables), the initialization is applicable for DAEs 
C              with index up to 2. It is also required that the mass 
C              matrix G_{Y'_d} have full row rank.
C
C              DDASPK has the following options for initialization: 
C
C              1.  Given Y_d, calculate Y_a and Y'_d. 
C                  This option is applicable for index-0 or index-1 
C                  problems where the mass matrix G_{Y'_d} is square.
C              2.  Given Y', calculate Y.
C                  This option is applicable for any index problem but
C                  could generate a singular Jacobian matrix if it is not 
C                  well-posed (i.e., if there is not a unique solution
C                  for the given problem).
C              3.  Given Y, Y', calculate a new pair (Y, Y').
C                  This option will take a small artificial step with
C                  the backward Euler method to find a consistent new
C                  pair (Y, Y'). It allows a subset of Y or Y' to be 
C                  fixed during the initialization. 
C                  This option is applicable for any index DAEs. 
C              4.  Option 4 is for the index-1 DAEs where the mass matrix 
C                  is not square or for index-2 DAEs. This option ensures
C                  that not only the original constraints but also the 
C                  derived constraints are satisfied at the initial values.
C                  It is similar to option 3 in fixing variables. But the 
C                  fixed variables MUST be differential variables in option 4.
C              5.  Option 5 is a special case for option 4. It is applicable
C                  when the original constraints are already satisfied, and 
C                  only the derived constraints need to be satisfied. Y_d, 
C                  and perhaps a subset of Y_a, are fixed. Y'_d and the rest 
C                  of Y_a are calculated. 
C 
C                  Option 5 works better than option 4 when you want to 
C                  fix all of the values of Y_d that satisfy the original
C                  constraints to calculate other unknown variables during 
C                  the initialization.
C
C              6.  Option 6 is a special case for option 3. When you are not 
C                  sure which variables (derivatives) should be fixed, you 
C                  can use this option to free all of them and do not need
C                  to tell DDASPK the variable properties in IWORK array.
C
C              In all cases, initial values for the given
C              components are input. Initial guesses for
C              the unknown components must also be provided as input.
C              
C          ****  Are the initial T, Y, YPRIME consistent? ...
C
C                 yes - set INFO(11) = 0
C                  no - set INFO(11) = I to invoke option I above.
C                       I is a integer from 1 to 6.
C
C                  If you have specified INFO(11) = 1, then you
C                  will also need to identify  which are the
C                  differential and which are the algebraic variables.  
C                  You must set:
C
C                  IWORK(LID+I) = 1, 2 or 3 if Y(I) is a differential 
C                                    variable;
C                  IWORK(LID+I) = -1,  if Y(I) is an index-1 algebraic 
C                                      variable;
C
C                  where LID = 40 if INFO(10) = 0 or 2 and LID = 40+NEQ
C                  if INFO(10) = 1 or 3.
C
C                  If you have specified INFO(11) = 3 or 4, then you must 
C                  identify which differential variables or their 
C                  derivatives should be fixed by setting:
C
C                  IWORK(LID+I) = 1 if Y(I) is a differential variable 
C                                    but not fixed;
C                  IWORK(LID+I) = 2 if Y(I) is a differential variable and 
C                                    fixed, but Y'(I) is free;
C                  IWORK(LID+I) = 3 if Y(I) is a differential variable and 
C                                    Y'(I) is fixed, but Y(I) is free;
C                  IWORK(LID+I) = -1 or -2,  if Y(I) is an algebraic variable;
C
C
C                  If you have specified INFO(11) = 4 or 5, you need to 
C                  specify which algebraic constraints need to be 
C                  differentiated to satisfy the derived constraints 
C                  when obtaining the values of the unknown components of 
C                  Y or YPRIME by setting:
C
C                  IWORK(LID+NY+I) = 1, if equation I needs to be 
C                                       differentiated;
C                  IWORK(LID+NY+I) = 0, otherwise.
C
C                  Usually, the equation that needs to be differentiated in an
C                  index-2 system is an index-2 constraint and/or an index-1 
C                  constraint that contains index-2 variables. For an index-1 
C                  DAE where the mass matrix is not square, all of the index-1 
C                  constraints need to be differentiated.
C
C                  You MUST scale by CJ those algebraic constraints that 
C                  require differentiation (multiply them by CJ) in the RES 
C                  routine. The constraints will be differentiated 
C                  automatically inside DDASPK by either the finite difference 
C                  method, or ADIFOR, or input by the user in RES 
C                  (with IRES=2), whichever is specified in INFO(20). 
C                  If INFO(20) > 2 (ADIFOR option), then you also need to 
C                  provide the ADIFOR-generated routine in T_RES (see the 
C                  documentation for T_RES for details).
C
C                  If you have specified INFO(11) = 4 or 5, you also need to 
C                  specify which algebraic variables should be fixed during
C                  the differentiation by setting:
C
C                  IWORK(LID+I) = -1, if Y(I) is an algebraic variable and 
C                                             fixed; 
C                  IWORK(LID+I) = -2, if Y(I) is an algebraic variable, but
C                                             not fixed.
C                  IWORK(LID+I) = 1, 2 or 3 if Y(I) is a differential variable;
C
C         **Note** If you have specified sensitivity analysis
C                  (INFO(19)>0), sensitivities of differential variables
C                  will automatically be classified as differential,
C                  and sensitivities of algebraic variables will be
C                  classified as algebraic. The sensitivity of a fixed variable
C                  is also classified as fixed (for the initial condition
C                  computation).  
C
C       INFO(12) - Except for the addition of the RES argument CJ,
C              DDASPK by default is downward-compatible with DDASSL,
C              which uses only direct (dense or banded) methods to solve 
C              the linear systems involved.  You must set INFO(12) to
C              indicate whether you want the direct methods or the
C              Krylov iterative method.
C          ****   Do you want DDASPK to use standard direct methods
C                 (dense or banded) or the Krylov (iterative) method ? ...
C                   direct methods - set INFO(12) = 0.
C                   Krylov method  - set INFO(12) = 1,
C                       and check the settings of INFO(13) and INFO(15).
C
C       INFO(13) - used when INFO(12) = 1 (Krylov methods).  
C              DDASPK uses scalars MAXL, KMP, NRMAX, and EPLI for the
C              iterative solution of linear systems.  INFO(13) allows 
C              you to override the default values of these parameters.  
C              These parameters and their defaults are as follows:
C              MAXL = maximum number of iterations in the SPIGMR 
C                 algorithm (MAXL .le. NEQ).  The default is 
C                 MAXL = MIN(5,NEQ).
C              KMP = number of vectors on which orthogonalization is 
C                 done in the SPIGMR algorithm.  The default is 
C                 KMP = MAXL, which corresponds to complete GMRES 
C                 iteration, as opposed to the incomplete form.  
C              NRMAX = maximum number of restarts of the SPIGMR 
C                 algorithm per nonlinear iteration.  The default is
C                 NRMAX = 5.
C              EPLI = convergence test constant in SPIGMR algorithm.
C                 The default is EPLI = 0.05.
C              Note that the length of RWORK depends on both MAXL 
C              and KMP.  See the definition of LRW below.
C          ****   Are MAXL, KMP, and EPLI to be given their
C                 default values? ...
C                  yes - set INFO(13) = 0
C                   no - set INFO(13) = 1,
C                        and set all of the following:
C                        IWORK(24) = MAXL (1 .le. MAXL .le. NEQ)
C                        IWORK(25) = KMP  (1 .le. KMP .le. MAXL)
C                        IWORK(26) = NRMAX  (NRMAX .ge. 0)
C                        RWORK(10) = EPLI (0 .lt. EPLI .lt. 1.0) ****
C
C        INFO(14) - used with INFO(11) > 0 (initial condition 
C               calculation is requested).  In this case, you may
C               request control to be returned to the calling program
C               immediately after the initial condition calculation,
C               before proceeding to the integration of the system
C               (e.g. to examine the computed Y and YPRIME).
C               If this is done, and if the initialization succeeded
C               (IDID = 4), you should reset INFO(11) to 0 for the
C               next call, to prevent the solver from repeating the 
C               initialization (and to avoid an infinite loop). 
C          ****   Do you want to proceed to the integration after
C                 the initial condition calculation is done? ...
C                 yes - set INFO(14) = 0
C                  no - set INFO(14) = 1                        ****
C 
C          **Note** When INFO(14) = 1 is chosen, the integration will
C               restart (INFO(1)=0) with the consistent initial conditions
C               on the next call of DDASPK. 
C
C        INFO(15) - used when INFO(12) = 1 (Krylov methods).
C               When using preconditioning in the Krylov method,
C               you must supply a routine, PSOL, which solves the
C               associated linear systems using P.
C               The usage of DDASPK is simpler if PSOL can carry out
C               the solution without any prior calculation of data.
C               However, if some partial derivative data is to be
C               calculated in advance and used repeatedly in PSOL,
C               then you must supply a JAC routine to do this,
C               and set INFO(15) to indicate that JAC is to be called
C               for this purpose.  For example, P might be an
C               approximation to a part of the matrix A which can be
C               calculated and LU-factored for repeated solutions of
C               the preconditioner system.  The arrays WP and IWP
C               (described under JAC and PSOL) can be used to
C               communicate data between JAC and PSOL.
C          ****   Does PSOL operate with no prior preparation ...
C                 yes - set INFO(15) = 0 (no JAC routine)
C                  no - set INFO(15) = 1
C                       and supply a JAC routine to evaluate and
C                       preprocess any required Jacobian data.  ****
C
C         INFO(16) - option to exclude algebraic variables from
C               the error test.  
C          ****   Do you wish to control errors locally on
C                 all the variables...
C                 yes - set INFO(16) = 0
C                  no - set INFO(16) = 1
C                       If you have specified INFO(16) = 1, then you
C                       will also need to identify  which are the
C                       differential and which are the algebraic
C                       components (algebraic components are components
C                       whose derivatives do not appear explicitly
C                       in the function G(T,Y,YPRIME)).  You must set:
C                       IWORK(LID+I) = 1, 2 or 3 if Y(I) is a differential 
C                                      variable, and
C                       IWORK(LID+I) = -1 if Y(I) is an algebraic variable,
C                       where LID = 40 if INFO(10) = 0 or 2 and 
C                       LID = 40 + NEQ if INFO(10) = 1 or 3.
C
C                       Note that if you have specified sensitivity
C                       analysis (INFO(19) > 0), the sensitivities
C                       of variables that have been marked as algebraic
C                       will also be excluded from the error test.  ****
C
C           *Note*: The quadrature/output variables will be automatically 
C                   included in the error control except INFO(18) is set
C                   differently.
C
C       INFO(17) - used when INFO(11) > 0 (DDASPK is to do an 
C              initial condition calculation).
C              DDASPK uses several heuristic control quantities in the
C              initial condition calculation.  They have default values,
C              but can  also be set by the user using INFO(17).
C              These parameters and their defaults are as follows:
C              MXNIT  = maximum number of Newton iterations
C                 per Jacobian or preconditioner evaluation.
C                 The default is:
C                 MXNIT =  5 in the direct case (INFO(12) = 0), and
C                 MXNIT = 15 in the Krylov case (INFO(12) = 1).
C              MXNJ   = maximum number of Jacobian or preconditioner
C                 evaluations.  The default is:
C                 MXNJ = 6 in the direct case (INFO(12) = 0), and
C                 MXNJ = 2 in the Krylov case (INFO(12) = 1).
C              MXNH   = maximum number of values of the artificial
C                 stepsize parameter H to be tried if INFO(11) = 1.
C                 The default is MXNH = 5.
C                 NOTE: the maximum number of Newton iterations
C                 allowed in all is MXNIT*MXNJ*MXNH if INFO(11) = 1,
C                 and MXNIT*MXNJ if INFO(11) = 2.
C              LSOFF  = flag to turn off the linesearch algorithm
C                 (LSOFF = 0 means linesearch is on, LSOFF = 1 means
C                 it is turned off).  The default is LSOFF = 0.
C              STPTOL = minimum scaled step in linesearch algorithm.
C                 The default is STPTOL = (unit roundoff)**(2/3).
C              EPINIT = swing factor in the Newton iteration convergence
C                 test.  The test is applied to the residual vector,
C                 premultiplied by the approximate Jacobian (in the
C                 direct case) or the preconditioner (in the Krylov
C                 case).  For convergence, the weighted RMS norm of
C                 this vector (scaled by the error weights) must be
C                 less than EPINIT*EPCON, where EPCON = .33 is the
C                 analogous test constant used in the time steps.
C                 The default is EPINIT = .01.
C          ****   Are the initial condition heuristic controls to be 
C                 given their default values...
C                  yes - set INFO(17) = 0
C                   no - set INFO(17) = 1,
C                        and set all of the following:
C                        IWORK(32) = MXNIT (.GT. 0)
C                        IWORK(33) = MXNJ (.GT. 0)
C                        IWORK(34) = MXNH (.GT. 0)
C                        IWORK(35) = LSOFF ( = 0 or 1)
C                        RWORK(14) = STPTOL (.GT. 0.0)
C                        RWORK(15) = EPINIT (.GT. 0.0)  ****
C
C        INFO(18) -- Error control option for quadrature/output variables.
C                    This option is only used when INFO(28) > 0. 
C
C        ****  Do you wish to exclude the quadrature/output variables
C              from the error control test?  
C              Yes - Set INFO(18) = 1.
C              No  - Set INFO(18) = 0.   **** 
C
C
C        INFO(19) --Sensitivity toggle.
C
C         ****   Do you wish to have a sensitivity analysis performed
C                on the given ODE/DAE?
C                Yes - Set INFO(19) = NP, where NP equals to the number of
C                      parameters involved in the system to be solved, 
C                      including parameters that appear only in the initial 
C                      conditions. The sensitivity parameters that appear 
C                      directly in RES are normally stored in SENPAR(*).
C                      The dimension of SENPAR(*) is stored in INFO(22), which
C                      must be specified by the user.
C                No  - Set INFO(19) = 0.  ****
C
C
C        INFO(20) --options for evaluating the sensitivity equations. 
C                 This option is used only when INFO(19)>0 or INFO(11)=5.
C                
C         ****   Do you wish to use a finite difference method (FDM)
C                for the approximation of the sensitivity analysis residuals?
C                Yes - In this case the routine RES should define the
C                      residual of the differential (or differential-
C                      algebraic) equations only, i.e. it should define Ny
C                      residuals. Set the value of INFO(20) to one of the
C                      following choices:
C                      INFO(20) = 0. A second order centered FDM is used.
C                      INFO(20) = 1. A first order forward FDM is used.
C
C                No  - There are three options:
C                   1. Set INFO(20) = 2 and provide a user-supplied RES. 
C                      In this case, the RES routine should also compute 
C                      the residuals of the sensitivity equations, which 
C                      means that it should define NEQ = Ny(Np + 1) residuals.
C                      IRES is used to determine whether to compute the 
C                      sensitivities or not in RES routine. (see the
C                      description of RES in this documentation)
C 
C                   2. Set INFO(20) = 3 or 4 and provide an ADIFOR-generated 
C                      routine in argument G_RES. 
C                      The residuals for the sensitivities will be evaluated 
C                      via ADIFOR. 
C
C                      INFO(20)=3 uses the seed matrix option in ADIFOR,
C                      INFO(20)=4 uses the matrix-vector product option in 
C                      ADIFOR. Usually INFO(20)=3 takes less computation time 
C                      but requires more work space than INFO(20)=4. When the 
C                      parallel option (INFO(27)>1) is selected, INFO(20)=3 
C                      cannot be chosen.
C
C                   3. Set INFO(5) = 5 and provide an ADIFOR-generated 
C                      routine in argument G_RES. This option is only for 
C                      the staggered corrector method (INFO(25)=1).
C                      The residuals for the sensitivities will be evaluated 
C                      by matrix times vector methods using the compressed
C                      sparse row format. G_RES should be a routine generated 
C                      by ADIFOR with the SparsLinC option and Y, YPRIME, 
C                      SENPAR as independent variables. The maximal percentage
C                      of the number of nonzero elements in the matrix 
C                      (dG/dy, dG/dy') should be set in RWORK(17).  
C                                                                ****
C
C       INFO(21) --Perturbation factor option. The selection in the case 
C              INFO(20) < 2 of the perturbation used in the finite 
C              difference approximation to the sensitivity equations is 
C              determined by the expression
C              delta_i = cnst*max(SENPAR(i), 1/vnorm_i),
C              where vnorm_i is the 2-norm of the vector v_i defined by 
C              v_i(j) = WT(j)/WT(i*Ny + j) for j = 1..Ny, 
C              where WT(*) is a vector of weights determined by RTOL, ATOL 
C              and Y, and cnst is the perturbation factor. The default value 
C              for cnst is 1.0d-3. This option allows the user to change this
C              value of cnst.
C
C       ****   Do you wish to alter the value of cnst?
C              No  - Set info(21) = 0.
C              Yes - Set info(21) = 1 and set RWORK(16) equal to the
C                    desired value.    ****
C
C        INFO(22) --Number of parameters that appear in RES. If INFO(20) is 
C              NOT equal to 2, INFO(22) must be equal to or greater than 
C              the dimension of SENPAR(*).
C                        
C        INFO(23) --Error control option for sensitivity variables.
C
C        ****  Do you wish to include the sensitivity variables
C              in the error control test?  Note that DDASPK  will
C              produce more reliable results for the sensitivities if they
C              are included in the error test.
C              Yes - Set INFO(23) = 0.
C              No  - Set INFO(23) = 1.   **** 
C
C         *** Note: all variables will be included in the Newton iteration
C             test.
C
C        INFO(24) --Sensitivity of a derived quantity.
C             In addition to computing the sensitivities of the solution 
C             and its derivative with respect to the parameters, you may 
C             want to compute the sensitivity of a quantity Q(t,y,y',p) 
C             with respect to the parameters.  To do this, you can call 
C             the auxiliary routine DSENSD
C                 CALL DSENSD(QRES, NEQ, T, Y, YPRIME, QSEN, INFO, 
C                *   RWORK, IWORK, RPAR, IPAR, SENPAR)
C             in between calls to DDASPK to find the sensitivity of the 
C             derived quantity Q at the outpout points of DDASPK.  
C             (For further details, see the documentation of DSENSD).
C
C         ****  Do you plan to be calling DSENSD to obtain sensitivities
C               of the derived quantity Q?
C               No  - Set INFO(24) = 0
C               Yes - Set INFO(24) = NQ,
C                     where NQ is the dimension of the vector function Q.
C                     Note that you must also augment the dimension
C                     of RWORK by at 2*NQ locations (see description
C                     of LRW).       ****
C
C       INFO(25) -- Method option for sensitivity analysis.
C
C         **** Do you plan to use the staggered method to obtain the 
C              sensitivities?
C              No  - Set INFO(25) = 0. The simultaneous corrector method
C                    will be used, where the DAE with the sensitivity 
C                    equations are solved simultaneously.
C              Yes - There are two options:
C                    Set INFO(25) = 1. The staggered corrector method will
C                        be used. On each time step, the DAEs are solved 
C                        first, then the sensitivity equations are solved 
C                        in each Newton iteration. 
C                    Set INFO(25) = 2. The staggered direct method will be
C                        used. On each time step, the DAEs are solved 
C                        first, then the linear system of sensitivity 
C                        equations are solved directly by one iteration.
C                        In this method, the Jacobian will be evaluated and
C                        factorized in every time step.   ****
C
C              ***Note***For the Krylov method, the staggered director and the 
C                        staggered corrector methods are the same.
C
C       INFO(26) -- This option is for parallel computing only. For serial 
C               computing, set INFO(26) = 0.
C
C               In parallel computing, each processor has an identity to
C               differentiate itself from the others. This identity is 
C               generated by a set-up routine in the parallel environment,
C               such as message passing interface. All of the identities 
C               form a consecutive integer array that begins with 0.
C
C       INFO(27) -- This option is for parallel computing only.
C               Set INFO(27) to 0 (DASPK will change it to 1) or to 1 for
C               serial computation. 
C
C               INFO(27) stores the total number of processors that participate
C               in the parallel computing. In parallel computing (INFO(27)>1),
C               INFO(27) is used to distribute the sensitivity equations
C               to each processor. The processor with ID i=INFO(26) will
C               compute the sensitivities with respect to parameter number
C               i+1, i+1+INFO(27), i+1+2*INFO(27), ..., i+1+j*INFO(27).
C                                      
C       INFO(28) -- Option for efficient computation of quadrature variables.
C
C               Does your problem include quadrature variables that should be
C               treated separately for efficiency?
C                No   -- Set INFO(28) = 0.
C                Yes  -- Set INFO(28) to the number of state variables to be 
C                        treated as quadrature variables and make sure that 
C                        the state variable vector is arranged so that the
C                        quadrature variables appear last.
C
C       INFO(29) and INFO(30) 
C                -- These two options are reserved for sensitivity 
C                   computation by the adjoint method, which is implemented
C                   in DDASPKADJOINT. Users who do not need the adjoint should
C                   set them to 0.
C
C   RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
C               error tolerances to tell the code how accurately you
C               want the solution to be computed.  They must be defined
C               as variables because the code may change them.
C               You have two choices --
C                     Both RTOL and ATOL are scalars (INFO(2) = 0), or
C                     both RTOL and ATOL are vectors (INFO(2) = 1) of length
C                     NEQ.
C               In either case all components must be non-negative.
C
C               The tolerances are used by the code in a local error
C               test at each step which requires roughly that
C                        abs(local error in Y(i)) .le. EWT(i) ,
C               where EWT(i) = RTOL*abs(Y(i)) + ATOL is an error weight 
C               quantity, for each vector component.
C               (More specifically, a maximum of root-mean-square norms over
C               the state variables and each of the sensitivity variables is 
C               used to measure the size of vectors, and the error test uses 
C               the magnitude of the solution at the beginning of the step.)
C
C               The true (global) error is the difference between the
C               true solution of the initial value problem and the
C               computed approximation.  Practically all present day
C               codes, including this one, control the local error at
C               each step and do not even attempt to control the global
C               error directly.
C
C               Usually, but not always, the true accuracy of
C               the computed Y is comparable to the error tolerances.
C               This code will usually, but not always, deliver a more
C               accurate solution if you reduce the tolerances and
C               integrate again.  By comparing two such solutions you 
C               can get a fairly reliable idea of the true error in the
C               solution at the larger tolerances.
C
C               Setting ATOL = 0 results in a pure relative error test
C               on that component.  Setting RTOL = 0 results in a pure
C               absolute error test on that component.  A mixed test
C               with non-zero RTOL and ATOL corresponds roughly to a
C               relative error test when the solution component is
C               much bigger than ATOL and to an absolute error test
C               when the solution component is smaller than the
C               threshold ATOL.
C
C               The code will not attempt to compute a solution at an
C               accuracy unreasonable for the machine being used.  It
C               will advise you if you ask for too much accuracy and
C               inform you as to the maximum accuracy it believes
C               possible.
C
C  RWORK(*) -- a real work array, which should be dimensioned in your
C               calling program with a length equal to the value of
C               LRW (or greater). If INFO(20) = 1 or 2 and INFO(21)=1, set
C               RWORK(16) = perturbation factor for finite difference.
C               If INFO(20) = 5 is chosen, set
C               RWORK(17) = maximum percentage of the number of nonzero
C                           elements in the matrix (dG/dy, dG'/dy).
C               RWORK(17) should be between 0 and 1.
C
C  LRW -- Set it to the declared length of the RWORK array.  The
C              minimum length depends on the options you have selected,
C              given by a base value plus additional storage as described
C              below. Some work space is used to pass the information 
C              between calls of DDASPK, while the rest, which we call ITMP,
C              is just temporally used by DDASPK and can be reused by other 
C              program between calls of DDASPK.
C
C              In the following, we denote NY = NEQ/(INFO(19)+1) to be the
C              number of state variables.
C              If INFO(12) = 0 (standard direct method), the base value is
C              base = 50 + max(MAXORD+4,7)*NEQ. The default value is 
C              MAXORD = 5 (see INFO(9)).  With the default MAXORD, 
C              base = 50 + 9*NEQ.
C              If INFO(12) = 1 (Krylov method), the base value is
C              base = 50 + max(MAXORD+5,9)*NEQ + (MAXL+3+MIN0(1,MAXL-KMP))*NY
C                      + (MAXL+3)*MAXL + 1 + LENWP.
C              See PSOL for description of LENWP.  The default values are:
C              MAXORD = 5 (see INFO(9)), MAXL = min(5,NEQ) and KMP = MAXL 
C              (see INFO(13)). 
C              With the default values for MAXORD, MAXL and KMP,
C              base = 91 + 10*NEQ + 8*NY + LENWP.
C
C              Additional storage must be added to the base value for
C              any or all of the following options:
C              If INFO(16) = 1, add NEQ.
C              If INFO(12) = 0 (standard direct method), and 
C                 if INFO(6) = 0 (dense matrix), add NY**2, 
C                    if INFO(5) = 3, ITMP1 = NY*(3*NY)
C                 if INFO(6) = 1 (banded matrix), then add (2*ML+MU+1)*NY.
C                    if INFO(5) = 0, ITMP1 = 2*(NY/(ML+MU+1)+1),
C                    if INFO(5) = 3, ITMP1 = (ML+MU+1)*(3*NY);
C              If INFO(12) = 1 (Krylov method), and
C                 if INFO(5) = 2, ITMP1 = NY;
C              add ITMP = max(ITMP1, ITMP2, ITMP3), where ITMP2 and ITMP3 are
C              defined by
C                 If INFO(11) = 4 (index-2 two-step process)
C                    if INFO(20) < 2 and INFO(19) = 0, ITMP2 = 2*NY.
C                    if INFO(20) < 2 and INFO(19) > 0, ITMP2 = 2*NEQ+4*NY+
C                                                              2*INFO(24).
C                    if INFO(20) > 2 and INFO(19) = 0, ITMP2 = NY.
C                    if INFO(20) > 2 and INFO(19) > 0, ITMP2 = 2*NY + INFO(22).
C                 Otherwise, ITMP2 = 0
C                 If INFO(19) > 0 (sensitivity analysis)
C                    if INFO(20) < 2, ITMP3 = 4*NY + 2*INFO(24).
C                    if INFO(20) = 3, ITMP3 = 
C                       INFO(19)*(2*NY+MAX(NY,INFO(24))+INFO(22)).
C                    if INFO(20) = 4, ITMP3 = INFO(22).
C                    if INFO(20) = 5, ITMP3 = 
C                       INFO(19)*(NY+1)+RWORK(17)*(2*NY*NY).
C                 Otherwise, ITMP3 = 0.
C
C  IWORK(*) -- An integer work array, which should be dimensioned in
C              your calling program with a length equal to the value
C              of LIW (or greater).
C
C  LIW -- Set it to the declared length of the IWORK array.  The
C             minimum length depends on the options you have selected,
C             given by a base value plus additional storage as described
C             below.
C
C             If INFO(12) = 0 (standard direct method), the base value is
C             base = 40 + NY.
C             If INFO(12) = 1 (Krylov method), the base value is
C             base = 40 + LENIWP (See PSOL for description of LENIWP).
C             
C             Additional storage must be added to the base value for
C             any or all of the following options:
C             IF INFO(10) = 1 or 3, add NY.
C             If INFO(11) = 1 or 3, or INFO(16) = 1, add NY.
C             If INFO(11) = 4 or 5, add 2*NY.
C             IF INFO(12) = 0 (direct method) and INFO(5) = 2, add NY.
C             IF INFO(19) > 0 and INFO(20) = 5, add 
C                NY + 1 + RWORK(17)*(2*NY*NY) + INFO(22).
C             IF (INFO(12) = 0 AND INFO(5) = 2) OR (INFO(20) = 5) 
C                (ADIFOR with SparsLinC options), add
C                3*NY + INFO(22)
C
C  RPAR, IPAR -- These are arrays of double precision and integer type,
C             respectively, which are available for you to use
C             for communication between your program that calls
C             DDASPK and the RES routine (and the JAC and PSOL
C             routines).  They are not altered by DDASPK.
C             If you do not need RPAR or IPAR, ignore these
C             parameters by treating them as dummy arguments.
C             If you do choose to use them, dimension them in
C             your calling program and in RES (and in JAC and PSOL)
C             as arrays of appropriate length.
C
C  JAC -- This is the name of a routine that you may supply
C         (optionally) that relates to the Jacobian matrix of the
C         nonlinear system that the code must solve at each time step.
C         The role of JAC (and its call sequence) depends on whether
C         a direct (INFO(12) = 0) or Krylov (INFO(12) = 1) method 
C         is selected, or whether you have opted to provide your own solver 
C         (INFO(6) = 2).
C
C         **** INFO(12) = 0 (direct methods):
C           If you are letting DDASPK generate partial derivatives
C           numerically (INFO(5) = 0), then JAC can be absent
C           (or perhaps a dummy routine to satisfy the loader).
C           If you choose INFO(5) = 1, then you must supply a 
C           JAC routine to compute the matrix A = dG/dY + CJ*dG/dYPRIME.  
C           It must have the form
C
C           SUBROUTINE JAC (T, Y, YPRIME, PD, CJ, RPAR, IPAR, SENPAR, IJAC)
C
C           The JAC routine must dimension Y, YPRIME, and PD (and RPAR
C           and IPAR if used).  CJ is a scalar which is input to JAC.
C           For the given values of T, Y, and YPRIME, the JAC routine
C           must evaluate the nonzero elements of the matrix A, and 
C           store these values in the array PD.  The elements of PD are 
C           set to zero before each call to JAC, so that only nonzero
C           elements need to be defined. IJAC is a flag used for 
C           initialization. Depending on the values of
C           IJAC, JAC should be defined as:
C             If IJAC = 0, JAC computes A = dG/dY + CJ*dG/dYPRIME 
C                          as usual. 
C             If IJAC = 1, JAC computes A = dG/dY + CJ*dG/dYPRIME
C                          where some variables or derivatives are fixed.
C                          If Y(I) is fixed, dG/dY(i)=0;
C                          If YPRIME(I) is fixed, dG/dYPRIME(I)=0.
C             If IJAC = 2, JAC computes A = dG/dY + CJ*DG/dYPRIME for
C                          the second stage of Hessenberg index-2 
C                          initialization, where the differential variables 
C                          are fixed in the differential equations. 
C           The way you store the elements into the PD array depends
C           on the structure of the matrix indicated by INFO(6).
C           *** INFO(6) = 0 (full or dense matrix) ***
C               Give PD a first dimension of NEQ.  When you evaluate the
C               nonzero partial derivatives of equation i (i.e. of G(i))
C               with respect to component j (of Y and YPRIME), you must
C               store the element in PD according to
C                  PD(i,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
C           *** INFO(6) = 1 (banded matrix with half-bandwidths ML, MU
C                            as described under INFO(6)) ***
C               Give PD a first dimension of 2*ML+MU+1.  When you 
C               evaluate the nonzero partial derivatives of equation i 
C               (i.e. of G(i)) with respect to component j (of Y and 
C               YPRIME), you must store the element in PD according to 
C                  IROW = i - j + ML + MU + 1
C                  PD(IROW,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
C
C           If INFO(5) = 2, JAC is the name of an ADIFOR-generated routine
C              with the SparsLinC option. 
C           If INFO(5) = 3, JAC is the name of an ADIFOR-generated routine
C              with the seed matrix option.
C   
C           A step by step procedure to generate a routine for JAC via ADIFOR:
C            1. Put all the codes related to
C               SUBROUTINE RES (
C              *      T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C               in a file called "res.f"
C            2. Create a file "res.cmp" with one line:
C                res.f
C            3. Create a file "resjac.adf" with the following lines if 
C               the SparsLinC option is selected,
C
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime
C                 AD_OVARS = delta
C                 AD_PREFIX = j
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_FLAVOR = sparse
C
C               or with the following lines if the seed matrix option is 
C               selected,
C
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime
C                 AD_OVARS = delta
C                 AD_PREFIX = j
C                 AD_SUPPRESS_LDG = true
C                 AD_PMAX = min(NY, MU+ML+1)
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_FLAVOR = dense
C
C            4. Run Adifor to generate J_RES with the command:
C               % Adifor AD_SCRIPT=resjac.adf
C
C               If the SparseLinC option is selected, the JAC must have 
C               the form
C
C               SUBROUTINE JAC (
C              *      T, Y, J_Y, YPRIME, J_YPRIME, CJ, DELTA, J_DELTA, 
C              *      IRES, RPAR, IPAR, SENPAR)
C
C               where J_Y, J_YPRIME, J_DELTA are the seed matrices for the
C               SparseLinC option. All of the variables except DELTA and 
C               J_DELTA are input. DELTA and J_DELTA are output variables.
C
C               If the seed matrix option is selected, the JAC must have 
C               the form
C
C               SUBROUTINE JAC (
C              *      NP, T, Y, J_Y, YPRIME, J_YPRIME, CJ, DELTA, J_DELTA,
C              *      IRES, RPAR, IPAR, SENPAR)
C
C               where J_Y, J_YPRIME, J_DELTA are the seed matrices. All 
C               of the variables except DELTA and J_DELTA are input. DELTA 
C               and J_DELTA are output variables.
C
C          **** INFO(12) = 1 (Krylov method):
C            If you are not calculating Jacobian data in advance for use
C            in PSOL (INFO(15) = 0), JAC can be absent (or perhaps a
C            dummy routine to satisfy the loader).  Otherwise, you may
C            supply a JAC routine to compute and preprocess any parts of
C            of the Jacobian matrix  A = dG/dY + CJ*dG/dYPRIME that are
C            involved in the preconditioner matrix P.
C            It is to have the form
C
C            SUBROUTINE JAC (RES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
C                            WK, H, CJ, WP, IWP, IER, RPAR, IPAR, SENPAR)
C
C           Where NEQ is the number of state equations excluding those of
C           quadrature form you have specified in INFO(28).
C           The JAC routine must dimension Y, YPRIME, REWT, SAVR, WK,
C           and (if used) WP, IWP, RPAR, IPAR and SENPAR.
C           The Y, YPRIME, and SAVR arrays contain the current values
C           of Y, YPRIME, and the residual G, respectively.  
C           The array WK is work space of length NEQ.  
C           H is the step size.  CJ is a scalar, input to JAC, that is
C           normally proportional to 1/H.  REWT is an array of 
C           reciprocal error weights, 1/EWT(i), where EWT(i) is
C           RTOL*abs(Y(i)) + ATOL (unless you supplied routine DDAWTS
C           instead), for use in JAC if needed.  For example, if JAC
C           computes difference quotient approximations to partial
C           derivatives, the REWT array may be useful in setting the
C           increments used.  The JAC routine should do any
C           factorization operations called for, in preparation for
C           solving linear systems in PSOL.  The matrix P should
C           be an approximation to the Jacobian,
C           A = dG/dY + CJ*dG/dYPRIME.
C
C           WP and IWP are real and integer work arrays which you may
C           use for communication between your JAC routine and your
C           PSOL routine.  These may be used to store elements of the 
C           preconditioner P, or related matrix data (such as factored
C           forms).  They are not altered by DDASPK.
C           If you do not need WP or IWP, ignore these parameters by
C           treating them as dummy arguments.  If you do use them,
C           dimension them appropriately in your JAC and PSOL routines.
C           See the PSOL description for instructions on setting 
C           the lengths of WP and IWP.
C
C           On return, JAC should set the error flag IER as follows..
C             IER = 0    if JAC was successful,
C             IER .ne. 0 if JAC was unsuccessful (e.g. if Y or YPRIME
C                        was illegal, or a singular matrix is found).
C           (If IER .ne. 0, a smaller stepsize will be tried.)
C           IER = 0 on entry to JAC, so it needs to be reset only on a 
C           failure.
C           If RES is used within JAC, then a nonzero value of IRES will
C           override any nonzero value of IER (see the RES description).
C
C         Regardless of the method type, routine JAC must not
C         alter T, Y(*), YPRIME(*), H, CJ, or REWT(*).
C         You must declare the name JAC in an EXTERNAL statement in
C         your program that calls DDASPK.
C
C       **** INFO(6) = 2 (provide your own linear solver):
C           If you opt to provide your own linear solver, you
C           must input JAC routine to compute the matrix 
C                  A = dG/dY + CJ*dG/dYPRIME
C           and factorize it. It must have the form
C
C            SUBROUTINE JAC (RES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
C                            WK, H, CJ, WP, IWP, IER, RPAR, IPAR, SENPAR,
C                            ICOPT, ID)
C
C           Where NEQ is the number of state equations excluding those of
C           quadrature form that you have specified in INFO(28).
C           The JAC routine must dimension Y, YPRIME, REWT, SAVR, WK,
C           and (if used) WP, IWP, RPAR, IPAR and SENPAR.
C           The Y, YPRIME, and SAVR arrays contain the current values
C           of Y, YPRIME, and the residual G, respectively.  
C           The array WK is work space of length NEQ.  
C           H is the step size.  CJ is a scalar, input to JAC, that is
C           normally proportional to 1/H.  REWT is an array of 
C           reciprocal error weights, 1/EWT(i), where EWT(i) is
C           RTOL*abs(Y(i)) + ATOL (unless you supplied routine DDAWTS
C           instead), for use in JAC if needed. The JAC routine should do 
C           any factorization operations called for, in preparation for
C           solving linear systems in PSOL.  
C
C           WP and IWP are real and integer work arrays which you may
C           use for communication between your JAC routine and your
C           PSOL routine.  These may be used to store elements of the 
C           Jacobian A, or related matrix data (such as factored
C           forms).  They are not altered by DDASPK.
C           If you do not need WP or IWP, ignore these parameters by
C           treating them as dummy arguments.  If you do use them,
C           dimension them appropriately in your JAC and PSOL routines.
C           See the PSOL description for instructions on setting 
C           the lengths of WP and IWP.
C
C           ICOPT is used for the Jacobian evaluation during the 
C           initialization. The Jacobian may be different from A during the 
C           initialization because some variables or derivatives may be
C           fixed. DDASPK works better if you evaluate the Jacobian differently
C           for the initializaton. ID is an integer work array to indicate 
C           which variables or derivatives are fixed (see the documentation of
C           INFO(11) for details). 
C    
C           IER is a error return flag as follows..
C             IER = 0    if JAC was successful,
C             IER .ne. 0 if JAC was unsuccessful (e.g. if Y or YPRIME
C                        was illegal, or a singular matrix is found).
C
C           If you have specified INFO(6) = 2, then INFO(5) will be ignored.
C
C  PSOL --  This is the name of a routine you must supply if you have
C         selected to provide your own linear solver (INFO(6)) or selected  
C         the Krylov method (INFO(12) = 1) with preconditioning.
C         In the direct case (INFO(12) = 0), PSOL can be absent 
C         (a dummy routine may have to be supplied to satisfy the 
C         loader) if INFO(6)<2.  
C
C         When supplied with INFO(12) = 0 and INFO(6) = 2, the PSOL 
C         routine is to have the form
C
C         SUBROUTINE PSOL (NEQ, WP, IWP, B, IER, RPAR, IPAR)
C
C         The PSOL routine must solve linear systems of the form 
C         A*x = b where  A = dG/dY + CJ*dG/dYPRIME.
C
C         When supplied with INFO(12) = 1, the PSOL routine is to 
C         have the form
C
C         SUBROUTINE PSOL (NEQ, T, Y, YPRIME, SAVR, WK, CJ, WGHT,
C                          WP, IWP, B, EPLIN, IER, RPAR, IPAR, SENPAR)
C
C         The PSOL routine must solve linear systems of the form 
C         P*x = b where P is the left preconditioner matrix.
C
C         The right-hand side vector b is in the B array on input, and
C         PSOL must return the solution vector x in B.
C         The Y, YPRIME, and SAVR arrays contain the current values
C         of Y, YPRIME, and the residual G, respectively.  
C
C         Work space required by JAC and/or PSOL, and space for data to
C         be communicated from JAC to PSOL is made available in the form
C         of arrays WP and IWP, which are parts of the RWORK and IWORK
C         arrays, respectively.  The lengths of these real and integer
C         work spaces WP and IWP must be supplied in LENWP and LENIWP,
C         respectively, as follows..
C           IWORK(27) = LENWP = length of real work space WP
C           IWORK(28) = LENIWP = length of integer work space IWP.
C
C         WK is a work array of length NEQ for use by PSOL.
C         CJ is a scalar, input to PSOL, that is normally proportional
C         to 1/H (H = stepsize).  If the old value of CJ
C         (at the time of the last JAC call) is needed, it must have
C         been saved by JAC in WP.
C
C         WGHT is an array of weights, to be used if PSOL uses an
C         iterative method and performs a convergence test.  (In terms
C         of the argument REWT to JAC, WGHT is REWT/sqrt(NEQ).)
C         If PSOL uses an iterative method, it should use EPLIN
C         (a heuristic parameter) as the bound on the weighted norm of
C         the residual for the computed solution.  Specifically, the
C         residual vector R should satisfy
C              SQRT (SUM ( (R(i)*WGHT(i))**2 ) ) .le. EPLIN
C
C         PSOL must not alter NEQ, T, Y, YPRIME, SAVR, CJ, WGHT, EPLIN.
C
C         On return, PSOL should set the error flag IER as follows..
C           IER = 0 if PSOL was successful,
C           IER .lt. 0 if an unrecoverable error occurred, meaning
C                 control will be passed to the calling routine,
C           IER .gt. 0 if a recoverable error occurred, meaning that
C                 the step will be retried with the same step size
C                 but with a call to JAC to update necessary data,
C                 unless the Jacobian data is current, in which case
C                 the step will be retried with a smaller step size.
C           IER = 0 on entry to PSOL so need be reset only on a failure.
C
C         You must declare the name PSOL in an EXTERNAL statement in
C         your program that calls DDASPK.
C
C
C  SENPAR -- This is an array of sensitivity parameters that appear in 
C         RES. They are not altered by DDASPK. 
C         If you are not doing sensitivity analysis, ignore it by 
C         treating it as a dummy argument. You must use it and dimension 
C         it in your calling program and in RES (and in JAC and PSOL) 
C         as an array of appropriate length if the sensitivity parameters 
C         appear in the RES routine. 
C         The sensitivity parameters related to RES must be stored in 
C         the first part of SENPAR if sensitivity equations are not 
C         input by the user (i.e., INFO(20) .ne. 2). 
C
C  G_RES -- Evaluation of sensitivity residuals by automatic differentiation.
C         This is the name of an ADIFOR-generated routine you must supply 
C         if you have chosen INFO(19) > 0 and INFO(20) > 2. See the 
C         description of INFO(20) for more details. 
C 
C         The procedures to generate a routine for G_RES via ADIFOR are 
C         very similar to those for JAC except that the script file is 
C         different. The contents of the script file used for G_RES is as 
C         follows:
C             1. For seed matrix option (INFO(20)=3):
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime, senpar
C                 AD_OVARS = delta
C                 AD_PREFIX = g
C                 AD_PMAX = INFO(19)  (number of sensitivity parameters)
C                 AD_SUPPRESS_LDG = true
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_FLAVOR = dense
C                
C              The G_RES must have the form
C
C              SUBROUTINE G_RES(NP,T,Y,G_Y,YPRIME,G_YPRIME,
C            *                  CJ,DELTA,G_DELTA,
C            *                  IRES,RPAR,IPAR,SENPAR,G_SENPAR)
C
C              where NP is the total number of parameters, G_Y, G_YPRIME, 
C              G_DELTA and G_SENPAR are the seed matrices. All of the 
C              variables except DELTA and G_DELTA are input. DELTA and 
C              G_DELTA are output variables. It calculates
C                DELTA   = G(t, y, y', senpar)
C                G_DELTA = G_y*G_y + G_y'*G_YPRIME + G_senpar*G_SENPAR
C
C             2. For matrix-vector product only option (INFO(20)=4):
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime, senpar
C                 AD_OVARS = delta
C                 AD_PREFIX = g
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_SCALAR_GRADIENTS = true
C
C              The G_RES must have the form
C
C              SUBROUTINE G_RES(T,Y,G_Y,YPRIME,G_YPRIME,CJ,DELTA,G_DELTA,
C            *                  IRES,RPAR,IPAR,SENPAR,G_SENPAR)
C
C              where G_Y, G_YPRIME and G_SENPAR are the seed vectors.
C              All of the variables except DELTA and G_DELTA are input. DELTA 
C              and G_DELTA are output variables. It calculates
C                DELTA   = G(t, y, y', senpar)
C                G_DELTA = G_y*G_y + G_y'*G_YPRIME + G_senpar*G_SENPAR
C            
C             3. For matrix times vector options (INFO(20)=5):
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime, senpar
C                 AD_OVARS = delta
C                 AD_FLAVOR = sparse
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_PREFIX = s
C                                    
C              The G_RES must have the form
C
C              SUBROUTINE G_RES(T,Y,S_Y,YPRIME,S_YPRIME,CJ,DELTA,S_DELTA,
C            *                  IRES,RPAR,IPAR,SENPAR,S_SENPAR)
C
C              where S_Y, S_YPRIME and S_SENPAR are the seed matrices.
C              All of the variables except DELTA and S_DELTA are input. DELTA 
C              and S_DELTA are output variables. It calculates
C                DELTA   = G(t, y, y', senpar)
C                G_DELTA = (G_y*S_y, G_y'*S_YPRIME, G_senpar*S_SENPAR)
C
C  K_RES -- Evaluation of matrix-vector product for Krylov iteration by 
C           analytic or autoamatic differentiation.
C           This is the name of a routine you must supply if you have 
C           chosen INFO(5) > 0 (analytic or ADIFOR evaluation of partial
C           derivatives in matrix-vector product) and INFO(12) = 1 (Krylov
C           option).
C
C           Depending on the value of INFO(5), K_RES is given by:
C      **** 1. INFO(5) = 1. The matrix-vector product in the Krylov linear
C              iteration is computed by a user-input routine, which has the
C              form
C
C            SUBROUTINE K_RES(T,Y,YPRIME,CJ,IRES,RPAR,IPAR,SENPAR,V,AV)
C
C              where V is the input vector and AV is the matrix-vector
C              product AV = (G_y + CJ* G_y') * V.
C
C           2. INFO(5) = 2. K_RES will be the name of an ADIFOR-generated
C              routine, with the matrix-vector product only option. The
C              script used to generate K_RES is as follows
C
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime
C                 AD_OVARS = delta
C                 AD_PREFIX = g
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_SCALAR_GRADIENTS = true            
C
C              If K_RES is generated by other methods, it must have the form
C
C              SUBROUTINE K_RES(T,Y,G_Y,YPRIME,G_YPRIME,CJ,DELTA,G_DELTA,
C            *                  IRES,RPAR,IPAR,SENPAR)
C
C              where G_Y and G_YPRIME are the seed vectors.
C              All of the variables except DELTA and G_DELTA are input.  
C             DELTA  and G_DELTA are output variables. It calculates
C                DELTA   = G(t, y, y')
C                G_DELTA = G_y*G_y + G_y'*G_YPRIME
C
C           3. INFO(5) = 3, a modified Newton will be used to solve the 
C              nonlinear system at each time step. The matrix will be the 
C              old Jacobian stored in WP and IWP. The matrix-vector product 
C              is calculated by K_RES, which has the form
C
C              SUBROUTINE K_RES(NEQ, WP, IWP, V, AV)
C
C              where V is the input vector and AV is the matrix-vector
C              product AV = J * V.              
C              
C  T_RES -- Evaluate the time derivatives of the index-2 constraints by
C           automatic differentiation.
C           This is the name of a routine you must supply if you choose
C           INFO(11) = 4 or 5 and INFO(20) > 2.
C           
C           Depending on the value of INFO(19), T_RES is given by:
C      **** 1. INFO(19) = 0. No sensitivity is considered. The T_RES can 
C              be generated by the following ADIFOR script
C
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, t
C                 AD_OVARS = delta
C                 AD_PREFIX = g
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_SCALAR_GRADIENTS = true            
C
C              Finally, T_RES must have the form
C
C              SUBROUTINE T_RES(T,G_T,Y,G_Y,YPRIME,CJ,DELTA,G_DELTA,
C            *                  IRES,RPAR,IPAR,SENPAR)
C
C              where G_T and G_Y are the seed vectors. DELTA and G_DELTA
C              are the output variables, the rest are the input variables.
C              It calculates
C                DELTA   = G(t, y, y')
C                G_DELTA = G_t*G_T + G_y*G_Y
C              
C           2. INFO(19) > 0. Sensitivity is considered. T_RES can be
C              generated by the following two ADIFOR scripts.
C
C              First generate G_RES by
C
C                 AD_PROG = res.cmp
C                 AD_TOP = res
C                 AD_IVARS = y, yprime, senpar
C                 AD_OVARS = delta
C                 AD_PREFIX = g
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_SCALAR_GRADIENTS = true
C
C              Then input G_RES to ADIFOR and generate T_G_RES by
C       
C                 AD_PROG = g_res.cmp
C                 AD_TOP = g_res
C                 AD_IVARS = y, g_y, t
C                 AD_OVARS = delta, g_delta
C                 AD_PREFIX = t
C                 AD_EXCEPTION_FLAVOR=performance 
C                 AD_SCALAR_GRADIENTS = true
C
C              Finally, rename T_G_RES to T_RES.       
C
C              Finally, T_RES must have the form 
C
C              SUBROUTINE T_RES(T,G_T,Y,G__Y,G_Y,G_G_Y,YPRIME,G_YPRIME, 
C            *                  CJ,DELTA,G__DELTA, G_DELTA, G_G_DELTA,
C            *                  IRES,RPAR,IPAR,SENPAR,G_SENPAR)      
C 
C              where G_T, G_Y, G__Y and G_G_Y are the seed vectors. 
C              DELTA, G_DELTA, G__DELTA and G_G_DELTA are the output 
C              variables, and the rest are the input variables.
C              It calculates
C                DELTA   = G(t, y, y', senpar)
C                G_DELTA = G_y*G_Y + G_y'*G_YPRIME + G_senpar*G_SENPAR
C                G__DELTA = G_t*G_T + G_y*G_Y
C                G_G_DELTA = d(G_DELTA)/dt
C              
C  A_RES --  This is the name of a routine which is reserved only for
C            DASPKADJOINT. Users who do not need adjoint do not have to
C            worry about this. You can simply omit it in the argument 
C            list or declare it as dummy real variable when calling DDASPK. 
C                   
C
C   ****Note*****
C         If the Adifor3.0 (or newer version ADIFOR) is used, the script
C         file, "spec.ad3", can be written as follows:
C
C         1. Adifor3.0 with seed matrix option
C
C                   Function
C                     JM
C                   MaxColsM
C                     number equal to INFO(19)
C                   Independent
C                     y
C                     yprime
C                     senpar
C                   Dependent
C                     delta
C                   Top
C                     res
C                   ExceptionDefaultMode
C                     Performance
C                   Cache
C                     JM_Cache
C                   OutputDir
C                     JM_Output  
C
C        2. Adifor3.0 with matrix-vector product option
C
C                   Function
C                     JM
C                   MaxColsM
C                     1
C                   Independent
C                     y
C                     yprime
C                     senpar
C                   Dependent
C                     delta
C                   Top
C                     res
C                   ScalarDerivatives
C                     true
C                   ExceptionDefaultMode
C                     Performance
C                   Cache
C                     JM_Cache
C                   OutputDir
C                     JM_Output  
C
C        3. Adifor3.0 with SparseLinC option
C                 
C                   Function
C                     SJ
C                   MaxIndependent
C                     2
C                   Independent
C                     y
C                     yprime
C                   Dependent
C                     delta
C                   Top
C                     res
C                   ExceptionDefaultMode
C                     Performance
C                   Cache
C                     SJ_Cache
C                   OutputDir
C                     SJ_Output  
C
C         4. To run Adifor3.0, type:
C               % Adifor3.0 -spec spec.ad3 res.f
C         
C        Note that the independent variables are different for 
C        G_RES, K_RES, J_RES, and T_RES.
C
C  OPTIONALLY REPLACEABLE SUBROUTINE:
C
C  DDASPK uses a weighted root-mean-square norm to measure the 
C  size of various error vectors.  The weights used in this norm
C  are set in the following routine:
C
C    SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, EWT, RPAR, IPAR)
C    DIMENSION RTOL(*), ATOL(*), Y(*), EWT(*), RPAR(*), IPAR(*)
C
C  A DDAWTS routine has been included with DDASPK which sets the
C  weights according to
C    EWT(I) = RTOL*ABS(Y(I)) + ATOL
C  in the case of scalar tolerances (IWT = 0) or
C    EWT(I) = RTOL(I)*ABS(Y(I)) + ATOL(I)
C  in the case of array tolerances (IWT = 1).  (IWT is INFO(2).)
C  In some special cases, it may be appropriate for you to define
C  your own error weights by writing a routine DDAWTS to be 
C  called instead of the version supplied.  However, this should 
C  be attempted only after careful thought and consideration. 
C  If you supply this routine, you may use the tolerances and Y 
C  as appropriate, but do not overwrite these variables.  You
C  may also use RPAR and IPAR to communicate data as appropriate.
C  ***Note: Aside from the values of the weights, the choice of 
C  norm used in DDASPK (weighted root-mean-square) is not subject
C  to replacement by the user.  In this respect, DDASPK is not
C  downward-compatible with the original DDASSL solver (in which
C  the norm routine was optionally user-replaceable).
C
C
C------OUTPUT - AFTER ANY RETURN FROM DDASPK----------------------------
C
C  The principal aim of the code is to return a computed solution at
C  T = TOUT, although it is also possible to obtain intermediate
C  results along the way.  To find out whether the code achieved its
C  goal or if the integration process was interrupted before the task
C  was completed, you must check the IDID parameter.
C
C
C   T -- The output value of T is the point to which the solution
C        was successfully advanced.
C
C   Y(*) --contains the computed solution (and sensitivity) approximation 
C          at T.
C
C   YPRIME(*) -- contains the computed derivative (and sensitivity 
C        derivative) approximation at T.
C
C   IDID -- reports what the code did, described as follows:
C
C                     *** TASK COMPLETED ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- a step was successfully taken in the
C                   intermediate-output mode.  The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- the integration to TSTOP was successfully
C                   completed (T = TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- the integration to TOUT was successfully
C                   completed (T = TOUT) by stepping past TOUT.
C                   Y(*) and YPRIME(*) are obtained by interpolation.
C
C           IDID = 4 -- the initial condition calculation, with
C                   INFO(11) > 0, was successful, and INFO(14) = 1.
C                   No integration steps were taken, and the solution
C                   is not considered to have been started.
C
C                    *** TASK INTERRUPTED ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- a large amount of work has been expended
C                     (about 500 steps).
C
C           IDID = -2 -- the error tolerances are too stringent.
C
C           IDID = -3 -- the local error test cannot be satisfied
C                     because you specified a zero component in ATOL
C                     and the corresponding computed solution component
C                     is zero.  Thus, a pure relative error test is
C                     impossible for this component.
C
C           IDID = -5 -- there were repeated failures in the evaluation
C                     or processing of the preconditioner (in JAC).
C
C           IDID = -6 -- DDASPK had repeated error test failures on the
C                     last attempted step.
C
C           IDID = -7 -- the nonlinear system solver in the time 
C                     integration could not converge.
C
C           IDID = -8 -- the matrix of partial derivatives appears
C                     to be singular (direct method).
C
C           IDID = -9 -- the nonlinear system solver in the time integration
C                     failed to achieve convergence, and there were repeated 
C                     error test failures in this step.
C
C           IDID =-10 -- the nonlinear system solver in the time integration 
C                     failed to achieve convergence because IRES was equal 
C                     to -1.
C
C           IDID =-11 -- IRES = -2 was encountered and control is
C                     being returned to the calling program.
C
C           IDID =-12 -- DDASPK failed to compute the initial Y, YPRIME.
C
C           IDID =-13 -- unrecoverable error encountered inside user's
C                     PSOL routine, and control is being returned to
C                     the calling program.
C
C           IDID =-14 -- the Krylov linear system solver could not 
C                     achieve convergence.
C
C           IDID =-15,..,-32 -- Not applicable for this code.
C
C                    *** TASK TERMINATED ***
C                reported by the value of IDID=-33
C
C           IDID = -33 -- the code has encountered trouble from which
C                   it cannot recover.  A message is output
C                   explaining the trouble and control is returned
C                   to the calling program.  For example, this occurs
C                   when invalid input is detected.
C
C   RTOL, ATOL -- these quantities remain unchanged except when
C               IDID = -2.  In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration.  However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C   RWORK, IWORK -- contain information which is usually of no interest
C               to the user but is necessary for subsequent calls. 
C               However, you may be interested in the performance data
C               listed below.  These quantities are accessed in RWORK 
C               and IWORK but have internal mnemonic names, as follows..
C
C               RWORK(3)--contains H, the step size h to be attempted
C                        on the next step.
C
C               RWORK(4)--contains TN, the current value of the
C                        independent variable, i.e. the farthest point
C                        integration has reached.  This will differ 
C                        from T if interpolation has been performed 
C                        (IDID = 3).
C
C               RWORK(7)--contains HOLD, the stepsize used on the last
C                        successful step.  If INFO(11) = INFO(14) = 1,
C                        this contains the value of H used in the
C                        initial condition calculation.
C
C               IWORK(7)--contains K, the order of the method to be 
C                        attempted on the next step.
C
C               IWORK(8)--contains KOLD, the order of the method used
C                        on the last step.
C
C               IWORK(11)--contains NST, the number of steps (in T) 
C                        taken so far.
C
C               IWORK(12)--contains NRE, the number of calls to RES 
C                        so far.  For the Direct solver, it also contains
C                        the number of calls to RES in the Jacobian 
C                        evaluations. In the case of sensitivity analysis,
C                        it counts only the number of evaluations
C                        of the state variable system.
C
C               IWORK(13)--contains NJE, the number of calls to JAC so
C                        far (Jacobian or preconditioner evaluations).
C
C               IWORK(14)--contains NETF, the total number of error test
C                        failures so far.
C
C               IWORK(15)--contains NCFN, the total number of nonlinear
C                        convergence failures so far (includes counts
C                        of singular iteration matrix or singular
C                        preconditioners).
C
C               IWORK(16)--contains NCFL, the number of convergence
C                        failures of the linear iteration so far.
C
C               IWORK(17)--contains LENIW, the length of IWORK actually
C                        required.  This is defined on normal returns 
C                        and on an illegal input return for
C                        insufficient storage.
C
C               IWORK(18)--contains LENRW, the length of RWORK actually
C                        required.  This is defined on normal returns 
C                        and on an illegal input return for
C                        insufficient storage.
C
C               IWORK(19)--contains NNI, the total number of nonlinear
C                        iterations so far (each of which calls a
C                        linear solver).
C
C               IWORK(20)--contains NLI, the total number of linear
C                        (Krylov) iterations for the state variables 
C                        so far.
C
C               IWORK(21)--contains NPS, the number of the PSOL calls so
C                        far, for preconditioning solve operations or
C                        for solutions with the user-supplied method.
C
C               IWORK(22)--contains NSE, the number of evaluations
C                        of the sensitivity system. For the staggered 
C                        corrector method, IWORK(22) is also the number
C                        of nonlinear iterations for the sensitivity
C                        variables.
C
C               IWORK(38)--number of linear (Krylov) iterations for the
C                        sensitivity variables
C
C               Note: The various counters in IWORK do not include 
C               counts during a call made with INFO(11) > 0 and
C               INFO(14) = 1.
C
C
C------INPUT - WHAT TO DO TO CONTINUE THE INTEGRATION  -----------------
C              (CALLS AFTER THE FIRST)
C
C     This code is organized so that subsequent calls to continue the
C     integration involve little (if any) additional effort on your
C     part.  You must monitor the IDID parameter in order to determine
C     what to do next.
C
C     Recalling that the principal task of the code is to integrate
C     from T to TOUT (the interval mode), usually all you will need
C     to do is specify a new TOUT upon reaching the current TOUT.
C
C     Do not alter any quantity not specifically permitted below.  In
C     particular do not alter NEQ, T, Y(*), YPRIME(*), RWORK(*), 
C     IWORK(*), or the differential equation in routine RES.  Any 
C     such alteration constitutes a new problem and must be treated 
C     as such, i.e. you must start afresh.
C
C     You cannot change from array to scalar error control or vice
C     versa (INFO(2)), but you can change the size of the entries of
C     RTOL or ATOL.  Increasing a tolerance makes the equation easier
C     to integrate.  Decreasing a tolerance will make the equation
C     harder to integrate and should generally be avoided.
C
C     You can switch from the intermediate-output mode to the
C     interval mode (INFO(3)) or vice versa at any time.
C
C     If it has been necessary to prevent the integration from going
C     past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C     code will not integrate to any TOUT beyond the currently
C     specified TSTOP.  Once TSTOP has been reached, you must change
C     the value of TSTOP or set INFO(4) = 0.  You may change INFO(4)
C     or TSTOP at any time but you must supply the value of TSTOP in
C     RWORK(1) whenever you set INFO(4) = 1.
C
C     Do not change INFO(5), INFO(6), INFO(12-17) or their associated
C     IWORK/RWORK locations unless you are going to restart the code.
C
C                    *** FOLLOWING A COMPLETED TASK ***
C
C     If..
C     IDID = 1, call the code again to continue the integration
C                  another step in the direction of TOUT.
C
C     IDID = 2 or 3, define a new TOUT and call the code again.
C                  TOUT must be different from T.  You cannot change
C                  the direction of integration without restarting.
C
C     IDID = 4, reset INFO(11) = 0 and call the code again to begin
C                  the integration.  (If you leave INFO(11) > 0 and
C                  INFO(14) = 1, you may generate an infinite loop.)
C                  In this situation, the next call to DDASPK is 
C                  considered to be the first call for the problem,
C                  in that all initializations have been completed.
C
C                    *** FOLLOWING AN INTERRUPTED TASK ***
C
C     To show the code that you realize the task was interrupted and
C     that you want to continue, you must take appropriate action and
C     set INFO(1) = 1.
C
C     If..
C     IDID = -1, the code has taken about 500 steps.  If you want to
C                  continue, set INFO(1) = 1 and call the code again.
C                  An additional 500 steps will be allowed.
C
C
C     IDID = -2, the error tolerances RTOL, ATOL have been increased
C                  to values the code estimates appropriate for
C                  continuing.  You may want to change them yourself.
C                  If you are sure you want to continue with relaxed
C                  error tolerances, set INFO(1) = 1 and call the code
C                  again.
C
C     IDID = -3, a solution component is zero and you set the
C                  corresponding component of ATOL to zero.  If you
C                  are sure you want to continue, you must first alter
C                  the error criterion to use positive values of ATOL 
C                  for those components corresponding to zero solution
C                  components, then set INFO(1) = 1 and call the code
C                  again.
C
C     IDID = -4  --- cannot occur with this code.
C
C     IDID = -5, your JAC routine failed with the Krylov method.  Check
C                  for errors in JAC and restart the integration.
C
C     IDID = -6, repeated error test failures occurred on the last
C                  attempted step in DDASPK.  A singularity in the
C                  solution may be present.  If you are absolutely
C                  certain you want to continue, you should restart
C                  the integration.  (Provide initial values of Y and
C                  YPRIME which are consistent.)
C
C     IDID = -7, repeated convergence test failures occurred on the last
C                  attempted step in DDASPK.  An inaccurate or ill-
C                  conditioned Jacobian or preconditioner may be the
C                  problem.  If you are absolutely certain you want
C                  to continue, you should restart the integration.
C
C
C     IDID = -8, the matrix of partial derivatives is singular, with
C                  the use of direct methods.  Some of your equations
C                  may be redundant.  DDASPK cannot solve the problem
C                  as stated.  It is possible that the redundant
C                  equations could be removed, and then DDASPK could
C                  solve the problem.  It is also possible that a
C                  solution to your problem either does not exist
C                  or is not unique.
C
C     IDID = -9, DDASPK had multiple convergence test failures, preceded
C                  by multiple error test failures, on the last
C                  attempted step.  It is possible that your problem is
C                  ill-posed and cannot be solved using this code.  Or,
C                  there may be a discontinuity or a singularity in the
C                  solution.  If you are absolutely certain you want to
C                  continue, you should restart the integration.
C
C     IDID = -10, DDASPK had multiple convergence test failures
C                  because IRES was equal to -1.  If you are
C                  absolutely certain you want to continue, you
C                  should restart the integration.
C
C     IDID = -11, there was an unrecoverable error (IRES = -2) from RES
C                  inside the nonlinear system solver.  Determine the
C                  cause before trying again.
C
C     IDID = -12, DDASPK failed to compute the initial Y and YPRIME
C                  vectors.  This could happen because the initial 
C                  approximation to Y or YPRIME was not very good, or
C                  because no consistent values of these vectors exist.
C                  The problem could also be caused by an inaccurate or
C                  singular iteration matrix, or a poor preconditioner.
C
C     IDID = -13, there was an unrecoverable error encountered inside 
C                  your PSOL routine.  Determine the cause before 
C                  trying again.
C
C     IDID = -14, the Krylov linear system solver failed to achieve
C                  convergence.  This may be due to ill-conditioning
C                  in the iteration matrix, or a singularity in the
C                  preconditioner (if one is being used).
C                  Another possibility is that there is a better
C                  choice of Krylov parameters (see INFO(13)).
C                  Possibly the failure is caused by redundant equations
C                  in the system, or by inconsistent equations.
C                  In that case, reformulate the system to make it
C                  consistent and non-redundant.
C
C     IDID = -15,..,-32 --- Cannot occur with this code.
C
C                       *** FOLLOWING A TERMINATED TASK ***
C
C     If IDID = -33, you cannot continue the solution of this problem.
C                  An attempt to do so will result in your run being
C                  terminated.
C
C  ---------------------------------------------------------------------
C
C***REFERENCES
C  1.  L. R. Petzold, A Description of DASSL: A Differential/Algebraic
C      System Solver, in Scientific Computing, R. S. Stepleman et al.
C      (Eds.), North-Holland, Amsterdam, 1983, pp. 65-68.
C  2.  K. E. Brenan, S. L. Campbell, and L. R. Petzold, Numerical 
C      Solution of Initial-Value Problems in Differential-Algebraic
C      Equations, Elsevier, New York, 1989 (second edition, SIAM 1996).
C  3.  P. N. Brown and A. C. Hindmarsh, Reduced Storage Matrix Methods
C      in Stiff ODE Systems, J. Applied Mathematics and Computation,
C      31 (1989), pp. 40-91.
C  4.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov
C      Methods in the Solution of Large-Scale Differential-Algebraic
C      Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.
C  5.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent
C      Initial Condition Calculation for Differential-Algebraic
C      Systems, SIAM J. Sci. Comp., 19 (1998), pp. 1495-1512.      
C  6.  T. Maly and L. R. Petzold, Numerical Methods and Software for
C      Sensitivity Analysis of Differential-Algebraic Systems, Applied
C      Numerical Mathematics 20 (1996), pp. 57-79.
C  7.  W. F. Feehery, J. E. Tolsma and P. I. Barton, Efficient 
C      sensitivity analysis of large-scale differential-algebraic 
C      systems,  Applied Numerical Mathematics 25 (1997), pp. 41-54.
C  8.  Shengtai Li and Linda R. Petzold, Design of new DASPK for 
C      sensitivity analysis, Technical Report TRC99-28, 
C      Department of Computer Science,
C      University of California, Santa Barbara, (1999).
C  9.  Shengtai Li and Linda R. Petzold, Software and Algorithms for 
C      Sensitivity Analysis of Large-Scale Differential Algebraic 
C      Systems, J. Comput and Appl. Math.  (2000), pp. 131-145.
C  10. Additional information is available via our website
C       http://www.engineering.ucsb.edu/~cse 
C
C***ROUTINES CALLED
C
C   The following are all the subordinate routines used by DDASPK.
C
C   DDASIC computes consistent initial conditions.
C   DYYPNW updates Y and YPRIME in linesearch for initial condition
C          calculation.
C   DDSTP  carries out one step of the integration.
C   DCNSTR/DCNST0 check the current solution for constraint violations.
C   DDAWTS sets error weight quantities.
C   DINVWT tests and inverts the error weights.
C   DDATRP performs interpolation to get an output solution.
C   DDWNRM computes the weighted root-mean-square norm of a vector.
C   D1MACH provides the unit roundoff of the computer.
C   XERRWD/XSETF/XSETUN/IXSAV is a package to handle error messages. 
C   DDASID nonlinear equation driver to initialize Y and YPRIME using
C          direct linear system solver methods.  Interfaces to Newton
C          solver (direct case).
C   DNSID  solves the nonlinear system for unknown initial values by
C          modified Newton iteration and direct linear system methods.
C   DMATID assembles the iteration matrix (direct case) for the 
C          initialization.
C   DLINSD carries out linesearch algorithm for initial condition
C          calculation (direct case).
C   DFNRMD calculates weighted norm of preconditioned residual in
C          initial condition calculation (direct case).
C   DNEDD  nonlinear equation driver for direct linear system solver
C          methods.  Interfaces to Newton solver (direct case).
C   DMATD  assembles the iteration matrix (direct case).
C   DNSD   solves the associated nonlinear system by modified
C          Newton iteration and direct linear system methods.
C   DSLVD  interfaces to linear system solver (direct case).
C   DDASIK nonlinear equation driver to initialize Y and YPRIME using
C          Krylov iterative linear system methods.  Interfaces to
C          Newton solver (Krylov case).
C   DNSIK  solves the nonlinear system for unknown initial values by
C          Newton iteration and Krylov iterative linear system methods.
C   DLINSK carries out linesearch algorithm for initial condition
C          calculation (Krylov case).
C   DFNRMK calculates weighted norm of preconditioned residual in
C          initial condition calculation (Krylov case).
C   DNEDK  nonlinear equation driver for iterative linear system solver
C          methods.  Interfaces to Newton solver (Krylov case).
C   DNSK   solves the associated nonlinear system by Inexact Newton
C          iteration and (linear) Krylov iteration.
C   DSLVK  interfaces to linear system solver (Krylov case).
C   DSPIGM solves a linear system by SPIGMR algorithm.
C   DATV   computes matrix-vector product in Krylov algorithm.
C   DORTH  performs orthogonalization of Krylov basis vectors.
C   DHEQR  performs QR factorization of Hessenberg matrix.
C   DHELS  finds least-squares solution of Hessenberg linear system.
C   DGEFA, DGESL, DGBFA, DGBSL are LINPACK routines for solving 
C          linear systems (dense or band direct methods).
C   DAXPY, DCOPY, DDOT, DNRM2, DSCAL are Basic Linear Algebra (BLAS)
C          routines.
C   DDSEN  computes the approximation to the sensitivity equations.
C   RMATV     evaluation of the sensitivity equations by matrix-vector 
C             product method.
C   RESIDX2   residual evaluation for the second stage initialization of
C             an index-2 system
C   RI2ADFMV  residual evaluation for the second stage initialization of
C             an index-2 system by ADIFOR with matrix-vector product only
C             option.
C   JDADFSP   Jacobian evaluation routine by ADIFOR with SparsLinC 
C             option.
C   JIDADFSP  Jacobian evaluation routine by ADIFOR with SparsLinC 
C             option for the initialization.
C   JIDADFSM  Jacobian evaluation routine by ADIFOR with seed matrix 
C             option for the initialization.
C   JDADFSM   Jacobian evaluation routine by ADIFOR with seed matrix 
C             option.
C   KMVBYADF  calculates the matrix-vector product in Krylov algorithm 
C             by ADIFOR with matrix-vector product only option.
C   RADFMV    residual evaluation routine by ADIFOR with matrix-vector
C             product only option.
C   RADFSM    residual evaluation routine by ADIFOR with seed matrix 
C             option.
C   JRADFSP   matrix evaluation routine by ADIFOR with SparsLinC option.
C   AMUVSP, AMUV1SP, AMUV2SP are routines to calculate the matrix-vector 
C             product in sparse format.
C   SD_JAC_ONLY  Jacobian evaluation for the stagger direct method.
C   ADFSP_INIT   initialization routine for the ADIFOR with SparsLinC 
C                option.
C   MPI_STEPSIZE the same stepsize control for all the processors in MPI
C                parallel environment.
C   DDRESAD      residual evaluation for the adjoint sensitivity system.
C   RESADIDX2    residual evaluation for the second stage initialization 
C                of the adjoint sensitivity system.
C   DCHI         cubic Hermite interpolation routine.
C
C The routines called directly by DDASPK are:
C   DCNST0, DDAWTS, DINVWT, D1MACH, DDWNRM, DDASIC, DDATRP, DDSTP,
C   XERRWD, MPI_STEPSIZE
C
C***END PROLOGUE DDASPK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL DONE, LAVL, LCFN, LCFL, LWARN
      DIMENSION Y(*),YPRIME(*)
      DIMENSION INFO(*)
      DIMENSION RWORK(LRW),IWORK(LIW)
      DIMENSION RTOL(*),ATOL(*)
      DIMENSION RPAR(*),IPAR(*),SENPAR(*)
      CHARACTER MSG*80
      EXTERNAL  RES, JAC, PSOL, DDASID, DDASIK, DNEDD, DNEDK
      EXTERNAL  G_RES, K_RES, T_RES, A_RES
C
C     Set pointers into IWORK.
C
      PARAMETER (LML=1, LMU=2, LMTYPE=4, 
     *   LIWM=1, LMXORD=3, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *   LNS=9, LNSTL=10, LNST=11, LNRE=12, LNJE=13, LETF=14, LNCFN=15,
     *   LNCFL=16, LNIW=17, LNRW=18, LNNI=19, LNLI=20, LNPS=21,
     *   LNSE=22, LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26, LLNWP=27,
     *   LLNIWP=28, LLOCWP=29, LLCIWP=30, LKPRIN=31, LMXNIT=32,
     *   LMXNJ=33, LMXNH=34, LLSOFF=35, LNPD=36, LNY=37, 
     *   LNLIS=38, LICNS=41)
C
C     Set pointers into RWORK.
C
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4, LCJ=5, LCJOLD=6,
     *   LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12,
     *   LEPCON=13, LSTOL=14, LEPIN=15, LPRT=16, LNZMX=17, 
     *   LTFN = 19,
     *   LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51)
C
      INTEGER LID, LENID, NONNEG, NPF, NPB, LASTADFSP
      SAVE NONNEG, NPF, NPB, LASTADFSP
C
C
C***FIRST EXECUTABLE STATEMENT  DDASPK
C
      MYID = INFO(26)
      IF (INFO(27) .EQ. 0) INFO(27) = 1 ! for the multiprocessing
      NUMPROCS = INFO(27)
      INFO19 = INFO(19)    ! Saved value of info(19)
      INFO24 = INFO(24)
C
      IF(INFO(1).NE.0) GO TO 100
C
C-----------------------------------------------------------------------
C     This block is executed for the initial call only.
C     It contains checking of inputs and initializations.
C-----------------------------------------------------------------------
C
C     First check INFO array to make sure all elements of INFO
C     Are within the proper range.  (INFO(1) is checked later, because
C     it must be tested on every call.) ITEMP holds the location
C     within INFO which may be out of range.
C
      DO 10 I=2,4
         ITEMP = I
         IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
 10   CONTINUE
      ITEMP = 5
      IF (INFO(5) .LT. 0. OR. INFO(5). GT. 3) GO TO 701
      ITEMP = 6
      IF (INFO(6) .LT. 0. OR. INFO(6). GT. 2) GO TO 701
      DO I=7,9
         ITEMP = I
         IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
      END DO
      ITEMP = 10
      IF(INFO(10).LT.0 .OR. INFO(10).GT.3) GO TO 701
      ITEMP = 11
      IF(INFO(11).LT.0 .OR. INFO(11).GT.6) GO TO 701
      DO 15 I=12,17
         ITEMP = I
         IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
 15      CONTINUE
      ITEMP = 18
      IF(INFO(18).LT.0 .OR. INFO(18).GT.2) GO TO 701
C     Check sensitivity inputs
      ITEMP = 19
      IF(INFO(19).LT.0)GO TO 701
C     IF(INFO(19).EQ.0)GO TO 17
      ITEMP = 20
      IF(INFO(20).LT.0 .OR. INFO(20).GT.5) GO TO 701
      ITEMP = 21
      IF(INFO(21).LT.0 .OR. INFO(21).GT.1) GO TO 701
      ITEMP = 22
C      IF(INFO(22).LT.0 .OR. INFO(22).GT.INFO(19)) GO TO 701
      IF(INFO(22).LT.0) GO TO 701
      ITEMP = 23
      IF(INFO(23).LT.0 .OR. INFO(23).GT.1) GO TO 701
      ITEMP = 24
C      IF(INFO(24) .LT. 0) GO TO 701
      ITEMP = 25
      IF (INFO(25) .LT. 0. OR. INFO(25) .GT. 2) GOTO 701
      IF (INFO(25) .EQ. 0 .AND. INFO(20) .EQ. 5) GOTO 730
      ITEMP = 26 
      IF (INFO(26) .LT. 0) GOTO 701 ! Processor ID
      ITEMP = 27 
      IF (INFO(27) .LT. 0) GOTO 701 ! number of processors
      IF (INFO(27) .GT. 1 .AND. INFO(20) .EQ. 3) GOTO 731
      ITEMP = 28 
      IF (INFO(28) .LT. 0 .OR. INFO(28) .GT. NEQ) GOTO 701
      ITEMP = 29 
      IF (INFO(29) .LT. 0 .OR. INFO(29) .GT. 5) GOTO 701
      ITEMP = 30 
      IF (INFO(30) .LT. 0 .OR. INFO(30) .GT. 7) GOTO 701
 17   CONTINUE
C
C     Check NEQ to see if it is positive.
C
      IF (NEQ .LE. 0) GO TO 702
C
C     Make sure that if we are computing sensitivities, NEQ is
C     at least equal to the number of problem parameters + 1.
C
C     Forward sensitivity integration
      IF (INFO(24) .GE. 0) THEN
         NY = NEQ / (INFO(19) + 1)
         IF (NY .LT. 1) GOTO 728
         IWORK(LNY) = NY
C.....Compute the number of equations for this processor
C
         NAVG = INFO(19)/NUMPROCS
         IF (MYID .LT. MOD(INFO(19), NUMPROCS)) THEN
            NPF = NAVG + 1
         ELSE
            NPF = NAVG
         END IF
         MYNEQ = IWORK(LNY)*(NPF+1)
         TFINAL = TOUT
      ELSE
C     Backward adjoint sensitivity integration
         NY = NEQ / (-INFO(24))
         IF (NY .LT. 1) GOTO 728
         IWORK(LNY) = NY
C.....Compute the number of equations for this processor
C
         NAVG = (-INFO(24))/NUMPROCS
         IF (MYID .LT. MOD(-INFO(24), NUMPROCS)) THEN
            NPB = NAVG + 1
         ELSE
            NPB = NAVG
         END IF
         MYNEQ = IWORK(LNY)*NPB
         NPF = 0
         TFINAL = RWORK(LTFN)
      END IF
      ITEMP =28
      IF (INFO(28) .GT. NY) GOTO 701
C     Check the size of RPAR(*) 
C
C     Check and compute maximum order.
C
      MXORD=5
      IF (INFO(9) .NE. 0) THEN
         MXORD=IWORK(LMXORD)
         IF (MXORD .LT. 1 .OR. MXORD .GT. 5) GO TO 703
         ENDIF
      IWORK(LMXORD)=MXORD
C
C     Set and/or check inputs for constraint checking (INFO(10) .NE. 0).
C     Set values for ICNFLG, NONNEG, and pointer LID.
C
      ICNFLG = 0
      NONNEG = 0
      LID = LICNS
      IF (INFO(10) .EQ. 0) GO TO 20
      IF (INFO(10) .EQ. 1) THEN
         ICNFLG = 1
         NONNEG = 0
         LID = LICNS + NY
      ELSEIF (INFO(10) .EQ. 2) THEN
         ICNFLG = 0
         NONNEG = 1
      ELSE
         ICNFLG = 1
         NONNEG = 1
         LID = LICNS + NY
      ENDIF
C
 20   CONTINUE
C
C     Set and/or check inputs for Krylov solver (INFO(12) .NE. 0).
C     If indicated, set default values for MAXL, KMP, NRMAX, and EPLI.
C     Otherwise, verify inputs required for iterative solver.
C
      IF (INFO(12) .EQ. 0) GO TO 25
C
      IWORK(LMITER) = INFO(12)
      IF (INFO(13) .EQ. 0) THEN
         IWORK(LMAXL) = MIN(5,MYNEQ)
         IWORK(LKMP) = IWORK(LMAXL)
         IWORK(LNRMAX) = 5
         RWORK(LEPLI) = 0.05D0
      ELSE
         IF(IWORK(LMAXL) .LT. 1 .OR. IWORK(LMAXL) .GT. MYNEQ)GOTO 720
         IF(IWORK(LKMP) .LT. 1 .OR. IWORK(LKMP) .GT. IWORK(LMAXL))
     1      GO TO 721
         IF(IWORK(LNRMAX) .LT. 0) GO TO 722
         IF(RWORK(LEPLI).LE.0.0D0 .OR. RWORK(LEPLI).GE.1.0D0)GO TO 723
         ENDIF
C
 25   CONTINUE
C
C     Set and/or check controls for the initial condition calculation
C     (INFO(11) .GT. 0).  If indicated, set default values.
C     Otherwise, verify inputs required for iterative solver.
C
      IF (INFO(11) .EQ. 0) GO TO 30
      IF (INFO(17) .EQ. 0) THEN
        IWORK(LMXNIT) = 5
        IF (INFO(12) .GT. 0) IWORK(LMXNIT) = 15
        IWORK(LMXNJ) = 6
        IF (INFO(12) .GT. 0) IWORK(LMXNJ) = 2
        IWORK(LMXNH) = 5
        IWORK(LLSOFF) = 0
        RWORK(LEPIN) = 0.01D0
      ELSE
        IF (IWORK(LMXNIT) .LE. 0) GO TO 725
        IF (IWORK(LMXNJ) .LE. 0) GO TO 725
        IF (IWORK(LMXNH) .LE. 0) GO TO 725
        LSOFF = IWORK(LLSOFF)
        IF (LSOFF .LT. 0 .OR. LSOFF .GT. 1) GO TO 725
        IF (RWORK(LEPIN) .LE. 0.0D0) GO TO 725
      ENDIF
C
 30   CONTINUE
C
C     Below is the computation and checking of the work array lengths
C     LENIW and LENRW, using direct methods (INFO(12) = 0) or
C     the Krylov methods (INFO(12) = 1).
C
      LENIC = 0
      IF (INFO(10) .EQ. 1 .OR. INFO(10) .EQ. 3) LENIC = NY
      LENID = 0
      IF (INFO(11).EQ.1 .OR. INFO(11).EQ.3 .OR. INFO(16).EQ.1)LENID = NY
      IF (INFO(11) .EQ.4 .OR. INFO(11).EQ.5) LENID = 2*NY
      ITMP1 = 0
      IF (INFO(12) .EQ. 0) THEN
C
C        Compute MTYPE, etc.  Check ML and MU.
C
         IF(INFO(6).EQ.0) THEN 
            LENPD = NY**2
            IF(INFO(5).EQ.0) THEN
               IWORK(LMTYPE)=2
            ELSE IF (INFO(5) .EQ. 1) THEN
               IWORK(LMTYPE)=1
            ELSE IF (INFO(5) .EQ. 2) THEN
               IWORK(LMTYPE)=3
            ELSE
               IWORK(LMTYPE)=4
               ITMP1 = 3*NY*NY
            ENDIF
         ELSE IF (INFO(6).EQ.1) THEN
            IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NY)GO TO 717
            IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NY)GO TO 718
            LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NY
            IF(INFO(5).EQ.0) THEN
               IWORK(LMTYPE)=6
               MBAND=IWORK(LML)+IWORK(LMU)+1
               MSAVE=(NY/MBAND)+1
               ITMP1 = 2*MSAVE
            ELSE IF (INFO(5) .EQ. 1) THEN
               IWORK(LMTYPE)=5
            ELSE IF (INFO(5) .EQ. 2) THEN
               IWORK(LMTYPE)=7
            ELSE
               IWORK(LMTYPE)=8
               MBAND = IWORK(LML)+IWORK(LMU)+1
               ITMP1 = MIN(MBAND,NY)*(3*NY)
            ENDIF
         ELSE 
            IWORK(LMTYPE) = 10
            LENPD = IWORK(LLNWP)            
         ENDIF
         IF (INFO(11) .EQ. 0) THEN
            LENRW = 50 + (IWORK(LMXORD)+4)*MYNEQ + LENPD
         ELSE
            LENRW = 50 + MAX0(IWORK(LMXORD)+4,7)*MYNEQ + LENPD
         END IF
C
C        Compute LENIW, LENWP, LENIWP.
C
         LENWP = 0
         IF (INFO(5) .EQ. 2) THEN
C
C        Temporary space for ADIFOR with SparseLinC option + Pivoting
            LENIWP = 2*NY    
         ELSE
C
C        Pivoting saving for LU solver
            LENIWP = NY
         END IF
         IF (INFO(6) .EQ. 2) THEN
            LENIWP = IWORK(LLNIWP)
            LENWP  = IWORK(LLNWP)
         END IF
         IWORK(LLNIWP) = LENIWP
         LENIW = 40 + LENIC + LENID + LENIWP
      ELSE IF (INFO(12) .EQ. 1)  THEN
         IWORK(LMTYPE) = INFO(5)
         IF (INFO(5) .EQ. 2) THEN
C
C     ADIFOR calculates the matrix-vector product in Krylov method
            ITMP1 = NY
         END IF
         MAXL = IWORK(LMAXL)
         LENWP = IWORK(LLNWP)
         LENIWP = IWORK(LLNIWP)
         LENPD = (MAXL+3+MIN0(1,MAXL-IWORK(LKMP)))*NY
     1         + (MAXL+3)*MAXL + 1 + LENWP
         IF (INFO(11) .EQ. 0) THEN
            LENRW = 50 + (IWORK(LMXORD)+5)*MYNEQ + LENPD
         ELSE
            LENRW = 50 + MAX0(IWORK(LMXORD)+5,9)*MYNEQ + LENPD
         END IF
         LENIW = 40 + LENIC + LENID + LENIWP
      ENDIF
      LISE = LENIW+1
      IF (NPF .GT. 0 .AND. INFO(20) .EQ. 5) THEN
C      iwork space for the Sparse Compressed Row format 
         LENIW = LISE + NY + 1 + INFO(22) + RWORK(LNZMX)*(2*NY*NY)
      END IF
C       
      IF(INFO(16) .NE. 0 .OR. 
     *     (INFO(18).EQ.1.AND.INFO(28).GT.0)) LENRW = LENRW + MYNEQ
      ITMP2 = 0
      IF (INFO(11) .EQ.4 .OR. INFO(11).EQ.5) THEN
C     Work space for the index-2 initialization
         IF (INFO(20) .LT. 2) THEN
            ITMP2 = 2*NY
            IF (NPF .GT. 0) ITMP2 = 2*MYNEQ + 4*NY + 2*INFO(24)
         ELSE
            ITMP2 = NY
            IF (NPF .GT. 0) ITMP2 = 2*NY+INFO(22)+3
         END IF
      END IF
      ITMP3 = 0
      IF(NPF.NE.0) THEN
C     Work space for the sensitivity evaluations
         IF (INFO(20) .LT. 2) THEN
            ITMP3 = 4*NY + 2*INFO(24)
         ELSE IF (INFO(20) .EQ. 3) THEN
            ITMP3 = NPF*(2*NY+MAX(INFO(24),NY)+INFO(22))
         ELSE IF(INFO(20) .EQ. 4) THEN
            ITMP3 = NPF  
         ELSE IF(INFO(20) .EQ. 5) THEN
C     Space for    (DF/DP), (DF/DY,DF/DY')
            ITMP3 = NPF*NY + RWORK(LNZMX)*(2*NY*NY)+NPF 
         END IF
      END IF
      ITMP = MAX0(ITMP1, ITMP2, ITMP3)
      LENRW = LENRW + ITMP
C     Integer work space for seed matrix of ADIFOR with sparseLinC 
C     option.
C
      IF ((INFO(12).EQ.0.AND.INFO(5).EQ.2) .OR. 
     *     (NPF.GT.0 .AND. INFO(20).EQ.5)) 
     *     LENIW = LENIW + 3*NY
C
C     For sensitivity parameters
      IF (INFO(20) .EQ. 5) LENIW = LENIW + INFO(22) 
C     Backward adjoint method
      IF (INFO(24) .LT. 0) THEN
C     The 2*NY space is for the seeded matrix, 3*NY space is used to 
C     store the information of state varibles and residual information.
C
         LENRW = LENRW + 6*NY
         IF (INFO(11) .GT. 0) THEN
C
C...  Require G_Y information in ADRES evaluation during initialization
            LENRW = LENRW + MYNEQ
         END IF
         IF (INFO(30) .EQ. 2 .OR. INFO(11) .EQ. 5) LENRW = LENRW + NY
         LENIW = LENIW + 2
C
C...  Index-2 needed more space for storage
C
C...  Store the forward ID(*) information
         IF (INFO(25) .EQ. 1 .OR. INFO(11) .GT. 0) LENIW = LENIW + 2*NY
      END IF
C
C...  Update the work space for the adjoint case of F_y'(t,y)
      IF (INFO(24) .LT. 0 .AND. (INFO(30).EQ.2 .OR. INFO(30).EQ.6)) THEN
         LENRW  = LENRW + (IWORK(LMXORD)+1)*MYNEQ + 3*MYNEQ
      END IF
C
C     Check lengths of RWORK and IWORK.
C
      IWORK(LNIW)=LENIW
      IWORK(LNRW)=LENRW
      IWORK(LNPD)=LENPD
      IWORK(LLOCWP) = LENPD-LENWP+1
      IWORK(LLCIWP) = LID + LENID
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     Check ICNSTR for legality.
C
      IF (LENIC .GT. 0) THEN
         DO 40 I = 1,NY
            ICI = IWORK(LICNS-1+I)
            IF (ICI .LT. -2 .OR. ICI .GT. 2) GO TO 726
 40      CONTINUE
      ENDIF
C
C     Check state variables for consistency with constraints.
C
      IF (LENIC .GT. 0) THEN
        CALL DCNST0(NY,Y,IWORK(LICNS),IRET)
        IF (IRET .NE. 0) GO TO 727
      ENDIF
C
C     Check ID for legality.
C
      IF (LENID .GT. 0) THEN
         DO 50 I = 1,NY
            IDI = IWORK(LID-1+I)
            IF (IDI.LT.-2 .OR. IDI.GT.3 .OR. IDI.EQ.0) GO TO 724
 50      CONTINUE
         IF ((INFO(11).EQ.4 .OR. INFO(11).EQ.5) .AND. 
     *        INFO(24) .GE. 0) THEN
            DO I = 1,NY
               IDI = IWORK(LID-1+I+NY)
               IF (IDI.LT.0 .OR. IDI.GT.1) GO TO 724
            END DO
         END IF
      ENDIF
C.....Redistribute the initial conditions for Y and YPRIME.
C
      IF (NUMPROCS .GT. 1) THEN
         IF (INFO(24) .LT. 0) THEN
            DO I = 1, NPB
               II = (I-1)*NY
               IPOS = ((I-1)*NUMPROCS + MYID)*NY
               DO J = 1, NY
                  Y(II + J ) = Y(IPOS + J)
                  YPRIME(II + J ) = YPRIME(IPOS + J)
               END DO
            END DO
         ELSE
            DO I = 1, NPF
               II = I*NY
               IPOS = ((I-1)*NUMPROCS + MYID + 1)*NY
               DO J = 1, NY
                  Y(II + J ) = Y(IPOS + J)
                  YPRIME(II + J ) = YPRIME(IPOS + J)
               END DO
            END DO
         END IF
      END IF 
C
C     
C
C     Check to see that TOUT is different from T.
C
      IF(TOUT .EQ. T)GO TO 719
C
C     Check HMAX.
C
      IF(INFO(7) .NE. 0) THEN
         HMAX = RWORK(LHMAX)
         IF (HMAX .LE. 0.0D0) GO TO 710
      END IF
C
C     Initialize counters and other flags.
C
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNSE)=0
      IWORK(LNJE)=0
      IWORK(LETF)=0
      IWORK(LNCFN)=0
      IWORK(LNNI)=0
      IWORK(LNLI)=0
      IWORK(LNLIS)=0
      IWORK(LNPS)=0
      IWORK(LNCFL)=0
C      IWORK(LKPRIN)=INFO(18)
      IDID=1
      IF (INFO(20) .EQ. 5 .OR. (INFO(5).EQ.2.AND.INFO(12).EQ.0)) THEN
         LIADF = LIW - 3*NY + 1
         IF (INFO(20) .EQ. 5)  LIADF = LIADF - INFO(22)
C
C     Initialize the ADIFOR work space for SparsLinC option
         DO I = LIADF, LIW
            IWORK(I) = 0
         END DO
         LASTADFSP = LIW - LIADF + 1
      ELSE
         LIADF = LIW
      END IF

      GO TO 200
C
C-----------------------------------------------------------------------
C     This block is for continuation calls only.
C     Here we check INFO(1), and if the last step was interrupted,
C     we check whether appropriate action was taken.
C-----------------------------------------------------------------------
C
100   CONTINUE
      NY = IWORK(LNY)
      TFINAL = RWORK(LTFN)
C.....set IWORK array pointers
      LID = LICNS
      IF (INFO(10).EQ.1 .OR. INFO(10).EQ.3) LID = LID + NY
C.....Pointer for the ISENWK work array
C
      LISE = IWORK(LLCIWP) + IWORK(LLNIWP)
C.....Compute the number of equations for this processors
C
      IF (INFO(24) .GE. 0) THEN
         MYNEQ = IWORK(LNY)*(NPF+1)
      ELSE
C.....Compute the number of equations for this processor
C
         NAVG = (-INFO(24))/NUMPROCS
         IF (MYID .LT. MOD(-INFO(24), NUMPROCS)) THEN
            NPB = NAVG + 1
         ELSE
            NPB = NAVG
         END IF
         MYNEQ = IWORK(LNY)*NPB
      END IF
      IF(INFO(1).EQ.1)GO TO 110
      ITEMP = 1
      IF(INFO(1).NE.-1)GO TO 701
C
C     If we are here, the last step was interrupted by an error
C     condition from DDSTP, and appropriate action was not taken.
C     This is a fatal error.
C
      MSG = 'DASPK--  THE LAST STEP TERMINATED WITH A NEGATIVE'
      CALL XERRWD(MSG,49,201,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  VALUE (=I1) OF IDID AND NO APPROPRIATE'
      CALL XERRWD(MSG,47,202,0,1,IDID,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  ACTION WAS TAKEN. RUN TERMINATED'
      CALL XERRWD(MSG,41,203,1,0,0,0,0,0.0D0,0.0D0)
      RETURN
110   CONTINUE
C
C-----------------------------------------------------------------------
C     This block is executed on all calls.
C
C     Counters are saved for later checks of performance.
C     Then the error tolerance parameters are checked, and the
C     work array pointers are set.
C-----------------------------------------------------------------------
C
200   CONTINUE
C
C     Save counters for use later.
C
      IWORK(LNSTL)=IWORK(LNST)
      NLI0 = IWORK(LNLI)
      NNI0 = IWORK(LNNI)
      NCFN0 = IWORK(LNCFN)
      NCFL0 = IWORK(LNCFL)
      NWARN = 0
C
C     Check RTOL and ATOL.
C
      NZFLG = 0
      IF (INFO(24) .LT. 0) THEN
         NP = NPB-1
      ELSE
         NP = NPF
      END IF
      IF (INFO(2) .EQ. 0) THEN
         RTOLI = RTOL(1)
         ATOLI = ATOL(1)
         IF (RTOLI .GT. 0.0D0 .OR. ATOLI .GT. 0.0D0) NZFLG = 1
         IF (RTOLI .LT. 0.0D0) GO TO 706
         IF (ATOLI .LT. 0.0D0) GO TO 707
      ELSE
         DO 210 I=1,NY
            RTOLI = RTOL(I)
            ATOLI = ATOL(I)
            IF (RTOLI .GT. 0.0D0 .OR. ATOLI .GT. 0.0D0) NZFLG = 1
            IF (RTOLI .LT. 0.0D0) GO TO 706
            IF (ATOLI .LT. 0.0D0) GO TO 707
 210     CONTINUE         
         DO IP =1, NP
            DO I = 1, NY
               II = IP*NY + I
               RTOLI = RTOL(II)
               ATOLI = ATOL(II)
               IF (RTOLI .LT. 0.0D0) GO TO 706
               IF (ATOLI .LT. 0.0D0) GO TO 707
            END DO
         END DO
      END IF
      IF (NZFLG .EQ. 0) GO TO 708
C
C     Set pointer for ADIFOR with SparcLinC option in IWORK segments
C
      IF (INFO(20) .EQ. 5 .OR. (INFO(5).EQ.2.AND.INFO(12).EQ.0)) THEN
         LIADF = LIW - 3*NY + 1
         IF (INFO(20) .EQ. 5)  LIADF = LIADF - INFO(22)
C
C     Initialize the ADIFOR work space
         IF (LASTADFSP .LT. LIW - LIADF + 1) THEN
            DO I = LIADF, LIW - LASTADFSP
               IWORK(I) = 0
            END DO
            LASTADFSP = LIW - LIADF + 1
         END IF
      ELSE
         LIADF = LIW
      END IF
C
C     Set pointers to RWORK segments.
C     For direct methods, SAVR is not used.
C
      IF (INFO(24).GE.0 .OR. (INFO(24).LT.0.AND.
     *     INFO(30).NE.2.AND.INFO(30).NE.6) ) THEN
         NPHI = MYNEQ
         LSAVR = LDELTA
         IF (INFO(12) .EQ. 1) LSAVR = LDELTA + MYNEQ
         LE  = LSAVR + MYNEQ
         LWT = LE + MYNEQ
         LVT = LWT
         IF (INFO(16) .NE. 0 .OR. 
     *        (INFO(18).EQ.1.AND.INFO(28).GT.0)) LVT = LWT + MYNEQ
         LPHI = LVT + MYNEQ
         IF (INFO(11) .EQ. 0) THEN
            LWM  = LPHI + (IWORK(LMXORD)+1)*MYNEQ
         ELSE
            IF (INFO(12) .EQ. 0) THEN
               LWM = LPHI + MAX0(IWORK(LMXORD)+1,4)*MYNEQ
            ELSE
               LWM = LPHI + MAX0(IWORK(LMXORD)+1,5)*MYNEQ
            END IF
         END IF
      ELSE
C     
C     For the adjoint DAE case (F_y'\lambda)' + F_y\lmabda = 0
C
         NPHI = 2*MYNEQ
         IF (INFO(12) .EQ. 0) THEN
            LSAVR = LDELTA
            LE  = LSAVR + NPHI
         ELSE
            LSAVR = LDELTA + NPHI
            LE  = LSAVR + MYNEQ
         END IF
         LWT = LE + 2*MYNEQ
         LVT = LWT
         IF (INFO(16) .NE. 0 .OR. 
     *        (INFO(18).EQ.1.AND.INFO(28).GT.0)) LVT = LWT + MYNEQ
         LPHI = LVT + NPHI
         IF (INFO(11) .EQ. 0) THEN
            LPHI_LEN = (IWORK(LMXORD)+1)*MYNEQ
         ELSE
            IF (INFO(12) .EQ. 0) THEN
               LPHI_LEN = MAX0(IWORK(LMXORD)+1,4)*MYNEQ
            ELSE
               LPHI_LEN = MAX0(IWORK(LMXORD)+1,5)*MYNEQ
            END IF
         END IF
         LWM  = LPHI + 2*LPHI_LEN
      END IF
      LSE = LWM + IWORK(LNPD)
C
C.....Saved information for Adjoint method
      IF (INFO(24) .LT. 0) THEN
         LADI = IWORK(LNRW) - 6*NY + 1
         IF (INFO(11) .GT. 0 .AND. INFO(11) .LT. 6) THEN
C
C...  Require G_Y information in ADRES evaluation during initialization
            LADI = LADI - MYNEQ
         END IF
         IF (INFO(30) .EQ. 2 .OR. INFO(11) .EQ. 5) LADI = LADI-NY
      ELSE
         LADI = IWORK(LNRW)
      END IF
C      
      NYMNQ = NY - INFO(28)
      IF (INFO(1) .EQ. 1) GO TO 400
C
C-----------------------------------------------------------------------
C     This block is executed on the initial call only.
C     Set the initial step size, the error weight vector, and PHI.
C     Compute unknown initial components of Y and YPRIME, if requested.
C-----------------------------------------------------------------------
C
300   CONTINUE
      TN=T
      IDID=1
C
C     Set error weight array WT and altered weight array VT.
C
      CALL DDAWTS(MYNEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      CALL DINVWT(MYNEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) GO TO 713
      IF (LVT .GT. LWT) CALL DCOPY(MYNEQ, RWORK(LWT), 1, RWORK(LVT), 1)
      IF (INFO(16) .NE. 0) THEN
        DO I = 1, NY
           IF (IWORK(LID+I-1) .LT. 0) THEN
C
C...  Exclude algebraic variables from the error test
              RWORK(LVT+I-1) = 0.0D0
           END IF
        END DO
        DO J = 1, NP
           DO I = 1, NY
              IF (IWORK(LID+I-1) .LT. 0) THEN
                 RWORK(J*NY+LVT+I-1) = 0.0D0
              END IF
           END DO
        END DO
      ENDIF
      IF (INFO(18) .EQ. 1 .AND. INFO(28) .GT. 0) THEN
C
C...  Exclude quadrature variables from the error test
         DO J = 1, NP + 1
            JJ = LVT + (J-1)*NY
            DO I = NY-INFO(28), NY-1
               RWORK(JJ + I) = 0.0D0
            END DO
         END DO
      END IF
C
C     Compute unit roundoff and HMIN.
C
      UROUND = D1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
C
C     Set/check STPTOL control for initial condition calculation.
C     
      IF (INFO(11) .NE. 0) THEN
         IF( INFO(17) .EQ. 0) THEN
            RWORK(LSTOL) = UROUND**.6667D0
         ELSE
            IF (RWORK(LSTOL) .LE. 0.0D0) GO TO 725
         ENDIF
      ENDIF
C
C     Compute EPCON and square root of NY and its reciprocal, used
C     inside iterative solver.
C
      RWORK(LEPCON) = 0.33D0
      FLOATN = NY
      RWORK(LSQRN) = SQRT(FLOATN)
      RWORK(LRSQRN) = 1.D0/RWORK(LSQRN)
      IF (INFO(21) .EQ. 0) RWORK(LPRT) = 1.0D-3
C
C     Check initial interval to see that it is long enough.
C
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     Check H0, if this was input.
C
      IF (INFO(8) .EQ. 0) GO TO 310
      H0 = RWORK(LH)
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 711
      IF (H0 .EQ. 0.0D0) GO TO 712
      GO TO 320
 310  CONTINUE
C
C     Compute initial stepsize, to be used by either
C     DDSTP or DDASIC, depending on INFO(11).
C
      H0 = 0.001D0*TDIST
      YPNORM = DDWNRM(NY,YPRIME,RWORK(LVT),RPAR,IPAR)
      IF (INFO(23) .EQ. 0) THEN
         DO I = 1, NP
            II = I*NY + 1     
            TNORM2 = DDWNRM(NY,YPRIME(II),RWORK(LVT+II-1),RPAR,IPAR)
            IF (TNORM2 .GT. YPNORM) YPNORM = TNORM2
         END DO
      END IF
      IF (YPNORM .GT. 0.5D0/H0) H0 = 0.5D0/YPNORM
      H0 = SIGN(H0,TOUT-T)
C
C     Adjust H0 if necessary to meet HMAX bound.
C
320   IF (INFO(7) .EQ. 0) GO TO 330
      RH = ABS(H0)/RWORK(LHMAX)
      IF (RH .GT. 1.0D0) H0 = H0/RH
C
C     Check against TSTOP, if applicable.
C
330   IF (INFO(4) .EQ. 0) GO TO 340
      TSTOP = RWORK(LTSTOP)
      TN = T
      IF ((TSTOP - T)*H0 .LE. 0.0D0) GO TO 715
      IF ((T + H0 - TSTOP)*H0 .GT. 0.0D0) H0 = TSTOP - T
      IF ((TSTOP - TOUT)*H0 .LT. 0.0D0) GO TO 709
C
340   CONTINUE
      IF (INFO(11) .EQ. 0) GO TO 370
C
C     Compute unknown components of initial Y and YPRIME, depending
C     on INFO(11) and INFO(12).  INFO(12) represents the nonlinear
C     solver type (direct/Krylov).  Pass the name of the specific 
C     nonlinear solver, depending on INFO(12).  The location of the work
C     arrays SAVR, YIC, YPIC, PWK also differ in the two cases.
C
      NWT = 1
      EPCONI = RWORK(LEPIN)*RWORK(LEPCON)
C
C     Detect the index-0 case
C
      INDEX0 = 0
      IF (INFO(11) .EQ. 1) THEN
         INDEX0 = 1
         DO I = 1, NY
            IF (IWORK(LID+I-1) .LT. 0) THEN
               INDEX0 = 0
               GOTO 350
            END IF
         END DO
      END IF
350   CONTINUE
      INFO(19) = NPF
      IF(INFO(24) .LT. 0) INFO(24) = -NPB
      ICOPT = INFO(11)
      IF (INFO(12) .EQ. 0) THEN
         LYIC = LPHI + 2*MYNEQ
         LYPIC = LYIC + MYNEQ
         LPWK = LYPIC
         CALL DDASIC(TN,Y,YPRIME,MYNEQ,ICOPT,IWORK(LID),
     *     RES,JAC,PSOL,H0,RWORK(LWT),NWT,IDID,RPAR,IPAR,
     *     RWORK(LPHI),NPHI,RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *     RWORK(LYIC),RWORK(LYPIC),RWORK(LPWK),RWORK(LWM),IWORK(LIWM),
     *     HMIN,RWORK(LROUND),RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *     EPCONI,RWORK(LSTOL),INFO(15),ICNFLG,IWORK(LICNS),DDASID,
     *     INFO(19),RWORK(LSE),IWORK(LISE),RWORK(LPRT),SENPAR,G_RES,
     *     A_RES, K_RES, T_RES, LIADF, RWORK(LADI), TFINAL, 
     *     INDEX0,TDIST)
      ELSE IF (INFO(12) .EQ. 1) THEN
C***         LYIC = LWM
         LYIC = LPHI + 2*MYNEQ
         LYPIC = LYIC + MYNEQ
         LPWK = LYPIC + MYNEQ
         CALL DDASIC(TN,Y,YPRIME,MYNEQ,ICOPT,IWORK(LID),
     *     RES,JAC,PSOL,H0,RWORK(LWT),NWT,IDID,RPAR,IPAR,
     *     RWORK(LPHI),NPHI,RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *     RWORK(LYIC),RWORK(LYPIC),RWORK(LPWK),RWORK(LWM),IWORK(LIWM),
     *     HMIN,RWORK(LROUND),RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *     EPCONI,RWORK(LSTOL),INFO(15),ICNFLG,IWORK(LICNS),DDASIK,
     *     INFO(19),RWORK(LSE),IWORK(LISE),RWORK(LPRT),SENPAR,G_RES,
     *     A_RES,K_RES, T_RES, LIADF, RWORK(LADI), TFINAL,
     *     INDEX0,TDIST)
      ENDIF
      INFO(19) = INFO19
      INFO(24) = INFO24
C
      IF (IDID .EQ. -12) GO TO 600
      IF (IDID .LT. 0) THEN
         MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
         CALL XERRWD(MSG,44,685,0,0,0,0,0,0.0D0,0.0D0)
         MSG = 'DASPK--  INITIAL (Y,YPRIME) COULD NOT BE COMPUTED'
         CALL XERRWD(MSG,49,686,0,0,0,0,2,TN,H0)
         GO TO 600
      END IF
C
C     DDASIC was successful.  If this was the first call to DDASIC,
C     update the WT array (with the current Y) and call it again.
C
      IF (NWT .EQ. 2) GO TO 355
      NWT = 2
      DO I = 1, MYNEQ
         YPRIME(I) = 0.0D0
      END DO
      CALL DDAWTS(MYNEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      CALL DINVWT(MYNEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) GO TO 713
      GO TO 350
C
C     If INFO(14) = 1, return now with IDID = 4.
C
355   CONTINUE
      IF (IDID .EQ. 5) THEN
        MSG = 'DASPK-- Warning: Second stages of initialization failed.'
        CALL XERRWD (MSG, 56, 501, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      END IF         
      IF (INFO(14) .EQ. 1) THEN
        IF (IDID .NE. 5) IDID = 4
        H = H0
        IF (INFO(11) .GT. 0) RWORK(LHOLD) = H0
        GO TO 590
      ENDIF
C
C     Update the WT and VT arrays one more time, with the new Y.
C
      CALL DDAWTS(MYNEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      CALL DINVWT(MYNEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) GO TO 713
      IF (LVT .GT. LWT) CALL DCOPY(MYNEQ, RWORK(LWT), 1, RWORK(LVT), 1)
      IF (INFO(16) .NE. 0) THEN
        DO I = 1, NY
           IF (IWORK(LID+I-1) .LT. 0) THEN
C
C...  Exclude algebraic variables from the error test
              RWORK(LVT+I-1) = 0.0D0
           END IF
        END DO
        DO J = 1, NP
           DO I = 1, NY
              IF (IWORK(LID+I-1) .LT. 0) THEN
                 RWORK(J*NY+LVT+I-1) = 0.0D0
              END IF
           END DO
        END DO
      ENDIF
      IF (INFO(18) .EQ. 1 .AND. INFO(28) .GT. 0) THEN
C
C...  Exclude quadrature variables from the error test
         DO J = 1, NP + 1
            JJ = LVT + (J-1)*NY
            DO I = NY-INFO(28), NY-1
               RWORK(JJ + I) = 0.0D0
            END DO
         END DO
      END IF
C
C     Reset the initial stepsize to be used by DDSTP.
C     Use H0, if this was input.  Otherwise, recompute H0,
C     and adjust it if necessary to meet HMAX bound.
C
      IF (INFO(8) .NE. 0) THEN
         H0 = RWORK(LH)
         GO TO 360
      ENDIF
C
      H0 = 0.001D0*TDIST
      YPNORM = DDWNRM(NY,YPRIME,RWORK(LVT),RPAR,IPAR)
      IF (INFO(23) .EQ. 0) THEN
         DO I = 1, NP
            II = I*NY + 1     
            TNORM2 = DDWNRM(NY,YPRIME(II),RWORK(LVT+II-1),RPAR,IPAR)
            IF (TNORM2 .GT. YPNORM) YPNORM = TNORM2
         END DO
      END IF
      IF (YPNORM .GT. 0.5D0/H0) H0 = 0.5D0/YPNORM
      H0 = SIGN(H0,TOUT-T)
C
360   IF (INFO(7) .NE. 0) THEN
         RH = ABS(H0)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H0 = H0/RH
      ENDIF
C
C     Check against TSTOP, if applicable.
C
      IF (INFO(4) .NE. 0) THEN
         TSTOP = RWORK(LTSTOP)
         IF ((T + H0 - TSTOP)*H0 .GT. 0.0D0) H0 = TSTOP - T
      ENDIF
C
C     Load H and RWORK(LH) with H0.
C
370   H = H0
C.....Possible same stepsize control on all processors
C
      IF (INFO(27).GT.1) CALL MPI_STEPSIZE(H)
      RWORK(LH) = H
      IF (NPHI .EQ. 2*MYNEQ) THEN
C
C     For the adjoint DAE case F_y'(t,y)
         DO I = 1, NP+1
            II = (I-1)*NY+1
            II2 = II + MYNEQ
            CALL FXPV(TFINAL-TN,Y(II),CJ,Y(II2),IRES,RPAR,IPAR,RES,
     *           NY, INFO(19),IWORK(LISE),SENPAR,RWORK(LADI), 0)
            CALL FXPV(TFINAL-TN,Y(II),CJ,YPRIME(II2),IRES,RPAR,
     *           IPAR,G_RES,NY,
     *           INFO(19),IWORK(LISE),SENPAR,RWORK(LADI), 1)
            DO J = 1, NY
               YPRIME(II2+J-1) = - YPRIME(II2+J-1)
            END DO
         END DO
C
C     Calculate the weight for the new variables
         CALL DDAWTS(MYNEQ,INFO(2),RTOL,ATOL,Y(MYNEQ+1),
     *        RWORK(LVT+MYNEQ),RPAR,IPAR)
         CALL DINVWT(MYNEQ,RWORK(LVT+MYNEQ),IER)
         IF (IER .NE. 0) GO TO 713
      END IF
C
C     Load Y and H*YPRIME into PHI(*,1) and PHI(*,2).
C
      ITEMP = LPHI + NPHI
      DO 380 I = 1,NPHI
         RWORK(LPHI  + I - 1) = Y(I)
         RWORK(ITEMP + I - 1) = H*YPRIME(I)
 380  CONTINUE
C
      GO TO 500
C
C-----------------------------------------------------------------------
C     This block is for continuation calls only.
C     Its purpose is to check stop conditions before taking a step.
C     Adjust H if necessary to meet HMAX bound.
C-----------------------------------------------------------------------
C
400   CONTINUE
      UROUND=RWORK(LROUND)
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
      RH = ABS(H)/RWORK(LHMAX)
      IF(RH .GT. 1.0D0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
      CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0D0) GO TO 490
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 425
      CALL DDATRP(TN,TN,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
425   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
430   IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
      IF((TN-TSTOP)*H.GE.0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
      CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GE. 0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
      IF((TN-T)*H .LE. 0.0D0) GO TO 450
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 445
      CALL DDATRP(TN,TN,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
445   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
450   CONTINUE
C
C     Check whether we are within roundoff of TSTOP.
C
      IF(ABS(TN-TSTOP).GT.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 460
      CALL DDATRP(TN,TSTOP,Y,YPRIME,MYNEQ,NPHI,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
      H=TSTOP-TN
      RWORK(LH)=H
C
490   IF (DONE) GO TO 590
C
C-----------------------------------------------------------------------
C     The next block contains the call to the one-step integrator DDSTP.
C     This is a looping point for the integration steps.
C     Check for too many steps.
C     Check for poor Newton/Krylov performance.
C     Update WT.  Check for too much accuracy requested.
C     Compute minimum stepsize.
C-----------------------------------------------------------------------
C
500   CONTINUE
C
C     Check for too many steps.
C
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500) GO TO 505
      IDID=-1
      GO TO 527
C
C Check for poor Newton/Krylov performance.
C
505   IF (INFO(12) .EQ. 0) GO TO 510
      NSTD = IWORK(LNST) - IWORK(LNSTL)
      NNID = IWORK(LNNI) - NNI0
      IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 510
      AVLIN = REAL(IWORK(LNLI) - NLI0)/REAL(NNID)
      RCFN = REAL(IWORK(LNCFN) - NCFN0)/REAL(NSTD)
      RCFL = REAL(IWORK(LNCFL) - NCFL0)/REAL(NNID)
      FMAXL = IWORK(LMAXL)
      LAVL = AVLIN .GT. FMAXL
      LCFN = RCFN .GT. 0.9D0
      LCFL = RCFL .GT. 0.9D0
      LWARN = LAVL .OR. LCFN .OR. LCFL
      IF (.NOT.LWARN) GO TO 510
      NWARN = NWARN + 1
      IF (NWARN .GT. 10) GO TO 510
      IF (LAVL) THEN
        MSG = 'DASPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 501, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
        MSG = '      at T = R1. Average no. of linear iterations = R2  '
        CALL XERRWD (MSG, 56, 501, 0, 0, 0, 0, 2, TN, AVLIN)
      ENDIF
      IF (LCFN) THEN
        MSG = 'DASPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 502, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
        MSG = '      at T = R1. Nonlinear convergence failure rate = R2'
        CALL XERRWD (MSG, 56, 502, 0, 0, 0, 0, 2, TN, RCFN)
      ENDIF
      IF (LCFL) THEN
        MSG = 'DASPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 503, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
        MSG = '      at T = R1. Linear convergence failure rate = R2   '
        CALL XERRWD (MSG, 56, 503, 0, 0, 0, 0, 2, TN, RCFL)
      ENDIF
C
C     Update WT and VT, if this is not the first call.
C
510   CALL DDAWTS(MYNEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),RWORK(LWT),
     *            RPAR,IPAR)
      CALL DINVWT(MYNEQ,RWORK(LWT),IER)
      IF (NPHI .EQ. 2*MYNEQ) THEN
C
C     Calculate the weight for the new variables
         CALL DDAWTS(MYNEQ,INFO(2),RTOL,ATOL,RWORK(LPHI+MYNEQ),
     *        RWORK(LVT+MYNEQ),RPAR,IPAR)
         CALL DINVWT(MYNEQ,RWORK(LVT+MYNEQ),IER)
      END IF
      IF (IER .NE. 0) THEN
        IDID = -3
        GO TO 527
      ENDIF
      IF (LVT .GT. LWT) CALL DCOPY(MYNEQ, RWORK(LWT), 1, RWORK(LVT), 1)
      IF (INFO(16) .NE. 0) THEN
        DO I = 1, NY
           IF (IWORK(LID+I-1) .LT. 0) THEN
C
C...  Exclude algebraic variables from the error test
              RWORK(LVT+I-1) = 0.0D0
           END IF
        END DO
        DO J = 1, NP
           DO I = 1, NY
              IF (IWORK(LID+I-1) .LT. 0) THEN
                 RWORK(J*NY+LVT+I-1) = 0.0D0
              END IF
           END DO
        END DO
      ENDIF
      IF (INFO(18) .EQ. 1 .AND. INFO(28) .GT. 0) THEN
C
C...  Exclude quadrature variables from the error test
         DO J = 1, NP + 1
            JJ = LVT + (J-1)*NY
            DO I = NY-INFO(28), NY-1
               RWORK(JJ + I) = 0.0D0
            END DO
         END DO
      END IF
C
C     Test for too much accuracy requested.
C
      R = DDWNRM(NY,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)
      DO I = 1, NP
         II = I*NY     
         R2 = DDWNRM(NY,RWORK(LPHI+II),RWORK(LWT+II),RPAR,IPAR)
         IF (R2 .GT. R) R = R2
      END DO
      R = R*100.0D0*UROUND
      IF (R .LE. 1.0D0) GO TO 525
C
C     Multiply RTOL and ATOL by R and return.
C
      IF(INFO(2).EQ.1)GO TO 523
      RTOL(1)=R*RTOL(1)
      ATOL(1)=R*ATOL(1)
      IDID=-2
      GO TO 527
 523  DO 524 I=1,MYNEQ
         RTOL(I)=R*RTOL(I)
 524     ATOL(I)=R*ATOL(I)
      IDID=-2
      GO TO 527
 525  CONTINUE
C
C     Compute minimum stepsize.
C
      HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
C
C     Test H vs. HMAX
      IF (INFO(7) .NE. 0) THEN
         RH = ABS(H)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H = H/RH
      ENDIF
C
C     Call the one-step integrator.
C     Note that INFO(12) represents the nonlinear solver type.
C     Pass the required nonlinear solver, depending upon INFO(12).
C
      INFO(19) = NPF
      IF (INFO(24) .LT. 0) INFO(24) = -NPB
      IF (INFO(12) .EQ. 0) THEN
         CALL DDSTP(TN,Y,YPRIME,MYNEQ,
     *      RES,JAC,PSOL,H,RWORK(LWT),RWORK(LVT),INFO(1),IDID,RPAR,IPAR,
     *      RWORK(LPHI),NPHI,RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *      RWORK(LWM),IWORK(LIWM),
     *      RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *      RWORK(LPSI),RWORK(LSIGMA),
     *      RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN,
     *      RWORK(LROUND), RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *      RWORK(LEPCON), IWORK(LPHASE),IWORK(LJCALC),INFO(15),
     *      IWORK(LK), IWORK(LKOLD),IWORK(LNS),NONNEG,INFO(12),
     *      DNEDD,
     *      INFO(19),RWORK(LSE),IWORK(LISE),RWORK(LPRT),SENPAR,G_RES,
     *      A_RES, K_RES, LIADF, RWORK(LADI), TFINAL)
      ELSE IF (INFO(12) .EQ. 1) THEN
         CALL DDSTP(TN,Y,YPRIME,MYNEQ,
     *      RES,JAC,PSOL,H,RWORK(LWT),RWORK(LVT),INFO(1),IDID,RPAR,IPAR,
     *      RWORK(LPHI),NPHI,RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *      RWORK(LWM),IWORK(LIWM),
     *      RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *      RWORK(LPSI),RWORK(LSIGMA),
     *      RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN,
     *      RWORK(LROUND), RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *      RWORK(LEPCON), IWORK(LPHASE),IWORK(LJCALC),INFO(15),
     *      IWORK(LK), IWORK(LKOLD),IWORK(LNS),NONNEG,INFO(12),
     *      DNEDK,
     *      INFO(19),RWORK(LSE),IWORK(LISE),RWORK(LPRT),SENPAR,G_RES,
     *      A_RES, K_RES, LIADF, RWORK(LADI), TFINAL)  
      ENDIF
      INFO(19) = INFO19
      INFO(24) = INFO24
C
527   IF(IDID.LT.0)GO TO 600
C
C-----------------------------------------------------------------------
C     This block handles the case of a successful return from DDSTP
C     (IDID=1).  Test for stop conditions.
C-----------------------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
      IF(INFO(3).NE.0)GO TO 530
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
      CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=3
      T=TOUT
      GO TO 580
 530  IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
      T=TN
      IDID=1
      GO TO 580
 535  CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=3
      T=TOUT
      GO TO 580
 540  IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LE.0.0D0)GO TO 542
      CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
      GO TO 580
542   IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
      H=TSTOP-TN
      GO TO 500
545   CALL DDATRP(TN,TSTOP,Y,YPRIME,MYNEQ,NPHI,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GT.0.0D0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   CALL DDATRP(TN,TSTOP,Y,YPRIME,MYNEQ,NPHI,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
555   CALL DDATRP(TN,TOUT,Y,YPRIME,MYNEQ,NPHI,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
580   CONTINUE
C
C-----------------------------------------------------------------------
C     All successful returns from DDASPK are made from this block.
C-----------------------------------------------------------------------
C
590   CONTINUE
C
C-----------------------------------------------------------------------
C     Collect the sensitivity for MPI parallel machine
C-----------------------------------------------------------------------
c$$$      NP = INFO(19)/numprocs
c$$$      NY = NEQ/(INFO(19)+1)
c$$$      do i = 1, np
c$$$        ips = i*NY+1
c$$$        ipr = NY + (i-1)*numprocs*NY+1
c$$$        call MPI_gather(u(ips),NY,MPI_double_precision,unew(ipr),NY,
c$$$     *                 MPI_double_precision, root, MPI_COMM_WORLD)
c$$$      end do
c$$$      if (np*numprocs .lt. INFO(19)) then
c$$$c
c$$$c...  for the rest of the sensitivities
c$$$        if (myid .lt. mod(INFO(19),numprocs)) then
c$$$           ips = (np+1)*NY+1
c$$$           lsend = NY
c$$$        else
c$$$           ips = NY+1
c$$$           lsend = 0
c$$$        end if
c$$$        ipr = NY + NP*numprocs*NY+1
c$$$        call MPI_gather(u(ips),lsend,MPI_double_precision,unew(ipr),
c$$$     *             lsend,MPI_double_precision, root, MPI_COMM_WORLD)
c$$$      end if

      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     This block handles all unsuccessful returns other than for
C     illegal input.
C-----------------------------------------------------------------------
C
600   CONTINUE
      ITEMP = -IDID
      GO TO (610,620,630,700,655,640,650,660,670,675,
     *  680,685,690,695), ITEMP
C
C     The maximum number of steps was taken before
C     reaching tout.
C
610   MSG = 'DASPK--  AT CURRENT T (=R1)  500 STEPS'
      CALL XERRWD(MSG,38,610,0,0,0,0,1,TN,0.0D0)
      MSG = 'DASPK--  TAKEN ON THIS CALL BEFORE REACHING TOUT'
      CALL XERRWD(MSG,48,611,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Too much accuracy for machine precision.
C
620   MSG = 'DASPK--  AT T (=R1) TOO MUCH ACCURACY REQUESTED'
      CALL XERRWD(MSG,47,620,0,0,0,0,1,TN,0.0D0)
      MSG = 'DASPK--  FOR PRECISION OF MACHINE. RTOL AND ATOL'
      CALL XERRWD(MSG,48,621,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  WERE INCREASED TO APPROPRIATE VALUES'
      CALL XERRWD(MSG,45,622,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     WT(I) .LE. 0.0D0 for some I (not at start of problem).
C
630   MSG = 'DASPK--  AT T (=R1) SOME ELEMENT OF WT'
      CALL XERRWD(MSG,38,630,0,0,0,0,1,TN,0.0D0)
      MSG = 'DASPK--  HAS BECOME .LE. 0.0'
      CALL XERRWD(MSG,28,631,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Error test failed repeatedly or with H=HMIN.
C
640   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,640,0,0,0,0,2,TN,H)
      MSG='DASPK--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN'
      CALL XERRWD(MSG,57,641,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Nonlinear solver failed to converge repeatedly or with H=HMIN.
C
650   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,650,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  NONLINEAR SOLVER FAILED TO CONVERGE'
      CALL XERRWD(MSG,44,651,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  REPEATEDLY OR WITH ABS(H)=HMIN'
      CALL XERRWD(MSG,40,652,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     The preconditioner had repeated failures.
C
655   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,655,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  PRECONDITIONER HAD REPEATED FAILURES.'
      CALL XERRWD(MSG,46,656,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     The iteration matrix is singular.
C
660   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,660,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  ITERATION MATRIX IS SINGULAR.'
      CALL XERRWD(MSG,38,661,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Nonlinear system failure preceded by error test failures.
C
670   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,670,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  NONLINEAR SOLVER COULD NOT CONVERGE.'
      CALL XERRWD(MSG,45,671,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  ALSO, THE ERROR TEST FAILED REPEATEDLY.'
      CALL XERRWD(MSG,49,672,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Nonlinear system failure because IRES = -1.
C
675   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,675,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  NONLINEAR SYSTEM SOLVER COULD NOT CONVERGE'
      CALL XERRWD(MSG,51,676,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  BECAUSE IRES WAS EQUAL TO MINUS ONE'
      CALL XERRWD(MSG,44,677,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Failure because IRES = -2.
C
680   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2)'
      CALL XERRWD(MSG,40,680,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  IRES WAS EQUAL TO MINUS TWO'
      CALL XERRWD(MSG,36,681,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Failed to compute initial YPRIME.
C
685   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,685,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  INITIAL (Y,YPRIME) COULD NOT BE COMPUTED'
      CALL XERRWD(MSG,49,686,0,0,0,0,2,TN,H0)
      GO TO 700
C
C     Failure because IER was negative from PSOL.
C
690   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2)'
      CALL XERRWD(MSG,40,690,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  IER WAS NEGATIVE FROM PSOL'
      CALL XERRWD(MSG,35,691,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Failure because the linear system solver could not converge.
C
695   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,695,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  LINEAR SYSTEM SOLVER COULD NOT CONVERGE.'
      CALL XERRWD(MSG,50,696,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C
700   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     This block handles all error returns due to illegal input,
C     as detected before calling DDSTP.
C     First the error message routine is called.  If this happens
C     twice in succession, execution is terminated.
C-----------------------------------------------------------------------
C
701   MSG = 'DASPK--  ELEMENT (=I1) OF INFO VECTOR IS NOT VALID'
      CALL XERRWD(MSG,50,1,0,1,ITEMP,0,0,0.0D0,0.0D0)
      GO TO 750
702   MSG = 'DASPK--  NEQ (=I1) .LE. 0'
      CALL XERRWD(MSG,25,2,0,1,NEQ,0,0,0.0D0,0.0D0)
      GO TO 750
703   MSG = 'DASPK--  MAXORD (=I1) NOT IN RANGE'
      CALL XERRWD(MSG,34,3,0,1,MXORD,0,0,0.0D0,0.0D0)
      GO TO 750
704   MSG='DASPK--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)'
      CALL XERRWD(MSG,60,4,0,2,LENRW,LRW,0,0.0D0,0.0D0)
      GO TO 750
705   MSG='DASPK--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)'
      CALL XERRWD(MSG,60,5,0,2,LENIW,LIW,0,0.0D0,0.0D0)
      GO TO 750
706   MSG = 'DASPK--  SOME ELEMENT OF RTOL IS .LT. 0'
      CALL XERRWD(MSG,39,6,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
707   MSG = 'DASPK--  SOME ELEMENT OF ATOL IS .LT. 0'
      CALL XERRWD(MSG,39,7,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
708   MSG = 'DASPK--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO'
      CALL XERRWD(MSG,47,8,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
709   MSG='DASPK--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)'
      CALL XERRWD(MSG,54,9,0,0,0,0,2,TSTOP,TOUT)
      GO TO 750
710   MSG = 'DASPK--  HMAX (=R1) .LT. 0.0'
      CALL XERRWD(MSG,28,10,0,0,0,0,1,HMAX,0.0D0)
      GO TO 750
711   MSG = 'DASPK--  TOUT (=R1) BEHIND T (=R2)'
      CALL XERRWD(MSG,34,11,0,0,0,0,2,TOUT,T)
      GO TO 750
712   MSG = 'DASPK--  INFO(8)=1 AND H0=0.0'
      CALL XERRWD(MSG,29,12,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
713   MSG = 'DASPK--  SOME ELEMENT OF WT IS .LE. 0.0'
      CALL XERRWD(MSG,39,13,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
714   MSG='DASPK-- TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION'
      CALL XERRWD(MSG,60,14,0,0,0,0,2,TOUT,T)
      GO TO 750
715   MSG = 'DASPK--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)'
      CALL XERRWD(MSG,49,15,0,0,0,0,2,TSTOP,TN)
      GO TO 750
717   MSG = 'DASPK--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ'
      CALL XERRWD(MSG,52,17,0,1,IWORK(LML),0,0,0.0D0,0.0D0)
      GO TO 750
718   MSG = 'DASPK--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ'
      CALL XERRWD(MSG,52,18,0,1,IWORK(LMU),0,0,0.0D0,0.0D0)
      GO TO 750
719   MSG = 'DASPK--  TOUT (=R1) IS EQUAL TO T (=R2)'
      CALL XERRWD(MSG,39,19,0,0,0,0,2,TOUT,T)
      GO TO 750
720   MSG = 'DASPK--  MAXL (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. NEQ'
      CALL XERRWD(MSG,54,20,0,1,IWORK(LMAXL),0,0,0.0D0,0.0D0)
      GO TO 750
721   MSG = 'DASPK--  KMP (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. MAXL'
      CALL XERRWD(MSG,54,21,0,1,IWORK(LKMP),0,0,0.0D0,0.0D0)
      GO TO 750
722   MSG = 'DASPK--  NRMAX (=I1) ILLEGAL. .LT. 0'
      CALL XERRWD(MSG,36,22,0,1,IWORK(LNRMAX),0,0,0.0D0,0.0D0)
      GO TO 750
723   MSG = 'DASPK--  EPLI (=R1) ILLEGAL. EITHER .LE. 0.D0 OR .GE. 1.D0'
      CALL XERRWD(MSG,58,23,0,0,0,0,1,RWORK(LEPLI),0.0D0)
      GO TO 750
724   MSG = 'DASPK--  ILLEGAL IWORK VALUE FOR INFO(11 or 16) .NE. 0'
      CALL XERRWD(MSG,48,24,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
725   MSG = 'DASPK--  ONE OF THE INPUTS FOR INFO(17) = 1 IS ILLEGAL'
      CALL XERRWD(MSG,54,25,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
726   MSG = 'DASPK--  ILLEGAL IWORK VALUE FOR INFO(10) .NE. 0'
      CALL XERRWD(MSG,48,26,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
727   MSG = 'DASPK--  Y(I) AND IWORK(40+I) (I=I1) INCONSISTENT'
      CALL XERRWD(MSG,49,27,0,1,IRET,0,0,0.0D0,0.0D0)
      GO TO 750
728   MSG = 'DASPK-- NEQ MUST EQUAL #OF STATE VAR.*(#OF SENS.VAR.+1)'
      CALL XERRWD(MSG,58,28,0,0,0,0,1,0.D0,0.0D0)
      GO TO 750
 730  MSG = 'DASPK-- When INFO(25)=0, INFO(20) SHOULD LESS THAN 5.' 
      CALL XERRWD(MSG,53,30,0,0,0,0,1,0.D0,0.0D0)
      GO TO 750
 731  MSG = 'DASPK-- When INFO(27)>1, INFO(20) SHOULD BE 4.'
      CALL XERRWD(MSG,46,31,0,0,0,0,1,0.D0,0.0D0)
      GO TO 750
750   IF(INFO(1).EQ.-1) GO TO 760
      INFO(1)=-1
      IDID=-33
      RETURN
760   MSG = 'DASPK--  REPEATED OCCURRENCES OF ILLEGAL INPUT'
      CALL XERRWD(MSG,46,701,0,0,0,0,0,0.0D0,0.0D0)
770   MSG = 'DASPK--  RUN TERMINATED. APPARENT INFINITE LOOP'
      CALL XERRWD(MSG,47,702,1,0,0,0,0,0.0D0,0.0D0)
      RETURN
C
C------END OF SUBROUTINE DDASPK-----------------------------------------
      END
      SUBROUTINE DDASIC (X, Y, YPRIME, NEQ, ICOPT, ID, RES, JAC, PSOL,
     *   H, WT, NIC, IDID, RPAR, IPAR, PHI,NPHI,SAVR, DELTA,E,YIC,YPIC,
     *   PWK, WM, IWM, HMIN, UROUND, EPLI, SQRTN, RSQRTN, EPCONI,
     *   STPTOL, JFLG, ICNFLG, ICNSTR, NLSIC,
     *   ISENFO, SENWRK, ISENWK, CNST, SENPAR, G_RES, A_RES, K_RES,
     *   T_RES, LIADF,ADI, TFINAL, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DDASIC
C***REFER TO  DDASPK
C***DATE WRITTEN   940628   (YYMMDD)
C***REVISION DATE  941206   (YYMMDD)
C***REVISION DATE  950714   (YYMMDD)
C***REVISION DATE  990503   (YYMMDD)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DDASIC is a driver routine to compute consistent initial values
C     for Y and YPRIME.  There are two different options:  
C     Denoting the differential variables in Y by Y_d, and
C     the algebraic variables by Y_a, the problem solved is either:
C     1.  Given Y_d, calculate Y_a and Y_d', or
C     2.  Given Y', calculate Y.
C     3.  Given Y, Y', calculate new pair (Y, Y').
C     In either case, initial values for the given components
C     are input, and initial guesses for the unknown components
C     must also be provided as input.
C
C     The external routine NLSIC solves the resulting nonlinear system.
C
C     The parameters represent
C
C     X  --        Independent variable.
C     Y  --        Solution vector at X.
C     YPRIME --    Derivative of solution vector.
C     NEQ --       Number of equations to be integrated.
C     ICOPT     -- Flag indicating initial condition option chosen.
C                    ICOPT = 1 for option 1 above.
C                    ICOPT = 2 for option 2.
C                    ICOPT = 3 for option 3.
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if option 1 is chosen.
C                    ID(i) = +1 if Y_i is a differential variable,
C                    ID(i) = -1 if Y_i is an algebraic variable. 
C     RES --       External user-supplied routine to evaluate the
C                  residual.  See RES description in DDASPK prologue.
C     JAC --       External user-supplied routine to update Jacobian
C                  or preconditioner information in the nonlinear solver
C                  (optional).  See JAC description in DDASPK prologue.
C     PSOL --      External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  See PSOL in DDASPK prologue.
C     H --         Scaling factor in iteration matrix.  DDASIC may 
C                  reduce H to achieve convergence.
C     WT --        Vector of weights for error criterion.
C     NIC --       Input number of initial condition calculation call 
C                  (= 1 or 2).
C     IDID --      Completion code.  See IDID in DDASPK prologue.
C     RPAR,IPAR -- Real and integer parameter arrays that
C                  are used for communication between the
C                  calling program and external user routines.
C                  They are not altered by DNSK
C     PHI --       Work space for DDASIC of length at least 2*NEQ.
C     SAVR --      Work vector for DDASIC of length NEQ.
C     DELTA --     Work vector for DDASIC of length NEQ.
C     E --         Work vector for DDASIC of length NEQ.
C     YIC,YPIC --  Work vectors for DDASIC, each of length NEQ.
C     PWK --       Work vector for DDASIC of length NEQ.
C     WM,IWM --    Real and integer arrays storing
C                  information required by the linear solver.
C     EPCONI --    Test constant for Newton iteration convergence.
C     ICNFLG --    Flag showing whether constraints on Y are to apply.
C     ICNSTR --    Integer array of length NEQ with constraint types.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST   --    Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied routine which is generated by
C                  ADIFOR for evaluations of sensitivity equations.
C
C     The other parameters are for use internally by DDASIC.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DCOPY, NLSIC
C
C***END PROLOGUE  DDASIC
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),ID(*),WT(*),PHI(NPHI,*), ADI(*)
      DIMENSION SAVR(*),DELTA(*),E(*),YIC(*),YPIC(*),PWK(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*), ICNSTR(*)
      DIMENSION ISENFO(*), ISENWK(*), SENWRK(*), SENPAR(*)
      EXTERNAL RES, JAC, PSOL, NLSIC, G_RES, A_RES, K_RES, T_RES
C
      PARAMETER (LCFN=15, LNY=37,LMTYPE=4)
      PARAMETER (LMXNH=34)
C
C The following parameters are data-loaded here:
C     RHCUT  = factor by which H is reduced on retry of Newton solve.
C     RATEMX = maximum convergence rate for which Newton iteration
C              is considered converging.
C
      SAVE RHCUT, RATEMX
      DATA RHCUT/0.1D0/, RATEMX/0.8D0/
C
C
C-----------------------------------------------------------------------
C     BLOCK 1.
C     Initializations.
C     JSKIP is a flag set to 1 when NIC = 2 and NH = 1, to signal that
C     the initial call to the JAC routine is to be skipped then.
C     Save Y and YPRIME in PHI.  Initialize IDID, NH, and CJ.
C-----------------------------------------------------------------------
C
C****special treatment for ICOPT = 3, 4
      IF (ICOPT.EQ.3 .OR. ICOPT.EQ.4) RHCUT = 2.0D0
      MXNH = IWM(LMXNH)
      IDID = 1
      NH = 1
      JSKIP = 0
      IF (NIC .EQ. 2) JSKIP = 1
      CALL DCOPY (NEQ, Y, 1, PHI(1,1), 1)
      CALL DCOPY (NEQ, YPRIME, 1, PHI(1,2), 1)
C
      HLOC = H
      IF (ICOPT .EQ. 2) THEN
         CJ = 0.0D0
      ELSE IF (ICOPT.EQ.3.OR.ICOPT.EQ.4.OR.ICOPT.EQ.6
     *        .OR.IWM(LMTYPE).EQ.1.OR.IWM(LMTYPE).EQ.5) THEN
         CJ = 1.0D0/HLOC
      ELSE
C
C    ICOPT = 1 or 5
         IF (H.LT.0.D0) THEN
            HLOC = -1.D0
         ELSE
            HLOC = 1.D0
         END IF
         CJ = HLOC       
      END IF
C
C-----------------------------------------------------------------------
C     BLOCK 2
C     Call the nonlinear system solver to obtain
C     consistent initial values for Y and YPRIME.
C-----------------------------------------------------------------------
C
      IF (ISENFO(6) .LT. 0) THEN
         XLOC = TFINAL - X
      ELSE
         XLOC = X
      END IF
 200  CONTINUE
      CALL NLSIC(XLOC,Y,YPRIME,NEQ,ICOPT,ID,RES,JAC,PSOL,HLOC,WT,JSKIP,
     *   IDID,RPAR,IPAR,SAVR,DELTA,E,YIC,YPIC,PWK,WM,IWM,CJ,UROUND,
     *   EPLI,SQRTN,RSQRTN,EPCONI,RATEMX,STPTOL,JFLG,ICNFLG,ICNSTR,
     *   IERNLS,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, 
     *   K_RES, T_RES, LIADF,ADI, INDEX0, TDIST)
C
      IF (IDID .LT. 0) RETURN
      IF (IERNLS .EQ. 0) THEN 
c.... for index-2 problems...........................
         IF (ICOPT .EQ. 4) THEN
            NY = IWM(LNY)
C     
C     Restore the values of YPRIME(*) for the state variables
C     set YPRIME to zero for sensitivity variables
C     
            ICOPT = 5
            DO I = 1, NY
               TEMP = PHI(I,2)
               PHI(I,2) = YPRIME(I)
               YPRIME(I) = TEMP
            END DO
            DO I = NY+1, NEQ
               PHI(I,2) = YPRIME(I)
               YPRIME(I) = 0.0D0
            END DO
            DO I = 1, NEQ
               PHI(I,1) = Y(I)
            END DO
            GOTO 200
         END IF
         RETURN
      END IF
C
C-----------------------------------------------------------------------
C     BLOCK 3
C     The nonlinear solver was unsuccessful.  Increment NCFN.
C     Return with IDID = -12 if either
C       IERNLS = -1: error is considered unrecoverable,
C       ICOPT = 2: we are doing initialization problem type 2, or
C       NH = MXNH: the maximum number of H values has been tried.
C     Otherwise (problem 1 with IERNLS .GE. 1), reduce H and try again.
C     If IERNLS > 1, restore Y and YPRIME to their original values.
C-----------------------------------------------------------------------
C
      IWM(LCFN) = IWM(LCFN) + 1
      JSKIP = 0
C
      IF (IERNLS .EQ. -1) GO TO 350
      IF (ICOPT .EQ. 2) GO TO 350
      IF (NH .EQ. MXNH) GO TO 350
C
      NH = NH + 1
      HLOC = HLOC*RHCUT
      CJ = 1.0D0/HLOC
C
      IF (IERNLS .EQ. 1) GO TO 200
C
      CALL DCOPY (NEQ, PHI(1,1), 1, Y, 1)
      CALL DCOPY (NEQ, PHI(1,2), 1, YPRIME, 1)
      GO TO 200
C
 350  IDID = -12
      IF (ICOPT .EQ. 5) THEN
         IDID = 5
      END IF
      IF (IERNLS .EQ. 3) IDID = -8
      RETURN
C
C------END OF SUBROUTINE DDASIC-----------------------------------------
      END
      SUBROUTINE DYYPNW (NEQ, Y, YPRIME, CJ, RL, P, ICOPT, ID, 
     *                   YNEW, YPNEW)
C
C***BEGIN PROLOGUE  DYYPNW
C***REFER TO  DLINSK
C***DATE WRITTEN   940830   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DYYPNW calculates the new (Y,YPRIME) pair needed in the
C     linesearch algorithm based on the current lambda value.  It is
C     called by DLINSK and DLINSD.  Based on the ICOPT and ID values,
C     the corresponding entry in Y or YPRIME is updated.
C
C     In addition to the parameters described in the calling programs,
C     the parameters represent
C
C     P      -- Array of length NEQ that contains the current
C               approximate Newton step.
C     RL     -- Scalar containing the current lambda value.
C     YNEW   -- Array of length NEQ containing the updated Y vector.
C     YPNEW  -- Array of length NEQ containing the updated YPRIME
C               vector.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED (NONE)
C
C***END PROLOGUE  DYYPNW
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(*), YPRIME(*), YNEW(*), YPNEW(*), ID(*), P(*)
C
      IF (ICOPT .EQ. 1) THEN
C.... ICOPT = 1
         DO I=1,NEQ
            IF(ID(I) .LT. 0) THEN
               YNEW(I) = Y(I) - RL*P(I)
               YPNEW(I) = YPRIME(I)
            ELSE
               YNEW(I) = Y(I)
               YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
            ENDIF
         END DO
      ELSE IF (ICOPT .EQ. 5) THEN
C.... ICOPT = 5
         DO I=1,NEQ
            IF(ID(I) .EQ. -2) THEN
               YNEW(I) = Y(I) - RL*P(I)
            ELSE
               YNEW(I) = Y(I)
            END IF
            YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
         END DO        
      ELSE IF (ICOPT.EQ.3 .OR. ICOPT.EQ.4) THEN
C.... ICOPT = 3,4
         DO I = 1, NEQ
            IF (ID(I) .LT. 0 .OR. ID(I) .EQ. 3) THEN
               YNEW(I) = Y(I) - RL*P(I)
               YPNEW(I) = YPRIME(I)
            ELSE IF (ID(I) .EQ. 2) THEN                              
               YNEW(I) = Y(I)
               YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
            ELSE 
               YNEW(I) = Y(I) - RL*P(I)
               YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
            END IF
         END DO
      ELSE IF (ICOPT .EQ. 2) THEN
C.... ICOPT = 2    
         DO 20 I = 1,NEQ
            YNEW(I) = Y(I) - RL*P(I)
            YPNEW(I) = YPRIME(I)
 20      CONTINUE
      ELSE IF (ICOPT .EQ. 6) THEN
         DO I = 1, NEQ
            YNEW(I) = Y(I) - RL*P(I)
            YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
         END DO
      ENDIF
      RETURN
C----------------------- END OF SUBROUTINE DYYPNW ----------------------
      END
      SUBROUTINE DDSTP(X,Y,YPRIME,NEQ,
     *  RES,JAC,PSOL,H,WT,VT,
     *  JSTART,IDID,RPAR,IPAR,PHI,NPHI,SAVR,DELTA,E,WM,IWM,
     *  ALPHA,BETA,GAMMA,PSI,SIGMA,CJ,CJOLD,HOLD,S,HMIN,UROUND,
     *  EPLI,SQRTN,RSQRTN,EPCON,IPHASE,JCALC,JFLG,K,KOLD,NS,NONNEG,
     *  NTYPE,NLS,
     *  ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES,A_RES, K_RES,
     *  LIADF, ADI, TFINAL)
C
C***BEGIN PROLOGUE  DDSTP
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940909   (YYMMDD) (Reset PSI(1), PHI(*,2) at 690)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DDSTP solves a system of differential/algebraic equations of 
C     the form G(X,Y,YPRIME) = 0, for one step (normally from X to X+H).
C
C     The methods used are modified divided difference, fixed leading 
C     coefficient forms of backward differentiation formulas.  
C     The code adjusts the stepsize and order to control the local error
C     per step.
C
C
C     The parameters represent
C     X  --        Independent variable.
C     Y  --        Solution vector at X.
C     YPRIME --    Derivative of solution vector
C                  after successful step.
C     NEQ --       Number of equations to be integrated.
C     RES --       External user-supplied routine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JAC --       External user-supplied routine to update
C                  Jacobian or preconditioner information in the
C                  nonlinear solver.  See JAC description in DDASPK
C                  prologue.
C     PSOL --      External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  (This is optional).  See PSOL in DDASPK prologue.
C     H --         Appropriate step size for next step.
C                  Normally determined by the code.
C     WT --        Vector of weights for error criterion used in Newton 
C                  test.
C     VT --        Masked vector of weights used in error test.
C     JSTART --    Integer variable set 0 for
C                  first step, 1 otherwise.
C     IDID --      Completion code returned from the nonlinear solver.
C                  See IDID description in DDASPK prologue.
C     RPAR,IPAR -- Real and integer parameter arrays that
C                  are used for communication between the
C                  calling program and external user routines.
C                  They are not altered by DNSK
C     PHI --       Array of divided differences used by
C                  DDSTP. The length is NEQ*(K+1), where
C                  K is the maximum order.
C     SAVR --      Work vector for DDSTP of length NEQ.
C     DELTA,E --   Work vectors for DDSTP of length NEQ.
C     WM,IWM --    Real and integer arrays storing
C                  information required by the linear solver.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied routine which is generated by
C                  ADIFOR for evaluations of sensitivity equations.
C
C     The other parameters are information
C     which is needed internally by DDSTP to
C     continue from step to step.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   NLS, DDWNRM, DDATRP
C
C***END PROLOGUE  DDSTP
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),VT(*)
      DIMENSION PHI(NPHI,*),SAVR(*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*)
      DIMENSION PSI(*),ALPHA(*),BETA(*),GAMMA(*),SIGMA(*)
      DIMENSION RPAR(*),IPAR(*), ADI(*)
      DIMENSION ISENFO(*),ISENWK(*), SENWRK(*), SENPAR(*)
      EXTERNAL  RES, JAC, PSOL, NLS, G_RES, A_RES, K_RES
C
      PARAMETER (LMXORD=3)
      PARAMETER (LNST=11, LETF=14, LCFN=15, LNY=37)
C
C
C-----------------------------------------------------------------------
C     BLOCK 1.
C     Initialize.  On the first call, set
C     the order to 1 and initialize
C     other variables.
C-----------------------------------------------------------------------
C
C     Initializations for all calls
C
      XOLD=X
      NCF=0
      NEF=0
      NY = IWM(LNY)
      NP = ISENFO(1)
      IF(JSTART .NE. 0) GO TO 120
C
C     If this is the first step, perform
C     other initializations
C
      K=1
      KOLD=0
      HOLD=0.0D0
      PSI(1)=H
      CJ = 1.D0/H
      IPHASE = 0
      NS=0
120   CONTINUE
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 2
C     Compute coefficients of formulas for
C     this step.
C-----------------------------------------------------------------------
200   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN0(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
C
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
C
C     Compute ALPHAS, ALPHA0
C
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/I
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     Compute leading coefficient CJ
C
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     Compute variable stepsize error coefficient CK
C
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     Change PHI to PHI STAR
C
      IF(KP1 .LT. NSP1) GO TO 280      
      DO 270 J=NSP1,KP1
         DO 260 I=1,NPHI
 260        PHI(I,J)=BETA(J)*PHI(I,J)
 270     CONTINUE
 280  CONTINUE
C
C     Update time
C
      X=X+H
C
C     Initialize IDID to 1
C
      IDID = 1
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 3
C     Call the nonlinear system solver to obtain the solution and
C     derivative.
C-----------------------------------------------------------------------
C
      IF (ISENFO(6) .LT. 0) THEN
         XLOC = TFINAL - X
      ELSE
         XLOC = X
      END IF
      CALL NLS(XLOC,Y,YPRIME,NEQ,
     *   RES,JAC,PSOL,H,WT,JSTART,IDID,RPAR,IPAR,PHI,NPHI,GAMMA,
     *   SAVR,DELTA,E,WM,IWM,CJ,CJOLD,CJLAST,S,
     *   UROUND,EPLI,SQRTN,RSQRTN,EPCON,JCALC,JFLG,KP1,
     *   NONNEG,NTYPE,IERNLS,
     *   CK,VT,ENORM,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES,A_RES, K_RES,
     *   LIADF,ADI)
C
      IF(IERNLS .NE. 0 .AND. IERNLS .NE. -2) GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 4
C     Estimate the errors at orders K,K-1,K-2
C     as if constant stepsize was used. Estimate
C     the local error at order K and test
C     whether the current step is successful.
C-----------------------------------------------------------------------
C
C     Estimate errors at orders K,K-1,K-2
C
*      ENORM = DDWNRM(NY,E,VT,RPAR,IPAR)
      IF (ISENFO(5) .EQ. 0 .AND. IERNLS .EQ. 0) THEN
         DO I = 1, ISENFO(1)
            II = I*NY+1
            ENORMTMP = DDWNRM(NY,E(II),VT(II),RPAR,IPAR)
            IF (ENORMTMP .GT. ENORM) ENORM = ENORMTMP
         END DO
      ENDIF
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF (ISENFO(6) .LT. 0) NPB = -ISENFO(6) - 1
      IF(K .EQ. 1)GO TO 430
      IF (ISENFO(6) .LT. 0) THEN
         DO I = 1,NPHI
            DELTA(I) = PHI(I,KP1) + E(I)
         END DO
         ERKM1=DDWNRM(NY,DELTA,VT,RPAR,IPAR)        
         DO I = 1, NPB
            II = I*NY+1
            ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
            IF (ENORMTMP .GT. ERKM1) ERKM1 = ENORMTMP
         END DO
         IF (NPHI .EQ. 2*NEQ) THEN
            DO I = 0, NPB
               II = I*NY+1+NEQ
               ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
               IF (ENORMTMP .GT. ERKM1) ERKM1 = ENORMTMP
            END DO
         END IF
         ERKM1=SIGMA(K)*ERKM1
      ELSE
         DO 405 I = 1,NEQ
 405        DELTA(I) = PHI(I,KP1) + E(I)
         ERKM1=DDWNRM(NY,DELTA,VT,RPAR,IPAR)
         IF (ISENFO(5) .EQ. 0 .AND. IERNLS .EQ. 0) THEN
            DO I = 1, ISENFO(1)
               II = I*NY+1
               ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
               IF (ENORMTMP .GT. ERKM1) ERKM1 = ENORMTMP
            END DO
         ENDIF
         ERKM1=SIGMA(K)*ERKM1
      END IF
      TERKM1 = K*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5D0*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      IF (ISENFO(6) .LT. 0) THEN
         ERKM2=DDWNRM(NY,DELTA,VT,RPAR,IPAR)
         DO I = 1, NPB
            II = I*NY+1
            ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
            IF (ENORMTMP .GT. ERKM2) ERKM2 = ENORMTMP
         END DO
         IF (NPHI .EQ. 2*NEQ) THEN
            DO I = 0, NPB
               II = I*NY+1+NEQ
               ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
               IF (ENORMTMP .GT. ERKM2) ERKM2 = ENORMTMP
            END DO
         END IF
         ERKM2=SIGMA(K-1)*ERKM2
      ELSE
         ERKM2=DDWNRM(NY, DELTA,VT,RPAR,IPAR)
         IF (ISENFO(5) .EQ. 0 .AND. IERNLS .EQ. 0) THEN
            DO I = 1, ISENFO(1)
               II = I*NY+1
               ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
               IF (ENORMTMP .GT. ERKM2) ERKM2 = ENORMTMP
            END DO
         ENDIF
         ERKM2=SIGMA(K-1)*ERKM2
      END IF
      TERKM2 = (K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
C
C     Lower the order
C
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C
C     Calculate the local error for the current step
C     to see if the step was successful
C
430   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 5
C     The step is successful. Determine
C     the best order and stepsize for
C     the next step. Update the differences
C     for the next step.
C-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     Estimate the error at order K+1 unless
C        already decided to lower order, or
C        already using maximum order, or
C        stepsize not constant, or
C        order raised in previous step
C
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NPHI
510      DELTA(I)=E(I)-PHI(I,KP2)
      IF (ISENFO(6) .LT. 0) THEN
         ERKP1 = DDWNRM(NY,DELTA,VT,RPAR,IPAR)
         DO I = 1, NPB
            II = I*NY+1
            ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
            IF (ENORMTMP .GT. ERKP1) ERKP1 = ENORMTMP
         END DO
         IF (NPHI .EQ. 2*NEQ) THEN
            DO I = 0, NPB
               II = I*NY+1+NEQ
               ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
               IF (ENORMTMP .GT. ERKP1) ERKP1 = ENORMTMP
            END DO
         END IF
      ELSE
         ERKP1 = DDWNRM(NY,DELTA,VT,RPAR,IPAR)
         IF (ISENFO(5) .EQ. 0) THEN         
            DO I = 1, ISENFO(1)
               II = I*NY+1
               ENORMTMP = DDWNRM(NY,DELTA(II),VT(II),RPAR,IPAR)
               IF (ENORMTMP .GT. ERKP1) ERKP1 = ENORMTMP
            END DO
         ENDIF
      END IF
      ERKP1 = (1.0D0/(K+2))*ERKP1
      TERKP1 = (K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     Raise order
C
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     Lower order
C
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     If IPHASE = 0, increase order by one and multiply stepsize by
C     factor two
C
545   K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
C
C
C     Determine the appropriate stepsize for
C     the next step.
C
550   HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
555   IF(R .GT. 1.0D0) GO TO 560
      R = MAX(0.5D0,MIN(0.9D0,R))
      HNEW = H*R
560   H=HNEW
C
C
C     Update differences for next step
C
575   CONTINUE
C.....Possible same stepsize control on all processors
C
      IF (ISENFO(9).GT.1) CALL MPI_STEPSIZE(H)
C
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NPHI
 580     PHI(I,KP2)=E(I)
 585  CONTINUE
      DO 590 I=1,NPHI
 590     PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO 595 J1=2,KP1
         J=KP1-J1+1
         DO 595 I=1,NPHI
 595        PHI(I,J)=PHI(I,J)+PHI(I,J+1)
      JSTART = 1
      RETURN
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 6
C     The step is unsuccessful. Restore X,PSI,PHI
C     Determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C-----------------------------------------------------------------------
600   IPHASE = 1
C
C     Restore X,PHI,PSI
C
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NPHI
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C
C     Test whether failure is due to nonlinear solver
C     or error test
C
      IF(IERNLS .EQ. -2 .OR. IERNLS .EQ. 0) GO TO 660  ! error test failure
      IWM(LCFN)=IWM(LCFN)+1
C
C
C     The nonlinear solver failed to converge.
C     Determine the cause of the failure and take appropriate action.
C     If IERNLS .LT. 0, then return.  Otherwise, reduce the stepsize
C     and try again, unless too many failures have occurred.
C
      IF (IERNLS .LT. 0) GO TO 675
      NCF = NCF + 1
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
      IF (IDID .EQ. 1) IDID = -7
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C
C     The nonlinear solver converged, and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
C
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
C
C     On first error test failure, keep current order or lower
C     order by one.  Compute new stepsize based on differences
C     of the solution.
C
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = MAX(0.25D0,MIN(0.9D0,R))
      H = H*R
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     On second error test failure, use the current order or
C     decrease order by one.  Reduce the stepsize by a factor of
C     one quarter.
C
665   IF (NEF .GT. 2) GO TO 670
      K = KNEW
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     On third and subsequent error test failures, set the order to
C     one, and reduce the stepsize by a factor of one quarter.
C
670   K = 1
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C
C
C
C     For all crashes, restore Y to its last value,
C     interpolate to find YPRIME at last X, and return.
C
C     Before returning, verify that the user has not set
C     IDID to a nonnegative value.  If the user has set IDID
C     to a nonnegative value, then reset IDID to be -7, indicating
C     a failure in the nonlinear system solver.
C
675   CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,NPHI,K,PHI,PSI)
      JSTART = 1
      IF (IDID .GE. 0) IDID = -7
      RETURN
C
C
C     Go back and try this step again.  
C     If this is the first step, reset PSI(1) and rescale PHI(*,2).
C
690   IF (KOLD .EQ. 0) THEN
        PSI(1) = H
        DO 695 I = 1,NPHI
695       PHI(I,2) = R*PHI(I,2)
        ENDIF
      GO TO 200
C
C------END OF SUBROUTINE DDSTP------------------------------------------
      END
      SUBROUTINE DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
C
C***BEGIN PROLOGUE  DCNSTR
C***DATE WRITTEN   950808   (YYMMDD)
C***REVISION DATE  950814   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine checks for constraint violations in the proposed 
C new approximate solution YNEW.
C If a constraint violation occurs, then a new step length, TAU,
C is calculated, and this value is to be given to the linesearch routine
C to calculate a new approximate solution YNEW.
C
C On entry:
C
C   NEQ    -- size of the nonlinear system, and the length of arrays
C             Y, YNEW and ICNSTR.
C
C   Y      -- real array containing the current approximate y.
C
C   YNEW   -- real array containing the new approximate y.
C
C   ICNSTR -- INTEGER array of length NEQ containing flags indicating
C             which entries in YNEW are to be constrained.
C             if ICNSTR(I) =  2, then YNEW(I) must be .GT. 0,
C             if ICNSTR(I) =  1, then YNEW(I) must be .GE. 0,
C             if ICNSTR(I) = -1, then YNEW(I) must be .LE. 0, while
C             if ICNSTR(I) = -2, then YNEW(I) must be .LT. 0, while
C             if ICNSTR(I) =  0, then YNEW(I) is not constrained.
C
C   RLX    -- real scalar restricting update, if ICNSTR(I) = 2 or -2,
C             to ABS( (YNEW-Y)/Y ) < FAC2*RLX in component I.
C
C   TAU    -- the current size of the step length for the linesearch.
C
C On return
C
C   TAU    -- the adjusted size of the step length if a constraint
C             violation occurred (otherwise, it is unchanged).  it is
C             the step length to give to the linesearch routine.
C
C   IRET   -- output flag.
C             IRET=0 means that YNEW satisfied all constraints.
C             IRET=1 means that YNEW failed to satisfy all the
C                    constraints, and a new linesearch step
C                    must be computed.
C
C   IVAR   -- index of variable causing constraint to be violated.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(NEQ), YNEW(NEQ), ICNSTR(NEQ)
      SAVE FAC, FAC2, ZERO
      DATA FAC /0.6D0/, FAC2 /0.9D0/, ZERO/0.0D0/
C-----------------------------------------------------------------------
C Check constraints for proposed new step YNEW.  If a constraint has
C been violated, then calculate a new step length, TAU, to be
C used in the linesearch routine.
C-----------------------------------------------------------------------
      IRET = 0
      RDYMX = ZERO
      IVAR = 0
      DO 100 I = 1,NEQ
C
         IF (ICNSTR(I) .EQ. 2) THEN
            RDY = ABS( (YNEW(I)-Y(I))/Y(I) )
            IF (RDY .GT. RDYMX) THEN
               RDYMX = RDY
               IVAR = I
            ENDIF
            IF (YNEW(I) .LE. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ELSEIF (ICNSTR(I) .EQ. 1) THEN
            IF (YNEW(I) .LT. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ELSEIF (ICNSTR(I) .EQ. -1) THEN
            IF (YNEW(I) .GT. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ELSEIF (ICNSTR(I) .EQ. -2) THEN
            RDY = ABS( (YNEW(I)-Y(I))/Y(I) )
            IF (RDY .GT. RDYMX) THEN
               RDYMX = RDY
               IVAR = I
            ENDIF
            IF (YNEW(I) .GE. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ENDIF
 100  CONTINUE

      IF(RDYMX .GE. RLX) THEN
         TAU = FAC2*TAU*RLX/RDYMX
         IRET = 1
      ENDIF
C
      RETURN
C----------------------- END OF SUBROUTINE DCNSTR ----------------------
      END
      SUBROUTINE DCNST0 (NEQ, Y, ICNSTR, IRET)
C
C***BEGIN PROLOGUE  DCNST0
C***DATE WRITTEN   950808   (YYMMDD)
C***REVISION DATE  950808   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine checks for constraint violations in the initial 
C approximate solution u.
C
C On entry
C
C   NEQ    -- size of the nonlinear system, and the length of arrays
C             Y and ICNSTR.
C
C   Y      -- real array containing the initial approximate root.
C
C   ICNSTR -- INTEGER array of length NEQ containing flags indicating
C             which entries in Y are to be constrained.
C             if ICNSTR(I) =  2, then Y(I) must be .GT. 0,
C             if ICNSTR(I) =  1, then Y(I) must be .GE. 0,
C             if ICNSTR(I) = -1, then Y(I) must be .LE. 0, while
C             if ICNSTR(I) = -2, then Y(I) must be .LT. 0, while
C             if ICNSTR(I) =  0, then Y(I) is not constrained.
C
C On return
C
C   IRET   -- output flag.
C             IRET=0    means that u satisfied all constraints.
C             IRET.NE.0 means that Y(IRET) failed to satisfy its
C                       constraint.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(NEQ), ICNSTR(NEQ)
      SAVE ZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C Check constraints for initial Y.  If a constraint has been violated,
C set IRET = I to signal an error return to calling routine.
C-----------------------------------------------------------------------
      IRET = 0
      DO 100 I = 1,NEQ
         IF (ICNSTR(I) .EQ. 2) THEN
            IF (Y(I) .LE. ZERO) THEN
               IRET = I
               RETURN
            ENDIF
         ELSEIF (ICNSTR(I) .EQ. 1) THEN
            IF (Y(I) .LT. ZERO) THEN
               IRET = I
               RETURN
            ENDIF 
         ELSEIF (ICNSTR(I) .EQ. -1) THEN
            IF (Y(I) .GT. ZERO) THEN
               IRET = I
               RETURN
            ENDIF 
         ELSEIF (ICNSTR(I) .EQ. -2) THEN
            IF (Y(I) .GE. ZERO) THEN
               IRET = I
               RETURN
            ENDIF 
        ENDIF
 100  CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCNST0 ----------------------
      END
      SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDAWTS
C***REFER TO  DDASPK
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***END PROLOGUE  DDAWTS
C-----------------------------------------------------------------------
C     This routine sets the error weight vector,
C     WT, according to WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
C     I = 1 to NEQ.
C     RTOL and ATOL are scalars if IWT = 0,
C     and vectors if IWT = 1.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RTOL(*),ATOL(*),Y(*),WT(*)
      DIMENSION RPAR(*),IPAR(*)
      IF (IWT .EQ. 0) THEN
         RTOLI=RTOL(1)
         ATOLI=ATOL(1)
         DO I = 1, NEQ
            WT(I) = RTOLI*ABS(Y(I))+ATOLI
         END DO
      ELSE
         DO I = 1, NEQ
            RTOLI=RTOL(I)
            ATOLI=ATOL(I)
            WT(I)=RTOLI*ABS(Y(I))+ATOLI
         END DO
      END IF
      RETURN
C
C------END OF SUBROUTINE DDAWTS-----------------------------------------
      END
      SUBROUTINE DINVWT(NEQ,WT,IER)
C
C***BEGIN PROLOGUE  DINVWT
C***REFER TO  DDASPK
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   950125   (YYMMDD)
C***END PROLOGUE  DINVWT
C-----------------------------------------------------------------------
C     This routine checks the error weight vector WT, of length NEQ,
C     for components that are .le. 0, and if none are found, it
C     inverts the WT(I) in place.  This replaces division operations
C     with multiplications in all norm evaluations.
C     IER is returned as 0 if all WT(I) were found positive,
C     and the first I with WT(I) .le. 0.0 otherwise.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION WT(*)
C
      DO 10 I = 1,NEQ
        IF (WT(I) .LE. 0.0D0) GO TO 30
 10     CONTINUE
      DO 20 I = 1,NEQ
 20     WT(I) = 1.0D0/WT(I)
      IER = 0
      RETURN
C
 30   IER = I
      RETURN
C
C------END OF SUBROUTINE DINVWT-----------------------------------------
      END
      SUBROUTINE DDATRP(X,XOUT,YOUT,YPOUT,NEQ,NPHI,KOLD,PHI,PSI)
C
C***BEGIN PROLOGUE  DDATRP
C***REFER TO  DDASPK
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***END PROLOGUE  DDATRP
C
C-----------------------------------------------------------------------
C     The methods in routine DDSTP use polynomials
C     to approximate the solution.  DDATRP approximates the
C     solution and its derivative at time XOUT by evaluating
C     one of these polynomials, and its derivative, there.
C     Information defining this polynomial is passed from
C     DDSTP, so DDATRP cannot be used alone.
C
C     The parameters are
C
C     X     The current time in the integration.
C     XOUT  The time at which the solution is desired.
C     YOUT  The interpolated approximation to Y at XOUT.
C           (This is output.)
C     YPOUT The interpolated approximation to YPRIME at XOUT.
C           (This is output.)
C     NEQ   Number of equations.
C     KOLD  Order used on last successful step.
C     PHI   Array of scaled divided differences of Y.
C     PSI   Array of past stepsize history.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION YOUT(*),YPOUT(*)
      DIMENSION PHI(NPHI,*),PSI(*)
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0D0
      C=1.0D0
      D=0.0D0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
30       CONTINUE
      RETURN
C
C------END OF SUBROUTINE DDATRP-----------------------------------------
      END
      DOUBLE PRECISION FUNCTION DDWNRM(NEQ,V,RWT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDWNRM
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***END PROLOGUE  DDWNRM
C-----------------------------------------------------------------------
C     This function routine computes the weighted
C     root-mean-square norm of the vector of length
C     NEQ contained in the array V, with reciprocal weights
C     contained in the array RWT of length NEQ.
C        DDWNRM=SQRT((1/NEQ)*SUM(V(I)*RWT(I))**2)
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(*),RWT(*)
      DIMENSION RPAR(*),IPAR(*)
      DDWNRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)*RWT(I)) .GT. VMAX) VMAX = ABS(V(I)*RWT(I))
10    CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)*RWT(I))/VMAX)**2
      DDWNRM = VMAX*SQRT(SUM/NEQ)
30    CONTINUE
      RETURN
C
C------END OF FUNCTION DDWNRM-------------------------------------------
      END
      DOUBLE PRECISION FUNCTION DDWNRMA(NEQ,V,RWT,RPAR,IPAR,ID)
C
C***BEGIN PROLOGUE  DDWNRM
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   990111   (YYMMDD)
C***REVISION DATE  990126   (YYMMDD)
C***END PROLOGUE  DDWNRM
C-----------------------------------------------------------------------
C     This function routine computes the weighted
C     root-mean-square norm of the vector of length
C     NEQ contained in the array V, with reciprocal weights
C     contained in the array RWT of length NEQ.
C        DDWNRM=SQRT((1/NEQ)*SUM(V(I)*RWT(I))**2)
C     for the algebraic variables.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(*),RWT(*)
      DIMENSION RPAR(*),IPAR(*),ID(NEQ)
      DDWNRMA = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)*RWT(I)) .GT. VMAX) VMAX = ABS(V(I)*RWT(I))
10    CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      NEQA = 0
      DO I = 1,NEQ
         IF (ID(I) .LT. 0) THEN
            SUM = SUM + ((V(I)*RWT(I))/VMAX)**2
            NEQA = NEQA+1
         END IF
      END DO
      DDWNRMA = VMAX*SQRT(SUM/NEQA)
30    CONTINUE
      RETURN
C
C------END OF FUNCTION DDWNRMI2-----------------------------------------
      END
      SUBROUTINE DDASID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JACD,PSOL,H,WT,
     *  JSDUM,IDID,RPAR,IPAR,DUMSVR,DELTA,R,YIC,YPIC,DUMPWK,WM,IWM,CJ,
     *  UROUND,DUME,DUMS,DUMR,EPCON,RATEMX,STPTOL,JFDUM,
     *  ICNFLG,ICNSTR,IERNLS,
     *  ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES,A_RES,DUMMY,
     *  T_RES, LIADF, ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DDASID
C***REFER TO  DDASPK
C***DATE WRITTEN   940701   (YYMMDD)
C***REVISION DATE  950808   (YYMMDD)
C***REVISION DATE  951110   Removed unreachable block 390.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C
C     DDASID solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
C     the initial conditions.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     ICOPT     -- Initial condition option chosen (1 or 2).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied routine to evaluate the
C                  residual.  See RES description in DDASPK prologue.
C     JACD      -- External user-supplied routine to evaluate the
C                  Jacobian.  See JAC description for the case
C                  INFO(12) = 0 in the DDASPK prologue.
C     PDUM      -- Dummy argument.
C     H         -- Scaling factor for this initial condition calc.
C     WT        -- Vector of weights for error criterion.
C     JSDUM     -- Dummy argument.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     DUMSVR    -- Dummy argument.
C     DELTA     -- Work vector for NLS of length NEQ.
C     R         -- Work vector for NLS of length NEQ.
C     YIC,YPIC  -- Work vectors for NLS, each of length NEQ.
C     DUMPWK    -- Dummy argument.
C     WM,IWM    -- Real and integer arrays storing matrix information
C                  such as the matrix of partial derivatives,
C                  permutation vector, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1,3) or 0 (ICOPT = 2).
C     UROUND    -- Unit roundoff.
C     DUME      -- Dummy argument.
C     DUMS      -- Dummy argument.
C     DUMR      -- Dummy argument.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     JFDUM     -- Dummy argument.
C     STPTOL    -- Tolerance used in calculating the minimum lambda
C                  value allowed.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length 
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0   ==> nonlinear solver converged.
C                   1,2 ==> recoverable error inside nonlinear solver.
C                           1 => retry with current Y, YPRIME
C                           2 => retry with original Y, YPRIME
C                  -1   ==> unrecoverable error in nonlinear solver.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied routine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C
C     All variables with "DUM" in their names are dummy variables
C     which are not used in this routine.
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DDSEN, DMATD, DNSID
C
C***END PROLOGUE  DDASID
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),ID(*),WT(*),ICNSTR(*)
      DIMENSION DELTA(*),R(*),YIC(*),YPIC(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION ISENFO(*),ISENWK(*),SENWRK(*),SENPAR(*),ADI(*)
      EXTERNAL  RES, JACD, G_RES, A_RES, T_RES, PSOL
C
      PARAMETER (LNRE=12, LNJE=13, LMXNIT=32, LMXNJ=33, LNSE=22)
      PARAMETER (LMTYPE=4, LNY=37, LLCIWP=30)
C
C
C     Perform initializations.
C
      MXNIT = IWM(LMXNIT)
      MXNJ = IWM(LMXNJ)
      NY = IWM(LNY)
      NYMNQ = NY - ISENFO(10)
      IERNLS = 0
      NJ = 0
C
C     Initialize all error flags to zero.
C
      IERJ = 0
      IERNEW = 0
C
C     Looping point for updating the Jacobian.
C
300   CONTINUE
C
C     Call DDSEN(RES) to initialize DELTA for the state variables
C
      IRES = 0
C
C     Evaluation only for the state variable
      NSTATE = 1
      IF (ISENFO(6) .LT. -1) THEN
         NSTATE = -ISENFO(6)
         ISENFO(6) = -1
      END IF
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR, G_RES, A_RES, ADI)
      IF (IRES .LT. 0) THEN
         IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
         GO TO 370
      END IF
C
C     Reevaluate the iteration matrix, J = dG/dY + CJ*dG/dYPRIME,
C     where G(X,Y,YPRIME) = 0.
C
      NJ = NJ + 1
      IWM(LNJE)=IWM(LNJE)+1
      IRES = 0
      IF (IWM(LMTYPE) .EQ. 10) THEN
C
C     User input LU solver
C
         CALL JACD(RES, IRES, NYMNQ, X, Y, YPRIME, WT, DELTA,
     *        R, H, CJ, WM, IWM(IWM(LLCIWP)), IERJ, RPAR,  
     *        IPAR, SENPAR, ICOPT, ID)
      ELSE 
         IF (ICOPT .EQ. 2 .OR. ICOPT .EQ. 6) THEN
            CALL DMATD(NY,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,R,
     *           WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR,SENPAR,LIADF,
     *           ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, ADI)
         ELSE 
            CALL DMATID(NY,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,R,ID,ICOPT,
     *           WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR,SENPAR,LIADF,
     *           ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, ADI)
         END IF
      END IF
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
C         
      IF (IRES .LT. 0 .OR. IERJ .NE. 0) GO TO 370
C
C     Call the nonlinear Newton solver for up to MXNIT iterations.
C
      CALL DNSID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,
     *     WT,RPAR,IPAR,DUMSVR,DELTA,R,
     *     YIC,YPIC,DUMPWK,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,RATEMX,MXNIT,
     *     STPTOL,ICNFLG,ICNSTR,IERNEW,
     *     ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES,
     *     A_RES, T_RES, LIADF, ADI, INDEX0, TDIST)
C
      IF (IERNEW .EQ. 1 .AND. NJ .LT. MXNJ) THEN
C
C        MXNIT iterations were done, the convergence rate is < 1,
C        and the number of Jacobian evaluations is less than MXNJ.
C        Call DDSEN (RES), reevaluate the Jacobian, and try again.
C
         GO TO 300
      END IF
C
      IF (IERNEW .NE. 0) GO TO 380
C
      IF (ISENFO(10) .GT. 0) THEN
         DO I = 1, NSTATE
            IPS = (I-1)*NY + 1
            IRES = 3
            IF(ISENFO(6) .LT. 0) ISENFO(6) = -I
C
C     Evaluate the residual for the quadrature only 
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,
     *           DELTA(IPS),IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT(IPS),CNST,SENPAR,G_RES, A_RES, ADI) 
            IF (IRES .LT. 0) THEN
               ISENFO(6) = -NSTATE
               GOTO 370
            END IF
            DO ISEN = 1, ISENFO(1) + 1
               IISEN = IPS + (ISEN-1)*NY
               DO IQ = NYMNQ, NY-1
                  IIQ = IISEN + IQ
                  YPRIME(IIQ) = YPRIME(IIQ) - DELTA(IIQ)
               END DO
            END DO
         END DO
      END IF
      RETURN
C
C
C     Unsuccessful exits from nonlinear solver.
C     Compute IERNLS accordingly.
C
370   IERNLS = 2
      IF (IERJ .NE. 0) THEN
         IF (ICOPT.EQ.1.OR.ICOPT.EQ.5) THEN
            IERNLS = 3
         ELSE
            IDID = -8
            IERNLS = -1
         END IF
         RETURN
      END IF
      IF (IRES .LE. -2) THEN
         IDID = -11
         IERNLS = -1
      END IF
      RETURN      
C
380   IERNLS = MIN(IERNEW,2)
      RETURN
C
C------END OF SUBROUTINE DDASID-----------------------------------------
      END
      SUBROUTINE DNSID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,LUSOL,
     *     WT,RPAR,IPAR,DUMSVR,DELTA,R,
     *     YIC,YPIC,DUMPWK,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,RATEMX,MAXIT,
     *     STPTOL,ICNFLG,ICNSTR,IERNEW,
     *     ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES,
     *     A_RES, T_RES, LIADF, ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DNSID
C***REFER TO  DDASPK
C***DATE WRITTEN   940701   (YYMMDD)
C***REVISION DATE  950713   (YYMMDD)
C***REVISION DATE  990503   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSID solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME
C     in the initial conditions.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     ICOPT     -- Initial condition option chosen (1,2 or 3).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine to evaluate the
C                  residual.  See RES description in DDASPK prologue.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     DELTA     -- Residual vector on entry, and work vector of
C                  length NEQ for DNSID.
C     WM,IWM    -- Real and integer arrays storing matrix information
C                  such as the matrix of partial derivatives,
C                  permutation vector, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1,3) or 0 (ICOPT=2).
C     R         -- Array of length NEQ used as workspace by the 
C                  linesearch routine DLINSD.
C     YIC,YPIC  -- Work vectors for DLINSD, each of length NEQ.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     STPTOL    -- Tolerance used in calculating the minimum lambda
C                  value allowed.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length 
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> failed to converge, but RATE .le. RATEMX.
C                   2  ==> failed to converge, RATE .gt. RATEMX.
C                   3  ==> other recoverable error (IRES = -1, or
C                          linesearch failed).
C                   4  ==> state variables converged but sensitivity
C                          variables failed to converge. 
C                  -1  ==> unrecoverable error (IRES = -2).
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DSLVD, DDWNRM, DLINSD, DCOPY, DDSEN
C
C***END PROLOGUE  DNSID
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),R(*)
      DIMENSION ID(*),DELTA(*), YIC(*), YPIC(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION ISENFO(*), ISENWK(*),SENWRK(*),SENPAR(*)
      DIMENSION ICNSTR(*), ADI(*)
      EXTERNAL  RES, G_RES, A_RES, T_RES, LUSOL
C      CHARACTER MSG*80
C
      PARAMETER (LNNI=19, LLSOFF=35, LNY=37)
      PARAMETER (LNSE=22, LNRE=21, LLNIWP=28, LKPRIN=31,LLCIWP=30)
C
C
C     Initializations.  M is the Newton iteration counter.
C
      LSOFF = IWM(LLSOFF)
      M = 0
      RATE = 1.0D0
      RLX = 0.4D0
      NY = IWM(LNY)
      NYMNQ = NY - ISENFO(10)
C>>>>>>>>>>>>>>>>>>>>>>staggered method<<<<<<<<<<<<<<<<<<<<
C     
      NSTATE = 1
      LOCNNI = 0
      IF (ISENFO(6) .LT. -1) NSTATE = -ISENFO(6)
      DO I = 1, NSTATE
         IPS = (I-1)*NY+1
C     
C     For the adjoint method
         IF (ISENFO(6).LT.0) ISENFO(6) = -I
         IF (NSTATE.GT.1 .AND. I.GT.1) THEN
            IRES = 0
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,
     *           DELTA(IPS),IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT(IPS),CNST,SENPAR,G_RES, A_RES, ADI)
         END IF
C     
         M = 0
         IF (ICOPT .EQ. 5) THEN
C.....evaluate the state variable first for the second
c     stage of index-2 case
C     
            IRES = 0
            CALL RESIDX2(X,Y(IPS),YPRIME(IPS),CJ,DELTA(IPS),
     *           IRES,RPAR,IPAR,SENPAR,
     *           WT(IPS),NY,ID(NY+1),RES,T_RES,CNST,
     *           SENWRK,IWM(LIADF),ISENFO,A_RES,ADI)
         END IF
         CALL DSLVD (LUSOL, NYMNQ, DELTA(IPS), WM, IWM, RPAR, IPAR)
         DELNRM = DDWNRM(NYMNQ,DELTA(IPS),WT(IPS),RPAR,IPAR)
         FNRM = DELNRM
         IF (INDEX0 .EQ. 1) FNRM = TDIST*ABS(CJ)*FNRM
         IF (FNRM .LE. EPCON) GOTO 450
C     
C     Newton iteration loop for state variables
C     
 400     CONTINUE
         LOCNNI = LOCNNI + 1
C     
C     Call linesearch routine for global strategy and set RATE
C     
         OLDFNM = FNRM
         CALL DLINSD (
     *        NEQ, Y(IPS), X, YPRIME(IPS), CJ, DELTA(IPS), 
     *        DELNRM, WT(IPS), LSOFF,
     *        STPTOL, IRET, RES, LUSOL, IRES, WM, IWM, FNRM, ICOPT, ID,
     *        R, YIC, YPIC, ICNFLG, ICNSTR, RLX, RPAR, IPAR,
     *        ISENFO, SENWRK, ISENWK,CNST, SENPAR, G_RES, A_RES, 
     *        T_RES, 0, LIADF, ADI, INDEX0, TDIST)
C     
         RATE = FNRM/OLDFNM
C     
C     Check for error condition from linesearch.
         IF (IRET .NE. 0) THEN
            IWM(LNNI) = IWM(LNNI) + LOCNNI/NSTATE
            IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
            GO TO 490
         END IF
C     
C     Test for convergence of the iteration, and return or loop.
C     
         IF (FNRM .LE. EPCON .OR. ISENFO(6).LT.0) GOTO 450
C     
C     The iteration has not yet converged.  Update M.
C     Test whether the maximum number of iterations have been tried.
C     
         M = M + 1
         IF (M .GE. MAXIT) THEN
            IWM(LNNI) = IWM(LNNI) + LOCNNI/NSTATE
            IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
            GO TO 480
         END IF
C     
C     Copy the residual to DELTA and its norm to DELNRM, and loop for
C     another iteration.
C     
         CALL DCOPY (NYMNQ, R, 1, DELTA(IPS), 1)
         DELNRM = FNRM      
         GO TO 400
 450     CONTINUE               ! state variables converge
      END DO
      IWM(LNNI) = IWM(LNNI) + LOCNNI/NSTATE
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
C>>>>>>>>>>>>>>>>>>>>>>>>>>state variable is done<<<<<<<<<<<<<<<<<<<<<<
C
      IF (ISENFO(1) .EQ. 0 .OR. ISENFO(6) .LT. 0) RETURN 
C     IF (IWM(LKPRIN) .GE. 2) THEN
C     MSG = '------ IN ROUTINE DNSID-- STATE VARIABLES END'
C     CALL XERRWD(MSG, 50, 900, 0, 0, 0, 0, 0, 0.0D0,0.0D0)
C     ENDIF
      M = 0
      RATE = 1.0D0
      RLX = 0.4D0
C     
C     Reevaluate the residual again for the sensitivites and solution
      IF (ISENFO(2) .EQ. 1) THEN
C
C     For the forward finite difference method, DELTA should be updated
C     for calculating the sensitivities.
         IRES = 0
         CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES, ADI)
         CALL DCOPY(NYMNQ, DELTA, 1, R, 1) 
         IF (IRES .LT. 0) GOTO 490
      ELSE IF (ISENFO(2) .EQ. 5) THEN
C....................matrix times vector methods.....................
C     
C     evaluate the Jacobian and DF/DP
C     
         LJACI = 1
         LJACJ = LJACI + NY + 1
         LJAC = 1 + ISENFO(4)*NY
         CALL JRADFSP (
     1        NY, X, Y, YPRIME, DELTA, G_RES, CJ,  
     2        SENWRK(LJAC), ISENWK(LJACJ), ISENWK(LJACI), 
     3        IPAR, RPAR, SENPAR, IRES, SENWRK, ISENFO(4), 
     4        IWM(LIADF))
      END IF
      IRES = 1
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      IF (IRES .LT. 0) GOTO 490
      IF (ICOPT .EQ. 5) THEN
C.....evaluate the state variable first for the second
c     stage of index-2 case
C     
         IRES = 1
         CALL RESIDX2(X,Y,YPRIME,CJ,DELTA,
     *        IRES,RPAR,IPAR,SENPAR,
     *        WT,NY,ID(NY+1),RES,T_RES,CNST,
     *        SENWRK,IWM(LIADF),ISENFO,A_RES,ADI)
         IF (IRES .LT. 0) GOTO 490
      END IF
      IF (ISENFO(7) .EQ. 2) THEN
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>staggered direct method<<<<<<<<<<<<<<<<<<
C     
      END IF
C     
C     Compute a new step vector DELTA by back-substitution.
C     
      DO I = 1, ISENFO(1)
         II = I*NY + 1
         CALL DSLVD (LUSOL, NYMNQ, DELTA(II), WM, IWM, RPAR, IPAR)
      END DO
C     
C     Get norm of DELTA.  Return now if norm(DELTA) .le. EPCON.
C     
      DELNRM = 0.0D0
      DO I = 1, ISENFO(1)
         II = I*NY + 1     
         TNORM2 = DDWNRM(NYMNQ,DELTA(II),WT(II),RPAR,IPAR)
         IF (TNORM2 .GT. DELNRM) DELNRM = TNORM2
      END DO
      FNRM = DELNRM
      IF (INDEX0 .EQ. 1) FNRM = TDIST*ABS(CJ)*FNRM
      IF (FNRM .LE. EPCON) RETURN
C     
C     Newton iteration loop for sensitivity variables
C     
 460  CONTINUE
C     
C     Call linesearch routine for global strategy and set RATE
C     
      OLDFNM = FNRM
C     
      CALL DLINSD (
     *     NEQ, Y, X, YPRIME, CJ, DELTA, DELNRM, WT, LSOFF,
     *     STPTOL, IRET, RES, LUSOL, IRES, WM, IWM, FNRM, ICOPT, ID,
     *     R, YIC, YPIC, ICNFLG, ICNSTR, RLX, RPAR, IPAR,
     *     ISENFO, SENWRK, ISENWK, CNST, SENPAR, G_RES, A_RES, 
     *     T_RES, -1, LIADF, ADI, INDEX0, TDIST)
C     
      RATE = FNRM/OLDFNM
C     
C     Check for error condition from linesearch.
      IF (IRET .NE. 0) GO TO 490
C     
C     Test for convergence of the iteration, and return or loop.
C     
      IF (FNRM .LE. EPCON) RETURN
C     
C     The iteration has not yet converged.  Update M.
C     Test whether the maximum number of iterations have been tried.
C     
      M = M + 1
      IF (M .GE. MAXIT) GO TO 480
C     
C     Copy the residual to DELTA and its norm to DELNRM, and loop for
C     another iteration.
C     
      CALL DCOPY (NEQ-NY, R(NY+1), 1, DELTA(NY+1), 1)
      DELNRM = FNRM      
      GO TO 460         
C     
C     The maximum number of iterations was done.  Set IERNEW and return.
C     
 480  IF (RATE .LE. RATEMX) THEN
         IERNEW = 1
      ELSE
         IERNEW = 2
      ENDIF
      RETURN
C     
 490  IF (IRES .LE. -2) THEN
         IERNEW = -1
      ELSE
         IERNEW = 3
      ENDIF
      RETURN
C
C------END OF SUBROUTINE DNSID------------------------------------------
      END
      SUBROUTINE DLINSD (NEQ, Y, X, YPRIME, CJ, P, PNRM, WT, LSOFF,
     *     STPTOL, IRET, RES, LUSOL, IRES, WM, IWM,
     *     FNRM, ICOPT, ID, R, YNEW, YPNEW, ICNFLG,
     *     ICNSTR, RLX, RPAR, IPAR,
     *     ISENFO, SENWRK, ISENWK, CNST,
     *     SENPAR,G_RES, A_RES, T_RES,INDEX, LIADF, 
     *     ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DLINSD
C***REFER TO  DNSID
C***DATE WRITTEN   941025   (YYMMDD)
C***REVISION DATE  941215   (YYMMDD)
C***REVISION DATE  960129   Moved line RL = ONE to top block.
C***REVISION DATE  990503   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DLINSD uses a linesearch algorithm to calculate a new (Y,YPRIME)
C     pair (YNEW,YPNEW) such that 
C
C     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME) ,
C
C     where 0 < RL <= 1.  Here, f(y,y') is defined as
C
C      f(y,y') = (1/2)*norm( (J-inverse)*G(t,y,y') )**2 ,
C
C     where norm() is the weighted RMS vector norm, G is the DAE
C     system residual function, and J is the system iteration matrix
C     (Jacobian).
C
C     In addition to the parameters defined elsewhere, we have
C
C     P       -- Approximate Newton step used in backtracking.
C     PNRM    -- Weighted RMS norm of P.
C     LSOFF   -- Flag showing whether the linesearch algorithm is
C                to be invoked.  0 means do the linesearch, and
C                1 means turn off linesearch.
C     STPTOL  -- Tolerance used in calculating the minimum lambda
C                value allowed.
C     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations
C                in the proposed new approximate solution will be
C                checked for, and the maximum step length will be
C                adjusted accordingly.
C     ICNSTR  -- Integer array of length NEQ containing flags for
C                checking constraints.
C     RLX     -- Real scalar restricting update size in DCNSTR.
C     YNEW    -- Array of length NEQ used to hold the new Y in
C                performing the linesearch.
C     YPNEW   -- Array of length NEQ used to hold the new YPRIME in
C                performing the linesearch.
C     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW).
C     YPRIME  -- Array of length NEQ containing the new YPRIME 
C                (i.e.,=YPNEW).
C     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the
C                current (Y,YPRIME) on input and output.
C     R       -- Work array of length NEQ, containing the scaled 
C                residual (J-inverse)*G(t,y,y') on return.
C     IRET    -- Return flag.
C                IRET=0 means that a satisfactory (Y,YPRIME) was found.
C                IRET=1 means that the routine failed to find a new
C                       (Y,YPRIME) that was sufficiently distinct from
C                       the current (Y,YPRIME) pair.
C                IRET=2 means IRES .LT. 0 from RES.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C     INDEX   -- Indicator for the corrector method.
C                INDEX=0, for the state variables only;
C                INDEX=1, for both the state variables and 
C                             sensitivity variables
C                INDEX=-1, for the sensitivity variables only
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DFNRMD, DYYPNW, DCOPY
C
C***END PROLOGUE  DLINSD
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL  RES, G_RES, A_RES, T_RES, LUSOL
      DIMENSION Y(*), YPRIME(*), WT(*), R(*), ID(*)
      DIMENSION WM(*), IWM(*)
      DIMENSION YNEW(*), YPNEW(*), P(*), ICNSTR(*)
      DIMENSION RPAR(*), IPAR(*), ADI(*)
      DIMENSION ISENFO(*), ISENWK(*), SENWRK(*), SENPAR(*)
C      CHARACTER MSG*80
C
      PARAMETER (LNRE=12, LKPRIN=31, LNY=37)
C
      SAVE ALPHA, ONE, TWO
      DATA ALPHA/1.0D-4/, ONE/1.0D0/, TWO/2.0D0/
C
C      KPRIN=IWM(LKPRIN)
C
      F1NRM = (FNRM*FNRM)/TWO
      RATIO = ONE
C      IF (KPRIN .GE. 2) THEN
C        MSG = '------ IN ROUTINE DLINSD-- PNRM = (R1) )'
C        CALL XERRWD(MSG, 40, 901, 0, 0, 0, 0, 1, PNRM, 0.0D0)
C        ENDIF
      TAU = PNRM
      IVIO = 0
      RL = ONE
      NY = IWM(LNY)
C-----------------------------------------------------------------------
C Check for violations of the constraints, if any are imposed.
C If any violations are found, the step vector P is rescaled, and the 
C constraint check is repeated, until no violations are found.
C ****This is done only for state variables
C-----------------------------------------------------------------------
      IF (INDEX .EQ. -1) GOTO 20
      IF (ICNFLG .NE. 0) THEN
 10      CONTINUE
         CALL DYYPNW (NY,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
         CALL DCNSTR (NY, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
         IF (IRET .EQ. 1) THEN
            IVIO = 1
            RATIO1 = TAU/PNRM
            RATIO = RATIO*RATIO1
            DO I = 1,NY
               P(I) = P(I)*RATIO1
            END DO
            PNRM = TAU
C            IF (KPRIN .GE. 2) THEN
C               MSG = '------ CONSTRAINT VIOL., PNRM = (R1),INDEX = (I1)'
C               CALL XERRWD(MSG, 50, 902, 0, 1, IVAR, 0, 1, PNRM, 0.0D0)
C            ENDIF
            IF (PNRM .LE. STPTOL) THEN
               IRET = 1
               RETURN
            ENDIF
            GO TO 10
         ENDIF
      ENDIF
C
 20   CONTINUE
      SLPI = (-TWO*F1NRM)*RATIO
      RLMIN = STPTOL/PNRM
C      IF (LSOFF .EQ. 0 .AND. KPRIN .GE. 2) THEN
C         MSG = '------ MIN. LAMBDA = (R1)'
C         CALL XERRWD(MSG, 25, 903, 0, 0, 0, 0, 1, RLMIN, 0.0D0)
C      ENDIF
C-----------------------------------------------------------------------
C Begin iteration to find RL value satisfying alpha-condition.
C If RL becomes less than RLMIN, then terminate with IRET = 1.
C-----------------------------------------------------------------------
      IF (INDEX .EQ.-1) THEN
         CALL DCOPY (NY, Y, 1, YNEW, 1)
         CALL DCOPY (NY, YPRIME, 1, YPNEW, 1)
      END IF
 100  CONTINUE
      IF (INDEX .GE. 0) 
     *     CALL DYYPNW (NY,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
      IF (INDEX .NE. 0) THEN
         NP = ISENFO(1)
         DO I = 1, NP
            II = I*NY + 1
            CALL DYYPNW (NY,Y(II),YPRIME(II),CJ,RL,
     *           P(II),ICOPT,ID,YNEW(II),YPNEW(II))
         END DO
      END IF
      CALL DFNRMD (
     *     NY, YNEW, X, YPNEW, R, CJ, WT, RES, LUSOL, IRES,
     *     FNRMP, WM, IWM, RPAR, IPAR,
     *     ISENFO, SENWRK, ISENWK, CNST, SENPAR, G_RES, A_RES, INDEX,
     *     ICOPT, ID, LIADF, T_RES, ADI, INDEX0, TDIST)
      IF (IRES .LT. 0) THEN
         IRET = 2
         RETURN
      ENDIF
      IF (LSOFF .EQ. 1) GO TO 150
C
      F1NRMP = FNRMP*FNRMP/TWO
C      IF (KPRIN .GE. 2) THEN
C        MSG = '------ LAMBDA = (R1)'
C        CALL XERRWD(MSG, 20, 904, 0, 0, 0, 0, 1, RL, 0.0D0)
C        MSG = '------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)'
C        CALL XERRWD(MSG, 43, 905, 0, 0, 0, 0, 2, F1NRM, F1NRMP)
C      ENDIF
      IF (F1NRMP .GT. F1NRM + ALPHA*SLPI*RL) GO TO 200
C-----------------------------------------------------------------------
C Alpha-condition is satisfied, or linesearch is turned off.
C Copy YNEW,YPNEW to Y,YPRIME and return.
C-----------------------------------------------------------------------
 150  IRET = 0
      IF (INDEX .GE. 0) THEN
         CALL DCOPY (NY, YNEW, 1, Y, 1)
         CALL DCOPY (NY, YPNEW, 1, YPRIME, 1)
      END IF
      IF (INDEX .NE. 0) THEN
         CALL DCOPY (NEQ-NY, YNEW(NY+1), 1, Y(NY+1), 1)
         CALL DCOPY (NEQ-NY, YPNEW(NY+1), 1, YPRIME(NY+1), 1)
      END IF
      FNRM = FNRMP
C      IF (KPRIN .GE. 1) THEN
C        MSG = '------ LEAVING ROUTINE DLINSD, FNRM = (R1)'
C        CALL XERRWD(MSG, 42, 906, 0, 0, 0, 0, 1, FNRM, 0.0D0)
C        ENDIF
      RETURN
C-----------------------------------------------------------------------
C Alpha-condition not satisfied.  Perform backtrack to compute new RL
C value.  If no satisfactory YNEW,YPNEW can be found sufficiently 
C distinct from Y,YPRIME, then return IRET = 1.
C-----------------------------------------------------------------------
 200  CONTINUE
      IF (RL .LT. RLMIN) THEN
        IRET = 1
        RETURN
      ENDIF
C
      RL = RL/TWO
      GO TO 100
C
C----------------------- END OF SUBROUTINE DLINSD ----------------------
      END
      SUBROUTINE DFNRMD (NY, Y, X, YPRIME, R, CJ, WT, RES, LUSOL,IRES,
     *                   FNORM, WM, IWM, RPAR, IPAR,
     *                   ISENFO, SENWRK, ISENWK, CNST, 
     *                   SENPAR, G_RES, A_RES, INDEX, ICOPT, ID, LIADF,
     *                   T_RES, ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DFNRMD
C***REFER TO  DLINSD
C***DATE WRITTEN   941025   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DFNRMD calculates the scaled preconditioned norm of the nonlinear
C     function used in the nonlinear iteration for obtaining consistent
C     initial conditions.  Specifically, DFNRMD calculates the weighted
C     root-mean-square norm of the vector (J-inverse)*G(T,Y,YPRIME),
C     where J is the Jacobian matrix.
C
C     In addition to the parameters described in the calling program
C     DLINSD, the parameters represent
C
C     R      -- Array of length NEQ that contains
C               (J-inverse)*G(T,Y,YPRIME) on return.
C     FNORM  -- Scalar containing the weighted norm of R on return.
C     INDEX   -- Indicator for the corrector method.
C                INDEX=0, for the state variables only;
C                INDEX=1, for both the state variables and 
C                             sensitivity variables
C                INDEX=-1, for the sensitivity variables only
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DDSEN, DSLVD, DDWNRM
C
C***END PROLOGUE  DFNRMD
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RES, G_RES, A_RES, T_RES, LUSOL
      DIMENSION Y(*), YPRIME(*), WT(*), R(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*), ID(*), ADI(*)
      DIMENSION ISENFO(*), ISENWK(*), SENWRK(*),SENPAR(*)
      PARAMETER (LNRE=12, LNSE=22, LLCIWP=30, LLNIWP=28, LMTYPE=4)
C-----------------------------------------------------------------------
C     Call DDSEN (RES) routine.
C-----------------------------------------------------------------------
      IF (INDEX .NE. 0) THEN
         IRES = 1
         CALL DDSEN(X,Y,YPRIME,CJ,R,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR, G_RES, A_RES, ADI)
      ELSE 
         IRES = 0
         CALL DDSEN(X,Y,YPRIME,CJ,R,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR, G_RES, A_RES, ADI)
      END IF
      IF (ICOPT .EQ. 5) THEN
         CALL RESIDX2(X,Y,YPRIME,CJ,R,IRES,RPAR,IPAR,SENPAR,
     *        WT,NY,ID(NY+1),RES,T_RES,CNST,SENWRK,
     *        IWM(LIADF),ISENFO,A_RES,ADI)
      END IF
      IF (IRES .LT. 0) RETURN
C-----------------------------------------------------------------------
C     Apply inverse of Jacobian to vector R.
C-----------------------------------------------------------------------
      NYMNQ = NY - ISENFO(10)
      IF (INDEX .GE. 0) CALL DSLVD (LUSOL,NYMNQ, R, WM, IWM, RPAR, IPAR)
      IF (INDEX .NE. 0) THEN
         NP = ISENFO(1)
         DO I = 1, NP
            II = I*NY + 1
            CALL DSLVD (LUSOL, NYMNQ, R(II), WM, IWM, RPAR, IPAR)
         END DO
      END IF
C-----------------------------------------------------------------------
C     Calculate norm of R.
C-----------------------------------------------------------------------
      FNORM = 0.0D0
      IF (INDEX .GE. 0) THEN
         FNORM = DDWNRM(NYMNQ,R,WT,RPAR,IPAR)
      END IF
      IF (INDEX .NE. 0) THEN
         DO I = 1, NP
            II = I*NY + 1
            TNORM2 = DDWNRM(NYMNQ,R(II),WT(II),RPAR,IPAR)
            IF (TNORM2 .GT. FNORM) FNORM = TNORM2
         END DO
      END IF
      IF (INDEX0 .EQ. 1) FNORM = TDIST*ABS(CJ)*FNORM
      RETURN
C----------------------- END OF SUBROUTINE DFNRMD ----------------------
      END
      SUBROUTINE DNEDD(X,Y,YPRIME,NEQ,RES,JACD,PDUM,H,WT,
     *   JSTART,IDID,RPAR,IPAR,PHI,NPHI,GAMMA,DUMSVR,DELTA,E,
     *   WM,IWM,CJ,CJOLD,CJLAST,S,UROUND,DUME,DUMS,DUMR,
     *   EPCON,JCALC,JFDUM,KP1,NONNEG,NTYPE,IERNLS,
     *   CK,VT,ENORM,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, DUMMY,
     *   LIADF,ADI)
C
C***BEGIN PROLOGUE  DNEDD
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNEDD solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JACD      -- External user-supplied routine to evaluate the
C                  Jacobian.  See JAC description for the case
C                  INFO(12) = 0 in the DDASPK prologue.
C     PDUM      -- Dummy argument.
C     H         -- Appropriate step size for next step.
C     WT        -- Vector of weights for error criterion.
C     JSTART    -- Indicates first call to this routine.
C                  If JSTART = 0, then this is the first call,
C                  otherwise it is not.
C     IDID      -- Completion flag, output by DNEDD.
C                  See IDID description in DDASPK prologue.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     PHI       -- Array of divided differences used by
C                  DNEDD.  The length is NEQ*(K+1),where
C                  K is the maximum order.
C     GAMMA     -- Array used to predict Y and YPRIME.  The length
C                  is MAXORD+1 where MAXORD is the maximum order.
C     DUMSVR    -- Dummy argument.
C     DELTA     -- Work vector for NLS of length NEQ.
C     E         -- Error accumulation vector for NLS of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Parameter always proportional to 1/H.
C     CJOLD     -- Saves the value of CJ as of the last call to DMATD.
C                  Accounts for changes in CJ needed to
C                  decide whether to call DMATD.
C     CJLAST    -- Previous value of CJ.
C     S         -- A scalar determined by the approximate rate
C                  of convergence of the Newton iteration and used
C                  in the convergence test for the Newton iteration.
C
C                  If RATE is defined to be an estimate of the
C                  rate of convergence of the Newton iteration,
C                  then S = RATE/(1.D0-RATE).
C
C                  The closer RATE is to 0., the faster the Newton
C                  iteration is converging; the closer RATE is to 1.,
C                  the slower the Newton iteration is converging.
C
C                  On the first Newton iteration with an up-dated
C                  preconditioner S = 100.D0, Thus the initial
C                  RATE of convergence is approximately 1.
C
C                  S is preserved from call to call so that the rate
C                  estimate from a previous step can be applied to
C                  the current step.
C     UROUND    -- Unit roundoff.
C     DUME      -- Dummy argument.
C     DUMS      -- Dummy argument.
C     DUMR      -- Dummy argument.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     JCALC     -- Flag used to determine when to update
C                  the Jacobian matrix.  In general:
C
C                  JCALC = -1 ==> Call the DMATD routine to update
C                                 the Jacobian matrix.
C                  JCALC =  0 ==> Jacobian matrix is up-to-date.
C                  JCALC =  1 ==> Jacobian matrix is out-dated,
C                                 but DMATD will not be called unless
C                                 JCALC is set to -1.
C     JFDUM     -- Dummy argument.
C     KP1       -- The current order(K) + 1;  updated across calls.
C     NONNEG    -- Flag to determine nonnegativity constraints.
C     NTYPE     -- Identification code for the NLS routine.
C                   0  ==> modified Newton; direct solver.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0  ==> nonlinear solver converged.
C                   1  ==> recoverable error inside nonlinear solver.
C                  -1  ==> unrecoverable error inside nonlinear solver.
C                  -2  ==> error test failure for state variables in
C                          staggered corrector method
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C
C     All variables with "DUM" in their names are dummy variables
C     which are not used in this routine.
C
C     Following is a list and description of local variables which
C     may not have an obvious usage.  They are listed in roughly the
C     order they occur in this subroutine.
C
C     The following group of variables are passed as arguments to
C     the Newton iteration solver.  They are explained in greater detail
C     in DNSD:
C        TOLNEW, MULDEL, MAXIT, IERNEW
C
C     IERTYP -- Flag which tells whether this subroutine is correct.
C               0 ==> correct subroutine.
C               1 ==> incorrect subroutine.
C 
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DDWNRM, RES, DMATD, DNSD, DDSEN
C
C***END PROLOGUE  DNEDD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ISENFO(*),ISENWK(*), SENWRK(*),SENPAR(*)
      DIMENSION Y(*),YPRIME(*),WT(*),VT(*)
      DIMENSION DELTA(*),E(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION PHI(NPHI,*),GAMMA(*), ADI(*)
      EXTERNAL  RES, JACD, G_RES, A_RES, PDUM
C
      PARAMETER (LNST=11, LNRE=12, LNJE=13, LNSE=22, LNY=37)
      PARAMETER (LLCIWP=30,LMTYPE=4)
C
      SAVE MULDEL, MAXIT, XRATE
      DATA MULDEL/1/, MAXIT/4/, XRATE/0.25D0/
C
C     Verify that this is the correct subroutine.
C
      IERTYP = 0
      NSTATE = 1
      IF (NTYPE .NE. 0) THEN
         IERTYP = 1
         GO TO 380
      ENDIF
C
C     If this is the first step, perform initializations.
C
      IF (JSTART .EQ. 0) THEN
         CJOLD = CJ
         JCALC = -1
      ENDIF
C
C     Perform all other initializations.
C
      IERNLS = 0
      NY = IWM(LNY)
C
C     Decide whether new Jacobian is needed.
C
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) THEN
         S = 100.D0
      END IF
C
C-----------------------------------------------------------------------
C     Entry point for updating the Jacobian with current
C     stepsize.
C-----------------------------------------------------------------------
300   CONTINUE
C
C     Initialize all error flags to zero.
C
      IERJ = 0
      IERNEW = 0
C
C     Predict the solution and derivative and compute the tolerance
C     for the Newton iteration.
C
      DO 310 I=1,NPHI
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NPHI
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
C        
      NYMNQ = NY - ISENFO(10)
      PNORM = DDWNRM (NYMNQ,Y,WT,RPAR,IPAR)
      IF (ISENFO(6) .LT. 0) THEN
         NP = -ISENFO(6) - 1
      ELSE
         NP = ISENFO(1)
      END IF
      DO I = 1, NP
         II = I*NY + 1     
         TNORM2 = DDWNRM(NYMNQ,Y(II),WT(II),RPAR,IPAR)
         IF (TNORM2 .GT. PNORM) PNORM = TNORM2
      END DO
      TOLNEW = 100.D0*UROUND*PNORM
C     
C     Call DDSEN (RES) to initialize DELTA only for the state variables,
C     which might be used for Jacobian or preconditioner evaluation
C
      IRES = 0
      IF (ISENFO(6) .LT. 0) THEN
         NSTATE = -ISENFO(6)
         ISENFO(6) = -1
         ISENFO(1) = NSTATE
      END IF
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR, G_RES, A_RES, ADI)
      IF (IRES .LT. 0) GO TO 380
C
C     If indicated, reevaluate the iteration matrix 
C     J = dG/dY + CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
C     Set JCALC to 0 as an indicator that this has been done.
C
      IF(JCALC .EQ. -1 .OR. ISENFO(7) .EQ. 2) THEN
         IWM(LNJE)=IWM(LNJE)+1
         IF (JCALC .EQ. -1) THEN
            S = 100.D0
         ELSE
            IF (MOD(IWM(LNST), 20) .EQ. 0) THEN
               S = 100.D0
            END IF
         END IF
         JCALC=0
         IF(ISENFO(7).EQ.2.AND.IWM(LNSE).GT.0.AND.ISENFO(2).EQ.5)THEN
            LJACI = 1
            LJACJ = LJACI + NY + 1
            LJAC = 1 + ISENFO(4)*NY
            CALL SD_JAC_ONLY(
     *           NY, SENWRK(LJAC),ISENWK(LJACJ),ISENWK(LJACI),
     *           CJ, WM, IWM, IERJ)
         ELSE
            IRES = 0
            IF (IWM(LMTYPE) .EQ. 10) THEN
C
C     User input LU solver
C
               CALL JACD(RES, IRES, NYMNQ, X, Y, YPRIME, WT, DELTA,
     *              E, H, CJ, WM, IWM(IWM(LLCIWP)), IERJ,  
     *              RPAR, IPAR, SENPAR, 0)
            ELSE 
               CALL DMATD(NY,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,E,WM,IWM,
     *              RES,IRES,UROUND,JACD,RPAR,IPAR,SENPAR,LIADF,
     *              ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, ADI)
            END IF
         END IF
         CJOLD=CJ
         IF (IRES .LT. 0) GO TO 380
         IF(IERJ .NE. 0) GO TO 380
      ENDIF
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
C
C     Call the nonlinear Newton solver.
C
      IF (NPHI .EQ. 2*NEQ) THEN
         DO I = 1, NY
            DELTA(I) = DELTA(I) - CJ*DELTA(I+NEQ)
         END DO
      END IF
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      CALL DNSD(X,Y,YPRIME,NEQ,NPHI,RES,PDUM,WT,RPAR,IPAR,DUMSVR,
     *          DELTA,E,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,S,TEMP1,
     *          TOLNEW,MULDEL,MAXIT,IRES,IDUM,IERNEW,
     *          CK,VT,ENORM,NONNEG,
     *          ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES,LIADF,
     *          ADI)
      IF (ISENFO(6) .LT. 0) ISENFO(1) = 0
C
      IF (IERNEW .GT. 0 .AND. JCALC .NE. 0) THEN
C
C        The Newton iteration had a recoverable failure with an old
C        iteration matrix.  Retry the step with a new iteration matrix.
C
         JCALC = -1
         GO TO 300
      ENDIF
C
      IF (IERNEW .EQ. 0) THEN
C     Newton iterations converge
         IF (ISENFO(10) .GT. 0 .AND. ISENFO(1) .GT. 0) THEN
            IRES = 3
C
C     Evaluate the residual only for the sensitivity of quadratures
            CALL DDSEN(X,Y,YPRIME,CJ,
     *           DELTA,IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT,CNST,SENPAR,G_RES, A_RES, ADI) 
            IF (IRES .LT. 0) GOTO 380
            CJINV = 1.0D0/CJ
            DO ISEN = 1, ISENFO(1)
               IISEN = 1 + ISEN*NY
               DO IQ = NYMNQ, NY-1
                  IIQ = IISEN + IQ
                  YPRIME(IIQ) = YPRIME(IIQ) - DELTA(IIQ)
                  Y(IIQ) = Y(IIQ) - CJINV*DELTA(IIQ)
                  E(IIQ) = -CJINV*DELTA(IIQ)
               END DO
            END DO
         END IF
         GO TO 390
      END IF
C
C
C     Exits from nonlinear solver.
C     No convergence with current iteration
C     matrix, or singular iteration matrix.
C     Compute IERNLS and IDID accordingly.
C
380   CONTINUE
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
      IF (IRES .LE. -2 .OR. IERTYP .NE. 0) THEN
         IERNLS = -1
         IF (IRES .LE. -2) IDID = -11
         IF (IERTYP .NE. 0) IDID = -15
      ELSE IF (IERNEW .EQ. -2) THEN
         IERNLS = -2
      ELSE 
         IERNLS = 1
         IF (IRES .LT. 0) IDID = -10
         IF (IERJ .NE. 0) IDID = -8
      ENDIF
C
390   JCALC = 1
      RETURN
C
C------END OF SUBROUTINE DNEDD------------------------------------------
      END
      SUBROUTINE DNSD(X,Y,YPRIME,NEQ,NPHI,RES,LUSOL,WT,RPAR,IPAR,
     *   DUMSVR,DELTA,E,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,
     *   S,CONFAC,TOLNEW,MULDEL,MAXIT,IRES,IDUM,IERNEW,
     *   CK,VT,ENORM,NONNEG,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES,LIADF,ADI)
C
C***BEGIN PROLOGUE  DNSD
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  950126   (YYMMDD)
C***REVISION DATE  990221   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSD solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     PDUM      -- Dummy argument.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     DUMSVR    -- Dummy argument.
C     DELTA     -- Work vector for DNSD of length NEQ.
C     E         -- Error accumulation vector for DNSD of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Parameter always proportional to 1/H (step size).
C     DUMS      -- Dummy argument.
C     DUMR      -- Dummy argument.
C     DUME      -- Dummy argument.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     S         -- Used for error convergence tests.
C                  In the Newton iteration: S = RATE/(1 - RATE),
C                  where RATE is the estimated rate of convergence
C                  of the Newton iteration.
C                  The calling routine passes the initial value
C                  of S to the Newton iteration.
C     CONFAC    -- A residual scale factor to improve convergence.
C     TOLNEW    -- Tolerance on the norm of Newton correction in
C                  alternative Newton convergence test.
C     MULDEL    -- A flag indicating whether or not to multiply
C                  DELTA by CONFAC.
C                  0  ==> do not scale DELTA by CONFAC.
C                  1  ==> scale DELTA by CONFAC.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     IRES      -- Error flag returned from RES.  See RES description
C                  in DDASPK prologue.  If IRES = -1, then IERNEW
C                  will be set to 1.
C                  If IRES < -1, then IERNEW will be set to -1.
C     IDUM      -- Dummy argument.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> recoverable error inside Newton iteration.
C                  -1  ==> unrecoverable error inside Newton iteration.
C                  -2  ==> error test failure for the state variables in
C                          staggered corrector method.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C
C     All arguments with "DUM" in their names are dummy arguments
C     which are not used in this routine.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DSLVD, DDWNRM, RES, DDSEN
C
C***END PROLOGUE  DNSD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ISENFO(*), ISENWK(*),SENWRK(*),SENPAR(*)
      DIMENSION Y(*),YPRIME(*),WT(*),DELTA(*),E(*),VT(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*), ADI(*)
      EXTERNAL  RES,G_RES, A_RES, LUSOL
C
      PARAMETER (LNRE=12, LNNI=19,LNSE=22, LNY=37)
      PARAMETER (LLNIWP=28, LLCIWP=30)
C
C     Initialize Newton counter M and accumulation vector E. 
C
      M = 0
      NY = IWM(LNY)
      NYMNQ = NY - ISENFO(10)
      NP = ISENFO(1)
      DO 100 I=1,NPHI
 100     E(I)=0.0D0
C
C     Corrector loop.
C
C>>>>>>>>>>>>>>>>>>>>>>>>>staggered method<<<<<<<<<<<<<<<<<<<
C
C     First compute the solution Y, not including the sensitivities
C
      NSTATE = 1
      IF (ISENFO(6) .LT. -1) NSTATE = - ISENFO(6)
      LOCNNI = 0
      IPS = 1
      ENORM = 0.0D0
      DO ISTATE = 1, NSTATE
C     
C     For the adjoint method
         IF (ISENFO(6) .LT. 0) S = 100.0D0
         IF (ISENFO(6) .LT. 0) ISENFO(6) = -ISTATE
         IF (NSTATE .GT. 1 .AND. ISTATE .GT. 1) THEN
            IRES = 0
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,DELTA(IPS),
     *           IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT(IPS),
     *           CNST,SENPAR,G_RES, A_RES,ADI)
            IF (NPHI .EQ. 2*NEQ) THEN
               DO I = 0, NYMNQ-1
                  II = I + IPS
                  DELTA(II) = DELTA(II) - CJ*DELTA(II+NEQ)
               END DO
            END IF
            IF (IRES .LT. 0) GOTO 575
         END IF
 400     CONTINUE
C     
C     If necessary, multiply residual by convergence factor.
C     
         IF (MULDEL .EQ. 1) THEN
            DO I = 0,NYMNQ-1
               DELTA(I+IPS) = DELTA(I+IPS) * CONFAC
            end do
         ENDIF
C     
C     Compute a new iterate (back-substitution).
C     Store the correction in DELTA.
C     
         LOCNNI = LOCNNI + 1
         CALL DSLVD(LUSOL,NYMNQ,DELTA(IPS),WM,IWM,RPAR,IPAR)
C     
C     Update Y, E, and YPRIME.
C     
         DO I=0,NYMNQ-1
            II = I+IPS
            Y(II)=Y(II)-DELTA(II)
            E(II)=E(II)-DELTA(II)
            YPRIME(II)=YPRIME(II)-CJ*DELTA(II)
         END DO
         IF (NPHI .EQ. 2*NEQ) THEN
            DO I = 0, NYMNQ-1
	       II = IPS + NEQ + I
               E(II) = E(II) - DELTA(II)
	       Y(II) = Y(II) - DELTA(II)
	       YPRIME(II) = YPRIME(II) - CJ*DELTA(II)
            END DO
C
C     Calculate the vector-matrix product -- v(F_y)
            CALL FXPV(X,DELTA(IPS),CJ,DELTA(NEQ+IPS),IRES,RPAR,
     *		 IPAR,RES,NY,
     *           ISENFO,ISENWK,SENPAR,ADI, 0)
            DO I = 0, NYMNQ-1
               II = NEQ+I+IPS
               E(II) = E(II) - DELTA(II)
	       Y(II) = Y(II) - DELTA(II)
	       YPRIME(II) =YPRIME(II) - CJ*DELTA(II) 
            END DO
         END IF
C     
C     Test for convergence of the iteration.
C     
         DELNRM=DDWNRM(NYMNQ,DELTA(IPS),WT(IPS),RPAR,IPAR)
         IF (M .EQ. 0) THEN
            OLDNRM = DELNRM
            IF (DELNRM .LE. TOLNEW) GO TO 470
         ELSE
            RATE = (DELNRM/OLDNRM)**(1.0D0/M)
            IF (RATE .GT. 0.9D0) GO TO 575
            S = RATE/(1.0D0 - RATE)
         END IF
         IF (S*DELNRM .LE. EPCON) GO TO 470
C     
C     The corrector has not yet converged.
C     Update M and test whether the
C     maximum number of iterations have
C     been tried.
C     
         M=M+1
         IF(M.GE.MAXIT) GO TO 575
C     
C     Evaluate the residual,
C     and go back to do another iteration.
C     
         IRES = 0
         CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,DELTA(IPS),
     *        IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT(IPS),CNST,
     *        SENPAR,G_RES, A_RES,ADI)
         IF (NPHI .EQ. 2*NEQ) THEN
            DO I = 0, NYMNQ-1
               II = I + IPS
               DELTA(II) = DELTA(II) - CJ*DELTA(II+NEQ)
            END DO
         END IF
         IF (IRES .LT. 0) GO TO 575
         GO TO 400
C     
C     The iteration for solution has converged.
C     
 470     M = 0
C     
C     The Newton iteration has converged.  If nonnegativity of
C     solution is required, set the solution nonnegative, if the
C     perturbation to do it is small enough.  If the change is too
C     large, then consider the corrector iteration to have failed.
C     
         IF(NONNEG .EQ. 0) GO TO 475
C.............for state variables only
         DO I = 0,NYMNQ-1
            DELTA(I+IPS) = MIN(Y(I+IPS),0.0D0)
         END DO
         DELNRM = DDWNRM(NYMNQ,DELTA(IPS),WT(IPS),RPAR,IPAR)
         IF(DELNRM .GT. EPCON) THEN
            IERNEW = -1
            RETURN
         END IF
         DO I = 0,NYMNQ-1
            E(I+IPS) = E(I+IPS) - DELTA(I+IPS)
         END DO
 475     CONTINUE
C
C     Evaluate the residual for the quadrature only
         IF (ISENFO(10) .GT. 0) THEN
            IRES = 3
            ISENFO(1) = 0
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,
     *           DELTA(IPS),IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT(IPS),CNST,SENPAR,G_RES, A_RES, ADI)
            ISENFO(1) = NP
            IF (IRES .LT. 0) GOTO 575
            CJINV = 1.0D0/CJ
            DO IQ = NYMNQ, NY-1
               IIQ = IPS + IQ
               YPRIME(IIQ) = YPRIME(IIQ) - DELTA(IIQ)
               Y(IIQ) = Y(IIQ) - CJINV*DELTA(IIQ)
               E(IIQ) = -CJINV*DELTA(IIQ)
            END DO
         END IF
C     Error test for the state variable
C        
         ENORM1 = DDWNRM(NY,E(IPS),VT(IPS),RPAR,IPAR)
         IF (NPHI .EQ. 2*NEQ) THEN
            ENORM2 = DDWNRM(NY,E(IPS+NEQ),VT(IPS+NEQ),RPAR,IPAR)
            IF (ENORM2 .GT. ENORM1) ENORM1 = ENORM2
         END IF
         IF (ENORM1 .GT. ENORM) ENORM = ENORM1
         IF (CK*ENORM1 .GT. 1.0D0) THEN
            IERNEW = -2
            IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
            IWM(LNNI) =IWM(LNNI) + LOCNNI/NSTATE
            RETURN
         END IF
         IPS = IPS + NY
      END DO
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
      IWM(LNNI) =IWM(LNNI) + LOCNNI/NSTATE
      IF (NY .EQ. NEQ .OR. ISENFO(6) .LT. 0) RETURN ! no sensitivity
C>>>>>>>>>>>>>>>>>>>>>>>>>>>state variable finished<<<<<<<<<<<<<<<<<<<<
C
C     Reevaluate the residual again for the sensitivites and solution
      IF (ISENFO(2) .EQ. 1) THEN
C
C     For the forward finite difference method, DELTA is used
C     for evaluating the sensitivity equations.
         IRES = 0
         CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES, ADI)
      END IF
      IF (ISENFO(2) .EQ. 5) THEN
C....................matrix times vector methods.....................
C     
C     evaluate the Jacobian and DF/DP
C     
         LJACI = 1
         LJACJ = LJACI + NY + 1
         LJAC = 1 + ISENFO(4)*NY
         CALL JRADFSP (
     1        NY, X, Y, YPRIME, DELTA, G_RES, CJ,  
     2        SENWRK(LJAC), ISENWK(LJACJ), ISENWK(LJACI), 
     3        IPAR, RPAR, SENPAR, IRES, SENWRK, ISENFO(4),
     4        IWM(LIADF))
      END IF
      IRES = 1
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
C     
      IF (ISENFO(7) .EQ. 2) THEN
C     
C>>>>>>>>>>>>>>>>>>>>Staggered direct method<<<<<<<<<<<<<<<<<<<<<<<<<<
C     
C     
         DO I = 1, NP
            II = I*NY + 1
            CALL DSLVD(LUSOL,NYMNQ,DELTA(II),WM,IWM,RPAR,IPAR)
            DO J = 0, NYMNQ-1
               JII = II + J
               Y(JII)=Y(JII)-DELTA(JII)
               E(JII)=E(JII)-DELTA(JII)
               YPRIME(JII)=YPRIME(JII)-CJ*DELTA(JII)
            END DO
         END DO
         RETURN
      END IF
C     
C>>>>>>>>>>>>>>>>>>>>>>>>>staggered corrector method<<<<<<<<<<<<<<<<<
C     
      SS = S
 500  CONTINUE
C     
C     Compute a new iterate (back-substitution).
C     Store the correction in DELTA.
C     
      IF (MULDEL .EQ. 1) THEN
         DO J = NY+1, NEQ
            DELTA(J) = DELTA(J) * CONFAC
         END DO
      END IF
      DO I = 1, NP
         II = I*NY + 1
C     
C     If necessary, multiply residual by convergence factor.
C     
         CALL DSLVD(LUSOL,NYMNQ,DELTA(II),WM,IWM,RPAR,IPAR)
         DO J = 0, NYMNQ-1
            JII = II + J
            Y(JII)=Y(JII)-DELTA(JII)
            E(JII)=E(JII)-DELTA(JII)
            YPRIME(JII)=YPRIME(JII)-CJ*DELTA(JII)
         END DO
      END DO
      DELNRM = 0.0D0
      DO I = 1, NP
         II = I*NY + 1
         ENORMTMP = DDWNRM(NYMNQ,DELTA(II),WT(II),RPAR,IPAR)
         IF (ENORMTMP .GT. DELNRM) DELNRM = ENORMTMP
      END DO
C     
C     Test for convergence of the iteration.
C     
      IF (M .EQ. 0) THEN
         OLDNRM = DELNRM
         IF (DELNRM .LE. TOLNEW) GO TO 570
      ELSE
         RATE = (DELNRM/OLDNRM)**(1.0D0/M)
         IF (RATE .GT. 0.9D0) GO TO 580
         SS = RATE/(1.0D0 - RATE)
      ENDIF
      IF (SS*DELNRM .LE. EPCON) GO TO 570
C     
C     The corrector has not yet converged.
C     Update M and test whether the
C     maximum number of iterations have
C     been tried.
C     
      M=M+1
      IF(M.GE.MAXIT) GO TO 580
C     
C     Evaluate the residual,
C     and go back to do another iteration.
C     
      IRES = 1
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      IF (IRES .LT. 0) GO TO 580
      GO TO 500
C     
C     The iteration has converged.
 570  CONTINUE
      RETURN     
C     
C     The iteration has not converged.  Set IERNEW appropriately.
C     
 575  CONTINUE
      IWM(LNNI) =IWM(LNNI) + LOCNNI/NSTATE                 
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
 580  CONTINUE
      IF (IRES .LE. -2 ) THEN
         IERNEW = -1
      ELSE
         IERNEW = 1
      ENDIF
      RETURN
C
C------END OF SUBROUTINE DNSD-------------------------------------------
      END
      SUBROUTINE DMATID(NEQ,X,Y,YPRIME,DELTA,CJ,H,IER,EWT,E,ID,ICOPT,
     *     WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR,SENPAR,LIADF,
     *     ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES,ADI)
C
C***BEGIN PROLOGUE  DMATID
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD) (new LIPVT)
C***REVISION DATE  981020   (YYMMDD) (Add Adifor options)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine computes the iteration matrix
C     J = dG/dY+CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
C     for the computation of consistent initial conditions.
C     Here J is computed by:
C       the user-supplied routine JACD if IWM(LMTYPE) is 1 or 5, or
C       by numerical difference quotients if IWM(LMTYPE) is 2 or 6 or
C       by ADIFOR with SparsLinC if IWM(LMTYPE) is 3 or 7, or
C       by ADIFOR with seed matrix option if IWM(MTYPE) is 4 or 8.
C
C     The parameters have the following meanings.
C     X        = Independent variable.
C     Y        = Array containing predicted values.
C     YPRIME   = Array containing predicted derivatives.
C     DELTA    = Residual evaluated at (X,Y,YPRIME).
C                (Used only if IWM(MTYPE)=2 or 6).
C     CJ       = Scalar parameter defining iteration matrix.
C     H        = Current stepsize in integration.
C     IER      = Variable which is .NE. 0 if iteration matrix
C                is singular, and 0 otherwise.
C     EWT      = Vector of error weights for computing norms.
C     E        = Work space (temporary) of length NEQ.
C     WM       = Real work space for matrices.  On output
C                it contains the LU decomposition
C                of the iteration matrix.
C     ID       = Array of dimension NEQ, which indicates which variable 
C                is algebraic or fixed during the computation.
C     IWM      = Integer work space containing
C                matrix information.
C     RES      = External user-supplied subroutine
C                to evaluate the residual.  See RES description
C                in DDASPK prologue.
C     IRES     = Flag which is equal to zero if no illegal values
C                in RES, and less than zero otherwise.  (If IRES
C                is less than zero, the matrix was not completed).
C                In this case (if IRES .LT. 0), then IER = 0.
C     UROUND   = The unit roundoff error of the machine being used.
C     JACD     = Name of the external user-supplied routine
C                to evaluate the iteration matrix, or the name of 
C                ADIFOR generated routine.
C                See JAC description for the case INFO(12) = 0
C                in DDASPK prologue.
C     RPAR,IPAR= Real and integer parameter arrays that
C                are used for communication between the
C                calling program and external user routines.
C                They are not altered by DMATD.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   JACD, JDADFSP, JDADFSM, DGEFA, DGBFA, DDSEN
C
C***END PROLOGUE  DMATD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),DELTA(*),EWT(*),E(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*),SENPAR(*),ID(*)
      DIMENSION ISENFO(*), SENWRK(*), ISENWK(*), ADI(*)
      EXTERNAL  RES, JACD, G_RES, A_RES
C
      PARAMETER (LML=1,LMU=2,LMTYPE=4,LNRE=12,LNPD=36,LLCIWP=30)
      PARAMETER (LNSE=22)
C
      LIPVT = IWM(LLCIWP)
      IER = 0
      MTYPE=IWM(LMTYPE)
      IJAC = 1
      IF (ICOPT .EQ. 2 .OR. ICOPT.EQ.6) IJAC = 0
      IF (ICOPT .EQ. 5) IJAC = 2
      NEQMNQ = NEQ - ISENFO(10)
      GO TO (100,200,300,310,400,500,600,610),MTYPE
C
C
C     Dense user-supplied matrix.
C
 100  LENPD=IWM(LNPD)
      DO 110 I=1,LENPD
 110     WM(I)=0.0D0
      IF (ISENFO(6) .LT. 0) THEN
         CALL JACD(X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      ELSE 
         CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      END IF        
      GO TO 330
C
C
C     Dense finite-difference-generated matrix.
C
 200  IRES=0
      NROW=0
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQMNQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),
     *        ABS(1.D0/EWT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         IF (ICOPT .EQ. 1) THEN  
C     for ICOPT = 1
            IF (ID(I) .GT. 0) THEN
               YPRIME(I)=YPRIME(I)+CJ*DEL
            ELSE 
               Y(I)=Y(I)+DEL
            END IF
         ELSE IF (ICOPT .EQ. 5) THEN
C     for ICOPT = 5
            IF (ID(I) .GT. -2) THEN
               YPRIME(I)=YPRIME(I)+CJ*DEL
            ELSE
               Y(I) = Y(I) + DEL
            END IF           
         ELSE IF (ICOPT.EQ.3 .OR. ICOPT.EQ.4) THEN
C     for ICOPT = 3, 4
            IF (ID(I) .EQ. 1) THEN
               Y(I)=Y(I)+DEL
               YPRIME(I)=YPRIME(I)+CJ*DEL
            ELSE IF (ID(I) .EQ. 2) THEN
               YPRIME(I)=YPRIME(I)+CJ*DEL
            ELSE IF (ID(I).LT.0 .OR. ID(I).EQ.3) THEN
               Y(I)=Y(I)+DEL
            END IF
         END IF
         CALL DDSEN(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR,RES,NEQ,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES,ADI)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         IF (ICOPT .LT. 5) THEN
            DO L=1,NEQMNQ
               WM(NROW+L)=(E(L)-DELTA(L))*DELINV
            END DO
         ELSE
            DO L=1,NEQMNQ
               IF (ID(NEQ+L) .NE. 1) THEN
                  WM(NROW+L)=(E(L)-DELTA(L))*DELINV
               END IF
            END DO
            Y(I) = YSAVE + DEL
            CALL DDSEN(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR,RES,NEQ,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *           SENPAR,G_RES, A_RES,ADI)
            IF (IRES .LT. 0) RETURN
            DO L=1,NEQMNQ
               IF (ID(NEQ+L) .EQ. 1) THEN
                  WM(NROW+L)=(E(L)-DELTA(L))*DELINV
               END IF
            END DO
         END IF
         NROW=NROW+NEQMNQ
         Y(I)=YSAVE
         YPRIME(I)=YPSAVE
 210  CONTINUE
      GO TO 330
C
C     ADIFOR routine to compute the Jacobian
C     IWM(MTYPE)=3.
C
 300  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      END DO
      IF (ISENFO(6) .LT. 0) THEN
C
C...  Jacobian for Adjoint equations      
         CALL JIDADFSP(
     1        X, ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,0,IWM(LML),
     2        IWM(LMU),ISENWK(3),ICOPT,JACD,
     3        IRES,E,IWM(LIPVT+NEQ),
     4        IWM(LIADF),ISENFO(10), ISENFO(6))
      ELSE
         CALL JIDADFSP(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,0, IWM(LMU),
     2        IWM(LML),ID,ICOPT,JACD,IRES,E,IWM(LIPVT+NEQ),
     3        IWM(LIADF), ISENFO(10), ISENFO(6))
      END IF
      GO TO 330
C
C     ADIFOR routine to compute the Jacobian, 
C     seed matrix option is selected
C     IWM(MTYPE)=4.
C
 310  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      END DO
      IF (ISENFO(6) .LT. 0) THEN
         CALL JIDADFSM(
     1        X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,0, IWM(LML),
     2        IWM(LMU),ISENWK(3),ICOPT,JACD,IRES,
     3        E, WM(LENPD+1),
     4        ISENFO(10), ISENFO(6))
      ELSE
         CALL JIDADFSM(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,0, IWM(LMU),
     2        IWM(LML),ID,ICOPT,JACD,IRES,E, WM(LENPD+1),
     3        ISENFO(10), ISENFO(6))
      END IF
C
C
C     Do dense-matrix LU decomposition on J.
C
      
 330  CALL DGEFA(WM,NEQMNQ,NEQMNQ,IWM(LIPVT),IER)
      RETURN
C
C
C     Banded user-supplied matrix.
C
 400  LENPD=IWM(LNPD)
      DO 410 I=1,LENPD
 410     WM(I)=0.0D0
      IF (ISENFO(6) .LT. 0) THEN
         CALL JACD(X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      ELSE 
         CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      END IF        
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 650
C     
C     
C     Banded finite-difference-generated matrix.
C
 500  MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN0(MBAND,NEQMNQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQMNQ/MBAND)+1
      ISAVE=IWM(LNPD)
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO N=J,NEQMNQ,MBAND
            K= (N-J)/MBAND + 1
            WM(ISAVE+K)=Y(N)
            WM(IPSAVE+K)=YPRIME(N)
            DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *           ABS(1.D0/EWT(N)))
            DEL=SIGN(DEL,H*YPRIME(N))
            DEL=(Y(N)+DEL)-Y(N)
            IF (ICOPT .EQ. 1) THEN  
C     for ICOPT = 1
               IF (ID(N) .GT. 0) THEN
                  YPRIME(N)=YPRIME(N)+CJ*DEL
               ELSE 
                  Y(N)=Y(N)+DEL
               END IF
            ELSE IF (ICOPT .EQ. 5) THEN
C     for ICOPT = 5
               IF (ID(N) .GT. -2) THEN
                  YPRIME(N)=YPRIME(N)+CJ*DEL
               ELSE 
                  Y(N)=Y(N)+DEL
               END IF               
            ELSE 
C     for ICOPT = 3, 4
               IF (ID(N) .EQ. 1) THEN
                  Y(N)=Y(N)+DEL
                  YPRIME(N)=YPRIME(N)+CJ*DEL
               ELSE IF (ID(N) .EQ. 2) THEN
                  YPRIME(N)=YPRIME(N)+CJ*DEL
               ELSE IF (ID(N).LT.0 .OR. ID(N).EQ.3) THEN
                  Y(N)=Y(N)+DEL
               END IF
            END IF
         END DO
         CALL DDSEN(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR,RES,NEQ,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES,ADI)
         IF (IRES .LT. 0) RETURN
         DO 530 N=J,NEQMNQ,MBAND
            K= (N-J)/MBAND + 1
            Y(N)=WM(ISAVE+K)
            YPRIME(N)=WM(IPSAVE+K)
            DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *           ABS(1.D0/EWT(N)))
            DEL=SIGN(DEL,H*YPRIME(N))
            DEL=(Y(N)+DEL)-Y(N)
            DELINV=1.0D0/DEL
            I1=MAX0(1,(N-IWM(LMU)))
            I2=MIN0(NEQMNQ,(N+IWM(LML)))
            II=N*MEB1-IWM(LML)
            IF (ICOPT .LT. 5) THEN
               DO I=I1,I2
                  WM(II+I)=(E(I)-DELTA(I))*DELINV
               END DO
            ELSE
               DO I=I1,I2
                  IF (ID(NEQ+I) .NE. 1) THEN
                     WM(II+I)=(E(I)-DELTA(I))*DELINV
                  END IF
               END DO
            END IF
 530     CONTINUE
         IF (ICOPT .EQ. 5) THEN
            DO N=J,NEQMNQ,MBAND
               K= (N-J)/MBAND + 1
               WM(ISAVE+K)=Y(N)
               DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *              ABS(1.D0/EWT(N)))
               DEL=SIGN(DEL,H*YPRIME(N))
               DEL=(Y(N)+DEL)-Y(N)
               Y(N)=Y(N)+DEL
            END DO
            CALL DDSEN(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR,RES,NEQ,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *           SENPAR,G_RES, A_RES,ADI)
            IF (IRES .LT. 0) RETURN
            DO N=J,NEQMNQ,MBAND
               K= (N-J)/MBAND + 1
               Y(N)=WM(ISAVE+K)
               DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *              ABS(1.D0/EWT(N)))
               DEL=SIGN(DEL,H*YPRIME(N))
               DEL=(Y(N)+DEL)-Y(N)
               DELINV=1.0D0/DEL
               I1=MAX0(1,(N-IWM(LMU)))
               I2=MIN0(NEQMNQ,(N+IWM(LML)))
               II=N*MEB1-IWM(LML)
               DO I=I1,I2
                  IF (ID(NEQ+I) .EQ. 1) THEN
                     WM(II+I)=(E(I)-DELTA(I))*DELINV
                  END IF
               END DO
            END DO
         END IF
 540  CONTINUE
      GO TO 650
 600  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      END DO
      IF (ISENFO(6) .LT. 0) THEN
         CALL JIDADFSP(
     1        X, ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,1,IWM(LML),
     2        IWM(LMU),ISENWK(3),ICOPT,JACD,IRES,
     3        E,IWM(LIPVT+NEQ),
     4        IWM(LIADF),ISENFO(10), ISENFO(6))
      ELSE
         CALL JIDADFSP(
     1        X, Y, YPRIME, WM,CJ,RPAR,IPAR,SENPAR,NEQ,1,IWM(LMU),
     2        IWM(LML),ID,ICOPT,JACD,IRES,E,IWM(LIPVT+NEQ),
     3        IWM(LIADF),ISENFO(10), ISENFO(6))
      END IF
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 650
C
C     ADIFOR routine to compute the Jacobian, 
C     seed matrix option is selected
C     IWM(MTYPE)=8.
C
 610  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      end do
      IF (ISENFO(6) .LT. 0) THEN
         CALL JIDADFSM(
     1        X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,1, IWM(LML),
     2        IWM(LMU),ISENWK(3),ICOPT,JACD,IRES,
     3        E, WM(LENPD+1),
     4        ISENFO(10), ISENFO(6))
      ELSE
         CALL JIDADFSM(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,1, IWM(LMU),
     2        IWM(LML),ID,ICOPT,JACD,IRES,E, WM(LENPD+1),
     3        ISENFO(10), ISENFO(6))
      END IF
      MEBAND=2*IWM(LML)+IWM(LMU)+1
C
C
C     Do LU decomposition of banded J.
C
 650  CALL DGBFA (WM,MEBAND,NEQMNQ,IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C
C------END OF SUBROUTINE DMATID------------------------------------------
      END
      SUBROUTINE DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IER,EWT,E,
     *     WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR,SENPAR,LIADF,
     *     ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES,ADI)
C
C***BEGIN PROLOGUE  DMATD
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD) (new LIPVT)
C***REVISION DATE  981020   (YYMMDD) (Add Adifor options)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine computes the iteration matrix
C     J = dG/dY+CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
C     Here J is computed by:
C       the user-supplied routine JACD if IWM(LMTYPE) is 1 or 5, or
C       by numerical difference quotients if IWM(LMTYPE) is 2 or 6 or
C       by ADIFOR with SparsLinC if IWM(MTYPE) is 3 or 7, or
C       by ADIFOR with seed matrix option if IWM(LMTYPE) is 4 or 8.
C
C     The parameters have the following meanings.
C     X        = Independent variable.
C     Y        = Array containing predicted values.
C     YPRIME   = Array containing predicted derivatives.
C     DELTA    = Residual evaluated at (X,Y,YPRIME).
C                (Used only if IWM(MTYPE)=2 or 6).
C     CJ       = Scalar parameter defining iteration matrix.
C     H        = Current stepsize in integration.
C     IER      = Variable which is .NE. 0 if iteration matrix
C                is singular, and 0 otherwise.
C     EWT      = Vector of error weights for computing norms.
C     E        = Work space (temporary) of length NEQ.
C     WM       = Real work space for matrices.  On output
C                it contains the LU decomposition
C                of the iteration matrix.
C     IWM      = Integer work space containing
C                matrix information.
C     RES      = External user-supplied subroutine
C                to evaluate the residual.  See RES description
C                in DDASPK prologue.
C     IRES     = Flag which is equal to zero if no illegal values
C                in RES, and less than zero otherwise.  (If IRES
C                is less than zero, the matrix was not completed).
C                In this case (if IRES .LT. 0), then IER = 0.
C     UROUND   = The unit roundoff error of the machine being used.
C     JACD     = Name of the external user-supplied routine
C                to evaluate the iteration matrix, or the name of 
C                ADIFOR generated routine.
C                See JAC description for the case INFO(12) = 0
C                in DDASPK prologue.
C     RPAR,IPAR= Real and integer parameter arrays that
C                are used for communication between the
C                calling program and external user routines.
C                They are not altered by DMATD.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   JACD, DDSEN, JDADFSP, JDADFSM, DGEFA, DGBFA
C
C***END PROLOGUE  DMATD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),DELTA(*),EWT(*),E(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*),SENPAR(*)
      DIMENSION ISENFO(*), SENWRK(*), ISENWK(*), ADI(*)
      EXTERNAL  RES, JACD, G_RES, A_RES
C
      PARAMETER (LML=1,LMU=2,LMTYPE=4,LNRE=12,LNPD=36,LLCIWP=30)
      PARAMETER (LNSE=22)
C
      LIPVT = IWM(LLCIWP)
      IER = 0
      MTYPE=IWM(LMTYPE)
      NEQMNQ = NEQ - ISENFO(10)
      GO TO (100,200,300,310,400,500,600,610),MTYPE
C
C
C     Dense user-supplied matrix.
C
 100  LENPD=IWM(LNPD)
      DO 110 I=1,LENPD
 110     WM(I)=0.0D0
      IJAC = 0
      IF (ISENFO(6) .LT. 0) THEN
         CALL JACD(X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      ELSE 
         CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      END IF        
      GO TO 330
C
C
C     Dense finite-difference-generated matrix.
C
 200  IRES=0
      NROW=0
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQMNQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(1.D0/EWT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL DDSEN(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR,RES,NEQ,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES, ADI)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         IF (ISENFO(6).GE.0 .OR. (ISENFO(6).LT.0.AND.ISENFO(12).NE.2.
     *        AND.ISENFO(12).NE.6)) THEN
            DO 220 L=1,NEQMNQ
               WM(NROW+L)=(E(L)-DELTA(L))*DELINV
 220        CONTINUE
         ELSE
C     For the adjoint equation (Ay)' + By = 0 case
            NEQT = NEQ*ISENFO(1)
            DO L = 1, NEQMNQ
               WM(NROW+L)= (E(L)-DELTA(L))*DELINV - 
     *              CJ*(E(L+NEQT)-DELTA(L+NEQT))*DELINV
            END DO
         END IF
         NROW=NROW+NEQMNQ
         Y(I)=YSAVE
         YPRIME(I)=YPSAVE
 210  CONTINUE
      GO TO 330
C
C     ADIFOR routine to compute the Jacobian
C     IWM(MTYPE)=3.
C
 300  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      END DO
      IF (ISENFO(6) .LT. 0) THEN
C
C...  Jacobian for Adjoint equations
         CALL JDADFSP(
     1        X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,0,IWM(LML),
     2        IWM(LMU), JACD, IRES, E, IWM(LIPVT+NEQ),
     3        IWM(LIADF), ISENFO(10), ISENFO(6), ISENFO(7), ISENWK(3))
      ELSE
         CALL JDADFSP(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,0, IWM(LMU),
     2        IWM(LML), JACD, IRES, E, IWM(LIPVT+NEQ),
     3        IWM(LIADF), ISENFO(10), ISENFO(6), ISENFO(7), ISENWK)
      END IF
      GO TO 330
C
C     ADIFOR routine to compute the Jacobian, 
C     seed matrix option is selected
C     IWM(MTYPE)=4.
C
 310  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      END DO
      IF (ISENFO(6) .LT. 0) THEN
C
C...  Jacobian for Adjoint equations
         CALL JDADFSM(
     1        X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,0,IWM(LML),
     2        IWM(LMU), JACD, IRES, E, WM(LENPD+1),
     3        ISENFO(10), ISENFO(6), ISENFO(7), ISENWK(3))
      ELSE
         CALL JDADFSM(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,0, IWM(LMU),
     2        IWM(LML), JACD, IRES, E, WM(LENPD+1),
     3        ISENFO(10), ISENFO(6), ISENFO(7), ISENWK)
      END IF
C
C
C     Do dense-matrix LU decomposition on J.
C
 330  CALL DGEFA(WM,NEQMNQ,NEQMNQ,IWM(LIPVT),IER)
      RETURN
C
C
C     Banded user-supplied matrix.
C
 400  LENPD=IWM(LNPD)
      DO 410 I=1,LENPD
 410     WM(I)=0.0D0
      IJAC = 0
      IF (ISENFO(6) .LT. 0) THEN
         CALL JACD(X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      ELSE 
         CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR,SENPAR,IJAC)
      END IF        
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 650
C     
C     
C     Banded finite-difference-generated matrix.
C
 500  MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN0(MBAND,NEQMNQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQMNQ/MBAND)+1
      ISAVE=IWM(LNPD)
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
        DO 510 N=J,NEQMNQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(1.D0/EWT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
 510      YPRIME(N)=YPRIME(N)+CJ*DEL
         CALL DDSEN(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR,RES,NEQ,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES,ADI)
        IF (IRES .LT. 0) RETURN
        DO 530 N=J,NEQMNQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(1.D0/EWT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX0(1,(N-IWM(LMU)))
          I2=MIN0(NEQMNQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)
         IF (ISENFO(6).GE.0 .OR. (ISENFO(6).LT.0.AND.ISENFO(12).NE.2.
     *        AND.ISENFO(12).NE.6)) THEN
            DO 520 I=I1,I2
               WM(II+I)=(E(I)-DELTA(I))*DELINV
 520        CONTINUE
         ELSE
C     For the adjoint equation (Ay)' + By = 0 case
            NEQT = NEQ*ISENFO(1)
            DO I = I1, I2
               WM(II+I)=(E(I)-DELTA(I))*DELINV -
     *              CJ*(E(I+NEQT)-DELTA(I+NEQT))*DELINV
            END DO
         END IF
 530  CONTINUE
 540   CONTINUE
      GO TO 650
 600  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      end do
      IF (ISENFO(6) .LT. 0) THEN
C
C...  Jacobian for Adjoint equations
         CALL JDADFSP(
     1        X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,1,IWM(LML),
     2        IWM(LMU),JACD, IRES, E, IWM(LIPVT+NEQ),
     3        IWM(LIADF),ISENFO(10), ISENFO(6), ISENFO(7), ISENWK(3))
      ELSE
         CALL JDADFSP(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,1,IWM(LMU),
     2        IWM(LML),JACD, IRES, E, IWM(LIPVT+NEQ),
     3        IWM(LIADF),ISENFO(10), ISENFO(6), ISENFO(7), ISENWK)
      END IF
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 650
C
C     ADIFOR routine to compute the Jacobian, 
C     seed matrix option is selected
C     IWM(MTYPE)=8.
C
 610  CONTINUE
      LENPD=IWM(LNPD)
      DO I=1,LENPD
         WM(I)=0.0D0
      end do
      IF (ISENFO(6) .LT. 0) THEN
C
C...  Jacobian for Adjoint equations
         CALL JDADFSM(
     1        X,ADI,ADI(NEQ+1),WM,CJ,RPAR,IPAR,SENPAR,NEQ,1, IWM(LML),
     2        IWM(LMU), JACD, IRES, E, WM(LENPD+1),
     3        ISENFO(10), ISENFO(6), ISENFO(7), ISENWK(3))
      ELSE
         CALL JDADFSM(
     1        X, Y, YPRIME,WM,CJ,RPAR,IPAR,SENPAR,NEQ,1, IWM(LMU),
     2        IWM(LML), JACD, IRES, E, WM(LENPD+1),
     3        ISENFO(10), ISENFO(6), ISENFO(7), ISENWK)
      END IF
      MEBAND=2*IWM(LML)+IWM(LMU)+1
C
C
C     Do LU decomposition of banded J.
C
 650  CALL DGBFA (WM,MEBAND,NEQMNQ,IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C
C------END OF SUBROUTINE DMATD------------------------------------------
      END
      SUBROUTINE DSLVD(LUSOL,NEQ,DELTA,WM,IWM,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DSLVD
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD) (new LIPVT)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine manages the solution of the linear
C     system arising in the Newton iteration.
C     Real matrix information and real temporary storage
C     is stored in the array WM.
C     Integer matrix information is stored in the array IWM.
C     For a dense matrix, the LINPACK routine DGESL is called.
C     For a banded matrix, the LINPACK routine DGBSL is called.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DGESL, DGBSL
C
C***END PROLOGUE  DSLVD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DELTA(*),WM(*),IWM(*),RPAR(*),IPAR(*)
      EXTERNAL LUSOL
C
      PARAMETER (LML=1, LMU=2, LMTYPE=4, LLCIWP=30)
C
      LIPVT = IWM(LLCIWP)
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,100,100,400,400,400,400,100,500),MTYPE
C
C     Dense matrix.
C
100   CONTINUE
      CALL DGESL(WM,NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     Banded matrix.
C
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL DGBSL(WM,MEBAND,NEQ,IWM(LML),
     *     IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
 500  CALL LUSOL(NEQ,WM,IWM(LIPVT),DELTA,IELS,RPAR,IPAR)
      RETURN
C
C------END OF SUBROUTINE DSLVD------------------------------------------
      END
      SUBROUTINE DDASIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JACK,PSOL,H,WT,
     *   JSKIP,IDID,RPAR,IPAR,SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,
     *   UROUND,EPLI,SQRTN,RSQRTN,EPCON,RATEMX,STPTOL,JFLG,
     *   ICNFLG,ICNSTR,IERNLS,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES,
     *   T_RES, LIADF,ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DDASIK
C***REFER TO  DDASPK
C***DATE WRITTEN   941026   (YYMMDD)
C***REVISION DATE  950808   (YYMMDD)
C***REVISION DATE  951110   Removed unreachable block 390.
C***REVISION DATE  990503   (Added sensitivity analysis)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C
C     DDASIK solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
C     the initial conditions.
C
C     An initial value for Y and initial guess for YPRIME are input.
C
C     The method used is a Newton scheme with Krylov iteration and a
C     linesearch algorithm.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector at x.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of equations to be integrated.
C     ICOPT     -- Initial condition option chosen (1,2 or 3).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JACK     --  External user-supplied routine to update
C                  the preconditioner.  (This is optional).
C                  See JAC description for the case
C                  INFO(12) = 1 in the DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning.
C                  (This is optional).  See explanation inside DDASPK.
C     H         -- Scaling factor for this initial condition calc.
C     WT        -- Vector of weights for error criterion.
C     JSKIP     -- input flag to signal if initial JAC call is to be
C                  skipped.  1 => skip the call, 0 => do not skip call.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     SAVR      -- Work vector for DDASIK of length NEQ.
C     DELTA     -- Work vector for DDASIK of length NEQ.
C     R         -- Work vector for DDASIK of length NEQ.
C     YIC,YPIC  -- Work vectors for DDASIK, each of length NEQ.
C     PWK       -- Work vector for DDASIK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information for linear system
C                  solvers, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1,3) or 0 (ICOPT = 2).
C     UROUND    -- Unit roundoff.
C     EPLI      -- convergence test constant.
C                  See DDASPK prologue for more details.
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     JFLG      -- Flag showing whether a Jacobian routine is supplied.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length 
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0   ==> nonlinear solver converged.
C                   1,2 ==> recoverable error inside nonlinear solver.
C                           1 => retry with current Y, YPRIME
C                           2 => retry with original Y, YPRIME
C                  -1   ==> unrecoverable error in nonlinear solver.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, DDSEN, JACK, DNSIK, DCOPY
C
C***END PROLOGUE  DDASIK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ISENFO(*), ISENWK(*),SENWRK(*),SENPAR(*)
      DIMENSION Y(*),YPRIME(*),ID(*),WT(*),ICNSTR(*)
      DIMENSION SAVR(*),DELTA(*),R(*),YIC(*),YPIC(*),PWK(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*),ADI(*)
      EXTERNAL RES, JACK, PSOL, G_RES, A_RES, K_RES, T_RES
C
      PARAMETER (LNRE=12, LNJE=13, LLOCWP=29, LLCIWP=30, LNSE=22)
      PARAMETER (LMXNIT=32, LMXNJ=33, LNY=37)
C
C
C     Perform initializations.
C
      LWP = IWM(LLOCWP)
      LIWP = IWM(LLCIWP)
      MXNIT = IWM(LMXNIT)
      MXNJ = IWM(LMXNJ)
      IERNLS = 0
      NJ = 0
      EPLIN = EPLI*EPCON
      NY = IWM(LNY)
      NYMNQ = NY - ISENFO(10)
C
C     Looping point for updating the preconditioner.
C
 300  CONTINUE
C
C     Call DDSEN (RES) to initialize DELTA.
C
      IRES = 0
C
C     Evaluation only for the state variables.
      NSTATE = 1
      IF (ISENFO(6) .LT. -1) THEN
         NSTATE = -ISENFO(6)
         ISENFO(6) = -1
      END IF
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR, G_RES, A_RES, ADI)
      IF (IRES .LT. 0) THEN
         IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
         GO TO 370
      END IF
C
C     Initialize all error flags to zero.
C
      IERPJ = 0
      IRES = 0
      IERNEW = 0
C
C     If a Jacobian routine was supplied, call it.
C
      IF (JFLG .EQ. 1 .AND. JSKIP .EQ. 0) THEN
         NJ = NJ + 1
         IWM(LNJE)=IWM(LNJE)+1
         IF (ISENFO(6) .LT. 0) THEN
            CALL JACKAD (JACK,
     *           RES, IRES, NYMNQ, X, Y, YPRIME, WT, DELTA, R,H, CJ,
     *           WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR,SENPAR, ADI,
     *           ISENFO,SENWRK,ISENWK,G_RES, A_RES, ICOPT, ID)
         ELSE
            CALL JACK (
     *           RES, IRES, NYMNQ, X, Y, YPRIME, WT, DELTA, R, H,CJ,
     *           WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR,SENPAR,ICOPT,ID)
         END IF            
         IF (IRES .LT. 0 .OR. IERPJ .NE. 0) GO TO 370
      ENDIF
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
      JSKIP = 0
C
C     Call the nonlinear Newton solver for up to MXNIT iterations.
C
      CALL DNSIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,WT,RPAR,IPAR,
     *   SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,SQRTN,RSQRTN,
     *   EPLIN,EPCON,RATEMX,MXNIT,STPTOL,ICNFLG,ICNSTR,IERNEW,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES,
     *   T_RES, LIADF, ADI, INDEX0, TDIST)
C
      IF (IERNEW .EQ. 1 .AND. NJ .LT. MXNJ .AND. JFLG .EQ. 1) THEN
C
C       Up to MXNIT iterations were done, the convergence rate is < 1,
C       a Jacobian routine is supplied, and the number of JACK calls
C       is less than MXNJ. Try again
C
         GO TO 300
      ENDIF
C
      IF (IERNEW .NE. 0) GO TO 380
C
C     Evaluate the residual for the quadrature only
      IF (ISENFO(10) .GT. 0) THEN
         DO I = 1, NSTATE
            IPS = (I-1)*NY + 1
            IRES = 3
            IF(ISENFO(6) .LT. 0) ISENFO(6) = -I
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,
     *           DELTA(IPS),IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT(IPS),CNST,SENPAR,G_RES, A_RES, ADI) 
            IF (IRES .LT. 0) THEN
               ISENFO(6) = -NSTATE
               GOTO 370
            END IF
            DO ISEN = 1, ISENFO(1) + 1
               IISEN = IPS + (ISEN-1)*NY
               DO IQ = 1, ISENFO(10)
                  IIQ = IISEN + NY - IQ
                  YPRIME(IIQ) = YPRIME(IIQ) - DELTA(IIQ)
               END DO
            END DO
         END DO
      END IF
      RETURN
C
C
C     Unsuccessful exits from nonlinear solver.
C     Set IERNLS accordingly.
C
 370  IERNLS = 2
      IF (IRES .LE. -2) IERNLS = -1
      RETURN
C
 380  IERNLS = MIN(IERNEW,2)
      RETURN
C
C----------------------- END OF SUBROUTINE DDASIK-----------------------
      END
      SUBROUTINE DNSIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,WT,RPAR,IPAR,
     *   SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON,
     *   RATEMX,MAXIT,STPTOL,ICNFLG,ICNSTR,IERNEW,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES,
     *   T_RES, LIADF,ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DNSIK
C***REFER TO  DDASPK
C***DATE WRITTEN   940701   (YYMMDD)
C***REVISION DATE  950714   (YYMMDD)
C***REVISION DATE  990503   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSIK solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
C     the initial conditions.
C
C     The method used is a Newton scheme combined with a linesearch
C     algorithm, using Krylov iterative linear system methods.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     ICOPT     -- Initial condition option chosen (1, 2 or 3).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  See explanation inside DDASPK.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     SAVR      -- Work vector for DNSIK of length NEQ.
C     DELTA     -- Residual vector on entry, and work vector of
C                  length NEQ for DNSIK.
C     R         -- Work vector for DNSIK of length NEQ.
C     YIC,YPIC  -- Work vectors for DNSIK, each of length NEQ.
C     PWK       -- Work vector for DNSIK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1,3) or 0 (ICOPT = 2).
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPLIN     -- Tolerance for linear system solver.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     STPTOL    -- Tolerance used in calculating the minimum lambda
C                  value allowed.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> failed to converge, but RATE .lt. 1., or the
C                          residual norm was reduced by a factor of .1.
C                   2  ==> failed to converge, RATE .gt. RATEMX.
C                   3  ==> other recoverable error.
C                  -1  ==> unrecoverable error inside Newton iteration.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference 
C                  sensitivity increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated 
C                  by ADIFOR for evaluations of sensitivity equations.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DFNRMK, DSLVK, DDWNRM, DLINSK, DCOPY
C
C***END PROLOGUE  DNSIK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ISENFO(*), ISENWK(*), SENWRK(*),SENPAR(*)
      DIMENSION Y(*),YPRIME(*),WT(*),ID(*),DELTA(*),R(*),SAVR(*)
      DIMENSION YIC(*),YPIC(*),PWK(*),WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION ICNSTR(*), ADI(*)
      EXTERNAL RES, PSOL, G_RES, A_RES, K_RES, T_RES
C
      PARAMETER (LNNI=19,LNPS=21,LLOCWP=29,LLCIWP=30,LNRE=12,LNSE=22)
      PARAMETER (LLSOFF=35, LNLI=20, LSTOL=14, LNY=37)
      PARAMETER (LLNIWP=28, LNLIS=38)
C
C
C     Initializations.  M is the Newton iteration counter.
C
      LSOFF = IWM(LLSOFF)
      M = 0
      RATE = 1.0D0
      LWP = IWM(LLOCWP)
      LIWP = IWM(LLCIWP)
      RLX = 0.4D0
      NY = IWM(LNY)
      NYMNQ = NY - ISENFO(10)
      NP = ISENFO(1)
C>>>>>>>>>>>>>>>>>>>>>>staggered method<<<<<<<<<<<<<<<<<<<<
C     
C.....for the state variable first
      NSTATE = 1
      LOCNNI = 0
      IF (ISENFO(6) .LT. -1) NSTATE = -ISENFO(6)
      DO I = 1, NSTATE
         IPS = (I-1)*NY+1
         IF (ISENFO(6).LT.0) ISENFO(6) = -I
C     
C     For the adjoint method
         IF (ISENFO(6).LT.-1 .AND. I .GT. 1) THEN
            IRES = 0
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,
     *           DELTA(IPS),IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT(IPS),CNST,SENPAR,G_RES, A_RES, ADI)
            IF (IRES .LT. 0) GOTO 490
         END IF
C     
         M = 0
         IF (ICOPT .EQ. 5) THEN
C.....evaluate the state variable first 
c     only for index-2 problem
C     
            IRES = 0
            CALL RESIDX2(X,Y(IPS),YPRIME(IPS),CJ,DELTA(IPS),
     *           IRES,RPAR,IPAR,SENPAR,
     *           WT(IPS),NY,ID(NY+1),RES,T_RES, 
     *           CNST,SENWRK,IWM(LIADF),ISENFO,A_RES,ADI)
            IF (IRES .LT. 0) GOTO 490
         END IF
C     
C     Save residual in SAVR.
C     
         CALL DCOPY (NYMNQ, DELTA(IPS), 1, SAVR, 1)
C     
C     Compute a new step vector DELTA.
C     
         CALL DSLVK (
     *        NYMNQ, Y(IPS),X,YPRIME(IPS),SAVR,DELTA(IPS),WT(IPS),
     *        WM,IWM,
     *        RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, 
     *        RHOK, RPAR, IPAR, SENPAR, ICOPT, ID, 
     *        ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
         IF (IRES .NE. 0 .OR. IERSL .NE. 0) GO TO 485
         DELNRM = DDWNRM(NY,DELTA(IPS),WT(IPS),RPAR,IPAR)
         FNRM = DELNRM
         IF (INDEX0 .EQ. 1) FNRM = TDIST*ABS(CJ)*FNRM
C     
C     Return now if residual norm is .le. EPCON.
C     
         IF (FNRM .LE. EPCON .OR. ISENFO(6) .LT. 0) GOTO 450
C     
C     Newton iteration loop.
C     
         FNRM0 = FNRM
 400     CONTINUE
         LOCNNI = LOCNNI + 1
C     
C     Call linesearch routine for global strategy and set RATE.
C     
         OLDFNM = FNRM
         CALL DLINSK (
     *        NEQ, Y(IPS), X, YPRIME(IPS), SAVR, CJ, DELTA(IPS), DELNRM, 
     *        WT(IPS),
     *        SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL,
     *        WM, IWM, RHOK, FNRM, ICOPT, ID,  
     *        R, EPLIN, YIC, YPIC, ICNFLG, ICNSTR, RLX, RPAR, 
     *        IPAR, ISENFO, SENWRK,ISENWK, CNST, SENPAR, G_RES, A_RES, 
     *        K_RES, T_RES, 0, LIADF, ADI, INDEX0, TDIST)
C     
         RATE = FNRM/OLDFNM
C     
C     Check for error condition from linesearch.
         IF (IRET .NE. 0) GO TO 485
C     
C     Test for convergence of the iteration, and return or loop.
C     
         IF (FNRM .LE. EPCON) GOTO 450
C     
C     The iteration has not yet converged.  Update M.
C     Test whether the maximum number of iterations have been tried.
C     
         M=M+1
         IF(M .GE. MAXIT) THEN
            IWM(LNNI) = IWM(LNNI) + LOCNNI/NSTATE
            IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
            GO TO 480
         END IF
C     
C     Copy the residual SAVR to DELTA and loop for another iteration.
C     
         CALL DCOPY (NYMNQ, R, 1, DELTA(IPS), 1)
         GO TO 400
C---------------------------------------------------------------------
 450     CONTINUE                  
      END DO
      IWM(LNNI) = IWM(LNNI) + LOCNNI/NSTATE
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
C>>>>>>>>>>>>>>>>>>>>>>>state variables are done<<<<<<<<<<<<<<<<<<<<<<
C     
      IF (ISENFO(1) .EQ. 0 .OR. ISENFO(6).LT.0) RETURN
      M = 0
      RATE = 1.0D0
      RLX = 0.4D0
C     
C     Reevaluate the residual again for the sensitivites and solution
      IF (ISENFO(2) .LT. 2) THEN
         IRES = 0
         CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES, ADI)
         CALL DCOPY(NYMNQ, DELTA, 1, R, 1) 
      END IF
      IF (ISENFO(2) .EQ. 5) THEN
C....................matrix times vector methods.....................
C     
C     evaluate the Jacobian and DF/DP
C     
         LJACI = 1
         LJACJ = LJACI + NY + 1
         LJAC = 1 + ISENFO(4)*NY
         CALL JRADFSP (
     1        NY, X, Y, YPRIME, DELTA, G_RES, CJ,  
     2        SENWRK(LJAC), ISENWK(LJACJ), ISENWK(LJACI), 
     3        IPAR, RPAR, SENPAR, IRES, SENWRK, ISENFO(4),
     4        IWM(LIADF))
      END IF
      IRES = 1
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      IF (IRES .LT. 0) GOTO 490     
      IF (ICOPT .EQ. 5) THEN
C.....evaluate the state variable first 
c     only for index-2 problem
C     
         IRES = 1
         CALL RESIDX2(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,
     *        WT,NY,ID(NY+1),RES,T_RES,
     *        CNST,SENWRK,IWM(LIADF),ISENFO,A_RES,ADI)
         IF (IRES .LT. 0) GOTO 490     
      END IF
C     
C     Save residual in SAVR.
C     
      CALL DCOPY (NEQ, DELTA, 1, SAVR, 1)
C     
C     Compute norm of residual.
C     
      NLI4STAT = IWM(LNLI)
      DO I = 1, NP
         II = I*NY + 1
         CALL DSLVK (
     *        NYMNQ, Y, X, YPRIME,SAVR,DELTA(II),WT(II),WM,IWM,
     *        RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, 
     *        RHOK, RPAR, IPAR, SENPAR, ICOPT, ID, 
     *        ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
         IF (IRES .NE. 0 .OR. IERSL .NE. 0) THEN
            IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-NLI4STAT+NP-1)/NP
            IWM(LNLI) = NLI4STAT
            GO TO 490
         END IF
      END DO
      IWM(LNLIS) = IWM(LNLIS) + (IWM(LNLI) - NLI4STAT + NP-1)/NP
      IWM(LNLI) = NLI4STAT
C     
C     Get norm of DELTA.  
C     
      DELNRM = 0.0D0
      DO I = 1, NP
         II = I*NY + 1     
         TNORM2 = DDWNRM(NYMNQ,DELTA(II),WT(II),RPAR,IPAR)
         IF (TNORM2 .GT. DELNRM) DELNRM = TNORM2
      END DO
      FNRM = DELNRM
      IF (INDEX0 .EQ. 1) FNRM = TDIST*ABS(CJ)*FNRM
C     
C     Return now if residual norm is .le. EPCON.
C     
      IF (FNRM .LE. EPCON) RETURN
C     
C     Newton iteration loop.
C     
 460  CONTINUE
C     
C     Compute a new step vector DELTA.
C     
      INDEX = -1
C     
C     Call linesearch routine for global strategy and set RATE.
C     
      OLDFNM = FNRM
C     
      CALL DLINSK (
     *     NEQ, Y, X, YPRIME, SAVR, CJ, DELTA, DELNRM, WT,
     *     SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL,
     *     WM, IWM, RHOK, FNRM, ICOPT, ID, 
     *     R, EPLIN, YIC, YPIC, ICNFLG, ICNSTR, RLX, RPAR, 
     *     IPAR, ISENFO, SENWRK,ISENWK,CNST, SENPAR, G_RES, A_RES, 
     *     K_RES, T_RES, INDEX, LIADF, ADI, INDEX0, TDIST)
C     
      RATE = FNRM/OLDFNM
C     
C     Check for error condition from linesearch.
      IF (IRET .NE. 0) GO TO 490
C     
C     Test for convergence of the iteration, and return or loop.
C     
      IF (FNRM .LE. EPCON) GOTO 500
C     
C     The iteration has not yet converged.  Update M.
C     Test whether the maximum number of iterations have been tried.
C     
      IF (ISENFO(7) .EQ. 2) THEN
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>staggered direct method<<<<<<<<<<<<<<<<<<
C     
         RETURN
      END IF         
      M=M+1
      IF(M .GE. MAXIT) GO TO 480
C     
C     Copy the residual SAVR to DELTA and loop for another iteration.
C     
      CALL DCOPY (NEQ-NY,  R(NY+1), 1, DELTA(NY+1), 1)
      GO TO 460
C
C     The maximum number of iterations was done.  Set IERNEW and return.
C
 480  IF (RATE .LE. RATEMX .OR. FNRM .LE. 0.1D0*FNRM0) THEN
         IERNEW = 1
      ELSE
         IERNEW = 2
      ENDIF
      RETURN
C
 485  CONTINUE
      IWM(LNNI) = IWM(LNNI) + LOCNNI/NSTATE
C     
 490  IF (IRES .LE. -2 .OR. IERSL .LT. 0) THEN
         IERNEW = -1
      ELSE
         IERNEW = 3
         IF (IRES .GE. 0 .AND. IERSL .EQ. 1 .AND. M .GE. 2 
     1        .AND. RATE .LT. 1.0D0) IERNEW = 1
      ENDIF
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
 500  CONTINUE
      RETURN
C
C
C----------------------- END OF SUBROUTINE DNSIK------------------------
      END
      SUBROUTINE DLINSK (NEQ, Y, X, YPRIME, SAVR, CJ, P, PNRM, WT,
     *   SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL, WM, IWM,
     *   RHOK, FNRM, ICOPT, ID, R, EPLIN, YNEW, YPNEW, 
     *   ICNFLG, ICNSTR, RLX, RPAR, IPAR,
     *   ISENFO, SENWRK, ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES, T_RES, 
     *   INDEX, LIADF,ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DLINSK
C***REFER TO  DNSIK
C***DATE WRITTEN   940830   (YYMMDD)
C***REVISION DATE  951006   (Arguments SQRTN, RSQRTN added.)
C***REVISION DATE  960129   Moved line RL = ONE to top block.
C***REVISION DATE  990505   Rewrited to be unified with DLINSD
C***REVISION DATE  000628   RHOK*RHOK term removed in alpha test.
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DLINSK uses a linesearch algorithm to calculate a new (Y,YPRIME)
C     pair (YNEW,YPNEW) such that 
C
C     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME)
C
C     where 0 < RL <= 1, and f(y,y') is defined as
C
C      f(y,y') = (1/2)*norm( (P-inverse)*G(t,y,y') )**2 ,
C
C     where norm() is the weighted RMS vector norm, G is the DAE
C     system residual function, and P is the preconditioner used
C     in the Krylov iteration.
C
C     In addition to the parameters defined elsewhere, we have
C
C     SAVR    -- Work array of length NEQ.
C     P       -- Approximate Newton step used in backtracking.
C     PNRM    -- Weighted RMS norm of P.
C     LSOFF   -- Flag showing whether the linesearch algorithm is
C                to be invoked.  0 means do the linesearch, 
C                1 means turn off linesearch.
C     STPTOL  -- Tolerance used in calculating the minimum lambda
C                value allowed.
C     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations
C                in the proposed new approximate solution will be
C                checked for, and the maximum step length will be
C                adjusted accordingly.
C     ICNSTR  -- Integer array of length NEQ containing flags for
C                checking constraints.
C     RHOK    -- Weighted norm of preconditioned Krylov residual.
C     RLX     -- Real scalar restricting update size in DCNSTR.
C     YNEW    -- Array of length NEQ used to hold the new Y in
C                performing the linesearch.
C     YPNEW   -- Array of length NEQ used to hold the new YPRIME in
C                performing the linesearch.
C     PWK     -- Work vector of length NEQ for use in PSOL.
C     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW).
C     YPRIME  -- Array of length NEQ containing the new YPRIME 
C                (i.e.,=YPNEW).
C     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the
C                current (Y,YPRIME) on input and output.
C     R       -- Work space length NEQ, containing the residual
C                vector G(t,y,y') on return.
C     IRET    -- Return flag.
C                IRET=0 means that a satisfactory (Y,YPRIME) was found.
C                IRET=1 means that the routine failed to find a new
C                       (Y,YPRIME) that was sufficiently distinct from
C                       the current (Y,YPRIME) pair.
C                IRET=2 means a failure in RES or PSOL.
C     ISENFO --  Sensitivity analysis information vector.
C     SENWRK --  Sensitivity analysis work vector.
C     ISENWK --  Sensitivity analysis integer work vector.
C     CNST --    Optional user constant for finite difference 
C                sensitivity increment.
C     SENPAR --  Sensitivity parameter array.
C     G_RES  --  External user-supplied subroutine which is generated by
C                ADIFOR for evaluations of sensitivity equations.
C     INDEX   -- Indicator for the corrector method.
C                INDEX=0, for the state variables only;
C                INDEX=1, for both the state variables and 
C                             sensitivity variables
C                INDEX=-1, for the sensitivity variables only
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DFNRMK, DYYPNW, DCOPY
C
C***END PROLOGUE  DLINSK
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL  RES, PSOL, G_RES, A_RES, K_RES, T_RES
      DIMENSION Y(*), YPRIME(*), P(*), WT(*), SAVR(*), R(*), ID(*)
      DIMENSION WM(*), IWM(*), YNEW(*), YPNEW(*), ICNSTR(*)
      DIMENSION RPAR(*), IPAR(*), ADI(*)
      DIMENSION ISENFO(*), ISENWK(*), SENWRK(*),SENPAR(*)
C      CHARACTER MSG*80
C
      PARAMETER (LNRE=12, LNPS=21, LKPRIN=31,LNSE=22,LNY=37)
C
      SAVE ALPHA, ONE, TWO
      DATA ALPHA/1.0D-4/, ONE/1.0D0/, TWO/2.0D0/
C
C      KPRIN=IWM(LKPRIN)
      F1NRM = (FNRM*FNRM)/TWO
      RATIO = ONE
C
C      IF (KPRIN .GE. 2) THEN
C         MSG = '------ IN ROUTINE DLINSK-- PNRM = (R1) )'
C         CALL XERRWD(MSG, 40, 921, 0, 0, 0, 0, 1, PNRM, 0.0D0)
C      ENDIF
      TAU = PNRM
      IVIO = 0
      RL = ONE
      NY = IWM(LNY)
C-----------------------------------------------------------------------
C Check for violations of the constraints, if any are imposed.
C If any violations are found, the step vector P is rescaled, and the 
C constraint check is repeated, until no violations are found.
C-----------------------------------------------------------------------
      IF (INDEX .LT. 0) GOTO 20
      IF (ICNFLG .NE. 0) THEN
 10      CONTINUE
         CALL DYYPNW (NY,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
         CALL DCNSTR (NY, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
         IF (IRET .EQ. 1) THEN
            IVIO = 1
            RATIO1 = TAU/PNRM
            RATIO = RATIO*RATIO1
            DO I = 1,NY
               P(I) = P(I)*RATIO1
            END DO
            PNRM = TAU
C            IF (KPRIN .GE. 2) THEN
C            MSG = '------ CONSTRAINT VIOL., PNRM = (R1), INDEX = (I1)'
C            CALL XERRWD(MSG, 50, 922, 0, 1, IVAR, 0, 1, PNRM, 0.0D0)
C            ENDIF
            IF (PNRM .LE. STPTOL) THEN
               IRET = 1
               RETURN
            ENDIF
            GO TO 10
         ENDIF
      ENDIF
C
 20   CONTINUE
      SLPI = (-TWO*F1NRM)*RATIO
      RLMIN = STPTOL/PNRM
C      IF (LSOFF .EQ. 0 .AND. KPRIN .GE. 2) THEN
C        MSG = '------ MIN. LAMBDA = (R1)'
C        CALL XERRWD(MSG, 25, 923, 0, 0, 0, 0, 1, RLMIN, 0.0D0)
C      ENDIF
C-----------------------------------------------------------------------
C Begin iteration to find RL value satisfying alpha-condition.
C Update YNEW and YPNEW, then compute norm of new scaled residual and
C perform alpha condition test.
C-----------------------------------------------------------------------
      IF (INDEX .EQ.-1) THEN
         CALL DCOPY (NY, Y, 1, YNEW, 1)
         CALL DCOPY (NY, YPRIME, 1, YPNEW, 1)
      END IF
 100  CONTINUE
      IF (INDEX .GE. 0) 
     *     CALL DYYPNW (NY,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
      IF (INDEX .NE. 0) THEN
         NP = ISENFO(1)
         DO I = 1, NP
            II = I*NY + 1
            CALL DYYPNW (NY,Y(II),YPRIME(II),CJ,RL,
     *           P(II),ICOPT,ID,YNEW(II),YPNEW(II))
         END DO 
      END IF
      CALL DFNRMK (
     *  NY, YNEW, X, YPNEW, SAVR, R, CJ, WT, SQRTN, RSQRTN,
     *  RES, IRES, PSOL, IER, FNRMP, EPLIN, WM, IWM, RHOK, 
     *  RPAR, IPAR,
     *  ISENFO, SENWRK, ISENWK, CNST, SENPAR, G_RES, A_RES, K_RES,INDEX,
     *  ICOPT, ID, LIADF, T_RES, ADI, INDEX0, TDIST)
      IF (IRES .LT. 0 .OR. IER .LT. 0) THEN
        IRET = 2
        RETURN
      ENDIF
      IF (IER .GT. 0) THEN
         IRET = 1
         RETURN
      END IF
      IF (LSOFF .EQ. 1) GO TO 150
C
      F1NRMP = FNRMP*FNRMP/TWO
C      IF (KPRIN .GE. 2) THEN
C        MSG = '------ LAMBDA = (R1)'
C        CALL XERRWD(MSG, 20, 924, 0, 0, 0, 0, 1, RL, 0.0D0)
C        MSG = '------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)'
C        CALL XERRWD(MSG, 43, 925, 0, 0, 0, 0, 2, F1NRM, F1NRMP)
C      ENDIF
      IF (F1NRMP .GT. F1NRM + ALPHA*SLPI*RL) GO TO 200
C-----------------------------------------------------------------------
C Alpha-condition is satisfied, or linesearch is turned off.
C Copy YNEW,YPNEW to Y,YPRIME and return.
C-----------------------------------------------------------------------
 150  IRET = 0
      IF (INDEX .GE. 0) THEN
         CALL DCOPY(NY, YNEW, 1, Y, 1)
         CALL DCOPY(NY, YPNEW, 1, YPRIME, 1)
      END IF
      IF (INDEX .NE. 0) THEN
         CALL DCOPY(NEQ-NY, YNEW(NY+1), 1, Y(NY+1), 1)
         CALL DCOPY(NEQ-NY, YPNEW(NY+1), 1, YPRIME(NY+1), 1)
      END IF
      FNRM = FNRMP
C      IF (KPRIN .GE. 1) THEN
C        MSG = '------ LEAVING ROUTINE DLINSK, FNRM = (R1)'
C        CALL XERRWD(MSG, 42, 926, 0, 0, 0, 0, 1, FNRM, 0.0D0)
C        ENDIF
      RETURN
C-----------------------------------------------------------------------
C Alpha-condition not satisfied.  Perform backtrack to compute new RL
C value.  If RL is less than RLMIN, i.e. no satisfactory YNEW,YPNEW can
C be found sufficiently distinct from Y,YPRIME, then return IRET = 1.
C-----------------------------------------------------------------------
 200  CONTINUE
      IF (RL .LT. RLMIN) THEN
        IRET = 1
        RETURN
      ENDIF
C
      RL = RL/TWO
      GO TO 100
C
C----------------------- END OF SUBROUTINE DLINSK ----------------------
      END
      SUBROUTINE DFNRMK (NY, Y, X, YPRIME, SAVR, R, CJ, WT,
     *                   SQRTN, RSQRTN, RES, IRES, PSOL, IER,
     *                   FNORM, EPLIN, WM, IWM, RHOK, RPAR, IPAR,
     *                   ISENFO, SENWRK, ISENWK,CNST, 
     *                   SENPAR, G_RES, A_RES, K_RES, INDEX, ICOPT, ID, 
     *                   LIADF, T_RES, ADI, INDEX0, TDIST)
C
C***BEGIN PROLOGUE  DFNRMK
C***REFER TO  DLINSK
C***DATE WRITTEN   940830   (YYMMDD)
C***REVISION DATE  951006   (SQRTN, RSQRTN, and scaling of WT added.)
C***REVISION DATE  990505   Rewrited to be unified with DFNRMD
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DFNRMK calculates the scaled preconditioned norm of the nonlinear
C     function used in the nonlinear iteration for obtaining consistent
C     initial conditions.  Specifically, DFNRMK calculates the weighted
C     root-mean-square norm of the vector (P-inverse)*G(T,Y,YPRIME),
C     where P is the preconditioner matrix.
C
C     In addition to the parameters described in the calling program
C     DLINSK, the parameters represent
C
C     R      -- Array of length NEQ that contains
C               (J-inverse)*G(T,Y,YPRIME) on return.
C     FNORM  -- Scalar containing the weighted norm of R on return.
C     INDEX  -- Indicator for the corrector method.
C               INDEX=0, for the state variables only;
C               INDEX=1, for both the state variables and 
C                            sensitivity variables
C               INDEX=-1, for the sensitivity variables only
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DDSEN, DCOPY, DSCAL, PSOL, DDWNRM
C
C***END PROLOGUE  DFNRMK
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RES, PSOL, G_RES, A_RES, K_RES, T_RES
      DIMENSION Y(*), YPRIME(*), WT(*), SAVR(*), R(*)
      DIMENSION WM(*), IWM(*), RPAR(*), IPAR(*), ADI(*)
      DIMENSION ISENFO(*), ISENWK(*),SENWRK(*),SENPAR(*),ID(*)
      PARAMETER (LNRE=12, LNLI=20, LNSE=22, LNLIS=38)
      PARAMETER (LLCIWP=30, LLNIWP=28)
C-----------------------------------------------------------------------
C     Call DDSEN (RES) routine 
C-----------------------------------------------------------------------
      NP = ISENFO(1)
      IF (INDEX .EQ. 0) THEN
         IRES = 0
         CALL DDSEN(X,Y,YPRIME,CJ,R,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST, 
     *        SENPAR,G_RES, A_RES, ADI)
      ELSE 
         IRES = 1
         CALL DDSEN(X,Y,YPRIME,CJ,R,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST, 
     *        SENPAR,G_RES, A_RES, ADI)
      END IF
      IF (ICOPT .EQ. 5) THEN
         CALL RESIDX2(X,Y,YPRIME,CJ,R,IRES,RPAR,IPAR,SENPAR,
     *        WT,NY,ID(NY+1),RES,T_RES, CNST,SENWRK,
     *        IWM(LIADF),ISENFO,A_RES,ADI)
      END IF
      IF (IRES .LT. 0) RETURN
C-----------------------------------------------------------------------
C     Apply inverse of Jacobian to vector R.
C-----------------------------------------------------------------------
      NYMNQ = NY - ISENFO(10)      
      IF (INDEX .GE. 0) CALL DCOPY(NY, R, 1, SAVR, 1)
      IF (INDEX .GE. 0) THEN
         CALL DSLVK (
     *        NYMNQ, Y, X, YPRIME,SAVR,R,WT,WM,IWM,
     *        RES, IRES, PSOL, IER, CJ, EPLIN, SQRTN, RSQRTN, 
     *        RHOK, RPAR, IPAR, SENPAR, ICOPT, ID, 
     *        ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
         IF (IRES .NE. 0 .OR. IER .NE. 0) RETURN
      END IF
      IF (INDEX .NE. 0) THEN
         NLI4STAT = IWM(LNLI)
         DO I = 1, NP
            II = I*NY + 1
            CALL DSLVK (
     *           NYMNQ, Y, X, YPRIME,SAVR,R(II),WT(II),WM,IWM,
     *           RES, IRES, PSOL, IER, CJ, EPLIN, SQRTN, RSQRTN, 
     *           RHOK, RPAR, IPAR, SENPAR, ICOPT, ID, 
     *           ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
            IF (IRES .NE. 0 .OR. IER .NE. 0) THEN
               IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-NLI4STAT+NP-1)/NP
               IWM(LNLI) = NLI4STAT
               RETURN
            END IF
         END DO
         IF (ISENFO(1) .GT. 0) THEN
            IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-NLI4STAT+NP-1)/NP
            IWM(LNLI) = NLI4STAT
         END IF
      END IF
C-----------------------------------------------------------------------
C     Calculate norm of R.
C-----------------------------------------------------------------------
      FNORM = 0.0D0
      IF (INDEX .GE. 0) FNORM = DDWNRM(NYMNQ,R,WT,RPAR,IPAR)
      IF (INDEX .NE. 0) THEN
         DO I = 1, NP
            II = I*NY + 1     
            TNORM2 = DDWNRM(NYMNQ,R(II),WT(II),RPAR,IPAR)
            IF (TNORM2 .GT. FNORM) FNORM = TNORM2
         END DO
      END IF
      IF (INDEX0. EQ. 1) FNORM = TDIST*ABS(CJ)*FNORM
C
      RETURN
C----------------------- END OF SUBROUTINE DFNRMK ----------------------
      END
      SUBROUTINE DNEDK(X,Y,YPRIME,NEQ,RES,JACK,PSOL,
     *   H,WT,JSTART,IDID,RPAR,IPAR,PHI,NPHI,GAMMA,SAVR,DELTA,E,
     *   WM,IWM,CJ,CJOLD,CJLAST,S,UROUND,EPLI,SQRTN,RSQRTN,
     *   EPCON,JCALC,JFLG,KP1,NONNEG,NTYPE,IERNLS,
     *   CK,VT,ENORM,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES, 
     *   LIADF,ADI)
C
C***BEGIN PROLOGUE  DNEDK
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD)
C***REVISION DATE  990505   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNEDK solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a matrix-free Newton scheme.
C
C     The parameters represent
C     X         -- Independent variable.
C     Y         -- Solution vector at x.
C     YPRIME    -- Derivative of solution vector
C                  after successful step.
C     NEQ       -- Number of equations to be integrated.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JACK     --  External user-supplied routine to update
C                  the preconditioner.  (This is optional).
C                  See JAC description for the case
C                  INFO(12) = 1 in the DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  (This is optional).  See explanation inside DDASPK.
C     H         -- Appropriate step size for this step.
C     WT        -- Vector of weights for error criterion.
C     JSTART    -- Indicates first call to this routine.
C                  If JSTART = 0, then this is the first call,
C                  otherwise it is not.
C     IDID      -- Completion flag, output by DNEDK.

C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     PHI       -- Array of divided differences used by
C                  DNEDK.  The length is NEQ*(K+1), where
C                  K is the maximum order.
C     GAMMA     -- Array used to predict Y and YPRIME.  The length
C                  is K+1, where K is the maximum order.
C     SAVR      -- Work vector for DNEDK of length NEQ.
C     DELTA     -- Work vector for DNEDK of length NEQ.
C     E         -- Error accumulation vector for DNEDK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information for linear system
C                  solvers, and various other information.
C     CJ        -- Parameter always proportional to 1/H.
C     CJOLD     -- Saves the value of CJ as of the last call to JACK.
C                  Accounts for changes in CJ needed to
C                  decide whether to call JACK.
C     CJLAST    -- Previous value of CJ.
C     S         -- A scalar determined by the approximate rate
C                  of convergence of the Newton iteration and used
C                  in the convergence test for the Newton iteration.
C
C                  If RATE is defined to be an estimate of the
C                  rate of convergence of the Newton iteration,
C                  then S = RATE/(1.D0-RATE).
C
C                  The closer RATE is to 0., the faster the Newton
C                  iteration is converging; the closer RATE is to 1.,
C                  the slower the Newton iteration is converging.
C
C                  On the first Newton iteration with an up-dated
C                  preconditioner S = 100.D0, Thus the initial
C                  RATE of convergence is approximately 1.
C
C                  S is preserved from call to call so that the rate
C                  estimate from a previous step can be applied to
C                  the current step.
C     UROUND    -- Unit roundoff.
C     EPLI      -- convergence test constant.
C                  See DDASPK prologue for more details.
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     JCALC     -- Flag used to determine when to update
C                  the Jacobian matrix.  In general:
C
C                  JCALC = -1 ==> Call the DITMD routine to update
C                                 the Jacobian matrix.
C                  JCALC =  0 ==> Jacobian matrix is up-to-date.
C                  JCALC =  1 ==> Jacobian matrix is out-dated,
C                                 but DITMD will not be called unless
C                                 JCALC is set to -1.
C     JFLG      -- Flag showing whether a Jacobian routine is supplied.
C     KP1       -- The current order + 1;  updated across calls.
C     NONNEG    -- Flag to determine nonnegativity constraints.
C     NTYPE     -- Identification code for the DNEDK routine.
C                   1 ==> modified Newton; iterative linear solver.
C                   2 ==> modified Newton; user-supplied linear solver.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0 ==> nonlinear solver converged.
C                   1 ==> recoverable error inside non-linear solver.
C                  -1 ==> unrecoverable error inside non-linear solver.
C                  -2 ==> error test failure for state variables in 
C                         staggered corrector method.
C     ISENFO --    Sensitivity analysis information vector.
C     SENWRK --    Sensitivity analysis work vector.
C     ISENWK --    Sensitivity analysis integer work vector.
C     CNST --      Optional user constant for finite difference sensitivity
C                  increment.
C     SENPAR --    Sensitivity parameter array.
C     G_RES  --    External user-supplied subroutine which is generated by
C                  ADIFOR for evaluations of sensitivity equations.
C
C     The following group of variables are passed as arguments to
C     the Newton iteration solver.  They are explained in greater detail
C     in DNSK:
C        TOLNEW, MULDEL, MAXIT, IERNEW
C
C     IERTYP -- Flag which tells whether this subroutine is correct.
C               0 ==> correct subroutine.
C               1 ==> incorrect subroutine.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   RES, JACK, DDWNRM, DNSK, DDSEN
C
C***END PROLOGUE  DNEDK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),VT(*)
      DIMENSION PHI(NPHI,*),SAVR(*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*)
      DIMENSION GAMMA(*),RPAR(*),IPAR(*), ADI(*)
      DIMENSION ISENFO(*),ISENWK(*),SENWRK(*),SENPAR(*)
      EXTERNAL  RES, JACK, PSOL, G_RES, A_RES, K_RES
C
      PARAMETER (LNRE=12,LNJE=13,LLOCWP=29,LLCIWP=30,LNSE=22,LNY=37)
      PARAMETER (LMTYPE=4)
CC
      SAVE MULDEL, MAXIT, XRATE
      DATA MULDEL/0/, MAXIT/4/, XRATE/0.25D0/
C
C     Verify that this is the correct subroutine.
C
      IERTYP = 0
      NSTATE = 1
      IF (NTYPE .NE. 1) THEN
         IERTYP = 1
         GO TO 380
      ENDIF
C
C     If this is the first step, perform initializations.
C
      IF (JSTART .EQ. 0) THEN
         CJOLD = CJ
         JCALC = -1
      ENDIF
C
C     Perform all other initializations.
C
      IERNLS = 0
      LWP = IWM(LLOCWP)
      LIWP = IWM(LLCIWP)
      NY = IWM(LNY)
      NP = ISENFO(1)
C
C     Decide whether to update the preconditioner.
C
      IF (JFLG .NE. 0) THEN
         TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
         TEMP2 = 1.0D0/TEMP1
         IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
         IF (CJ .NE. CJLAST) THEN
            S = 100.D0
         END IF
      ELSE
         JCALC = 0
      ENDIF
C
C     Looping point for updating preconditioner with current stepsize.
C
300   CONTINUE
C
C     Initialize all error flags to zero.
C
      IERPJ = 0
      IERSL = 0
      IERNEW = 0
C
C     Predict the solution and derivative and compute the tolerance
C     for the Newton iteration.
C
      DO 310 I=1,NPHI
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NPHI
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
C
      EPLIN = EPLI*EPCON
      TOLNEW = EPLIN
C
C     Call DDSEN(RES) to initialize DELTA.
C
      IRES = 0
      NYMNQ = NY - ISENFO(10)
      IF (ISENFO(6) .LT. 0) THEN
         NSTATE = -ISENFO(6)
         ISENFO(6) = -1
         ISENFO(1) = NSTATE
      END IF
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      IF (IRES .LT. 0) GO TO 380
C
C
C     If indicated, update the preconditioner.
C     Set JCALC to 0 as an indicator that this has been done.
C
      IF(JCALC .EQ. -1)THEN
         IWM(LNJE) = IWM(LNJE) + 1
         JCALC=0
         IRES = 0
         IF (ISENFO(6) .LT. 0) THEN
            CALL JACKAD (JACK,
     *           RES, IRES, NYMNQ, X, Y, YPRIME, WT, DELTA,E,H,CJ,
     *           WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR,SENPAR, ADI,
     *           ISENFO,SENWRK,ISENWK,G_RES, A_RES,0,ID)
         ELSE
            CALL JACK (RES,IRES,NYMNQ,X,Y, YPRIME, WT, DELTA, E, H, CJ,
     *           WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR, SENPAR,0)
         END IF
         CJOLD=CJ
         S = 100.D0
         IF (IRES .LT. 0)  GO TO 380
         IF (IERPJ .NE. 0) GO TO 380
      ENDIF
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
C
C     Call the nonlinear Newton solver.
C
      IF (NPHI .EQ. 2*NEQ) THEN
         DO I = 1, NY
            DELTA(I) = DELTA(I) - CJ*DELTA(I+NEQ)
         END DO
      END IF
      IF (IWM(LMTYPE) .EQ. 3) THEN
         TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
         MULDEL = 1
      END IF
      CALL DNSK(X,Y,YPRIME,NEQ,NPHI,RES,PSOL,WT,RPAR,IPAR,SAVR,
     *   DELTA,E,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON,
     *   S,TEMP1,TOLNEW,MULDEL,MAXIT,IRES,IERSL,IERNEW,
     *   CK,VT,ENORM,NONNEG,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES,
     *   LIADF,ADI)
      IF (ISENFO(6) .LT. 0) ISENFO(1) = 0
C
      IF (IERNEW .GT. 0 .AND. JCALC .NE. 0) THEN
C
C     The Newton iteration had a recoverable failure with an old
C     preconditioner.  Retry the step with a new preconditioner.
C
         JCALC = -1
         GO TO 300
      ENDIF
C
      IF (IERNEW .EQ. 0) THEN
C
C     Evaluate the residual for the sensitivity of quadrature only
         IF (ISENFO(10) .GT. 0 .AND. ISENFO(1) .GT. 0) THEN
            IRES = 3
            CALL DDSEN(X,Y,YPRIME,CJ,
     *           DELTA,IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT,CNST,SENPAR,G_RES, A_RES, ADI) 
            IF (IRES .LT. 0) GOTO 380
            CJINV = 1.0D0/CJ
            DO ISEN = 1, ISENFO(1)
               IISEN = 1 + ISEN*NY
               DO IQ = NYMNQ, NY-1
                  IIQ = IISEN + IQ
                  YPRIME(IIQ) = YPRIME(IIQ) - DELTA(IIQ)
                  Y(IIQ) = Y(IIQ) - CJINV*DELTA(IIQ)
                  E(IIQ) = -CJINV*DELTA(IIQ)
               END DO
            END DO
         END IF
         GO TO 390
      END IF
C
C
C     Exits from nonlinear solver.
C     No convergence with current preconditioner.
C     Compute IERNLS and IDID accordingly.
C
380   CONTINUE
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
      IF (IRES .LE. -2 .OR. IERSL .LT. 0 .OR. IERTYP .NE. 0) THEN
         IERNLS = -1
         IF (IRES .LE. -2) IDID = -11
         IF (IERSL .LT. 0) IDID = -13
         IF (IERTYP .NE. 0) IDID = -15
      ELSE IF (IERNEW .EQ. -2) THEN
         IERNLS = -2
      ELSE 
         IERNLS =  1
         IF (IRES .EQ. -1) IDID = -10
         IF (IERPJ .NE. 0) IDID = -5
         IF (IERSL .GT. 0) IDID = -14
      ENDIF
C
C
390   JCALC = 1
      RETURN
C
C------END OF SUBROUTINE DNEDK------------------------------------------
      END
      SUBROUTINE DNSK (X,Y,YPRIME,NEQ,NPHI,RES,PSOL,WT,RPAR,IPAR,
     *   SAVR,DELTA,E,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON,
     *   S,CONFAC,TOLNEW,MULDEL,MAXIT,IRES,IERSL,IERNEW,
     *   CK, VT, ENORM,NONNEG,
     *   ISENFO,SENWRK,ISENWK,CNST,SENPAR,G_RES, A_RES, K_RES, 
     *   LIADF,ADI)
C
C***BEGIN PROLOGUE  DNSK
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  950126   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSK solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  See explanation inside DDASPK.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DDASPK.
C     SAVR      -- Work vector for DNSK of length NEQ.
C     DELTA     -- Work vector for DNSK of length NEQ.
C     E         -- Error accumulation vector for DNSK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Parameter always proportional to 1/H (step size).
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPLIN     -- Tolerance for linear system solver.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     S         -- Used for error convergence tests.
C                  In the Newton iteration: S = RATE/(1.D0-RATE),
C                  where RATE is the estimated rate of convergence
C                  of the Newton iteration.
C
C                  The closer RATE is to 0., the faster the Newton
C                  iteration is converging; the closer RATE is to 1.,
C                  the slower the Newton iteration is converging.
C
C                  The calling routine sends the initial value
C                  of S to the Newton iteration.
C     CONFAC    -- A residual scale factor to improve convergence.
C     TOLNEW    -- Tolerance on the norm of Newton correction in
C                  alternative Newton convergence test.
C     MULDEL    -- A flag indicating whether or not to multiply
C                  DELTA by CONFAC.
C                  0  ==> do not scale DELTA by CONFAC.
C                  1  ==> scale DELTA by CONFAC.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     IRES      -- Error flag returned from RES.  See RES description
C                  in DDASPK prologue.  If IRES = -1, then IERNEW
C                  will be set to 1.
C                  If IRES < -1, then IERNEW will be set to -1.
C     IERSL     -- Error flag for linear system solver.
C                  See IERSL description in subroutine DSLVK.
C                  If IERSL = 1, then IERNEW will be set to 1.
C                  If IERSL < 0, then IERNEW will be set to -1.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> recoverable error inside Newton iteration.
C                  -1  ==> unrecoverable error inside Newton iteration.
C                  -2  ==> error test failure for state variables in
C                          staggered corrector method.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, DDSEN, DSLVK, DDWNRM
C
C***END PROLOGUE  DNSK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),DELTA(*),E(*),SAVR(*),VT(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION ISENFO(*),ISENWK(*),SENWRK(*),SENPAR(*), ADI(*)
      EXTERNAL  RES, PSOL, G_RES, A_RES, K_RES
C
      PARAMETER (LNRE=12, LNNI=19, LNLI=20,LNSE=22,LNY=37)
      PARAMETER (LLCIWP=30, LLNIWP=28, LNLIS=38)
      
C
C     Initialize Newton counter M and accumulation vector E.
C
      M = 0
      NY = IWM(LNY)
      NYMNQ = NY - ISENFO(10)
      NP = ISENFO(1)
      DO 100 I=1,NPHI
100     E(I) = 0.0D0
C
C     Corrector loop.
C
C>>>>>>>>>>>>>>>>>>>>>Staggered method<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C
C     First compute the solution Y, not including the sensitivities
C
      NSTATE = 1
      IF (ISENFO(6) .LT. -1) NSTATE = - ISENFO(6)
      LOCNNI = 0
      IPS = 1
      ENORM = 0.0D0
      DO ISTATE = 1, NSTATE
C     
C     For the adjoint method
         IF (ISENFO(6).LT.0) ISENFO(6) = -ISTATE
         IF (NSTATE .GT. 1 .AND. ISTATE .GT. 1) THEN
            IRES = 0
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,DELTA(IPS),
     *           IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT(IPS),
     *           CNST,SENPAR,G_RES, A_RES,ADI)
            IF (NPHI .EQ. 2*NEQ) THEN
               DO I = 0, NY-1
                  IIPS = I + IPS
                  DELTA(IIPS) = DELTA(IIPS) - CJ*DELTA(IIPS+NEQ)
               END DO
            END IF
            IF (IRES .LT. 0) GOTO 575
         END IF
 400     CONTINUE
C
C     If necessary, multiply residual by convergence factor.
C
         IF (MULDEL .EQ. 1) THEN
            DO I = 0,NYMNQ-1
               DELTA(I+IPS) = DELTA(I+IPS) * CONFAC
            end do
         ENDIF
         DO I = 1,NYMNQ
            SAVR(I) = DELTA(I+IPS-1)
         END DO
         LOCNNI = LOCNNI + 1
C
C     Compute a new iterate.  Store the correction in DELTA.
C
         IRES = 0
         CALL DSLVK (
     *        NYMNQ, Y(IPS), X, YPRIME(IPS), SAVR, DELTA(IPS), 
     *        WT(IPS), WM, IWM,
     *        RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, 
     *        RSQRTN, RHOK, RPAR, IPAR, SENPAR, 0, ID, 
     *        ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
         IF (IRES .NE. 0 .OR. IERSL .NE. 0) GO TO 575
C
C     Update Y, E, and YPRIME.
C
         DO I=0,NYMNQ-1
            II = I+IPS
            Y(II)=Y(II)-DELTA(II)
            E(II)=E(II)-DELTA(II)
            YPRIME(II)=YPRIME(II)-CJ*DELTA(II)
         END DO
         IF (NPHI .EQ. 2*NEQ) THEN
            DO I = 0, NYMNQ-1
	       II = IPS + NEQ + I
               E(II) = E(II) - DELTA(II)
	       Y(II) = Y(II) - DELTA(II)
	       YPRIME(II) = YPRIME(II) - CJ*DELTA(II)
            END DO
C
C     Calculate the vector-matrix product -- v(F_y')
            CALL FXPV(X,DELTA(IPS),CJ,DELTA(NEQ+IPS),IRES,RPAR,
     *		 IPAR,RES,NY,
     *           ISENFO,ISENWK,SENPAR,ADI, 0)
            DO I = 0, NYMNQ-1
               II = NEQ+I+IPS
               E(II) = E(II) - DELTA(II)
	       Y(II) = Y(II) - DELTA(II)
	       YPRIME(II) =YPRIME(II) - CJ*DELTA(II)
            END DO
         END IF
C
C     For Linear system, it should converge in one Newton iteration
C
         IF (ISENFO(6) .LT. 0) GOTO 470
C
C     Test for convergence of the iteration.
C
         DELNRM=DDWNRM(NYMNQ,DELTA(IPS),WT(IPS),RPAR,IPAR)
         IF (M .EQ. 0) THEN
            OLDNRM = DELNRM
            IF (DELNRM .LE. TOLNEW) GO TO 470
         ELSE
            RATE = (DELNRM/OLDNRM)**(1.0D0/M)
            IF (RATE .GT. 0.9D0) GO TO 575
            S = RATE/(1.0D0 - RATE)
         ENDIF
         IF (S*DELNRM .LE. EPCON) GO TO 470
C
C     The corrector has not yet converged.
C     Update M and test whether the
C     maximum number of iterations have
C     been tried.
C
         M=M+1
         IF(M.GE.MAXIT) GO TO 575
C
C     Evaluate the residual,
C     and go back to do another iteration.
C
         IRES = 0
         CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,DELTA(IPS),
     *        IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT(IPS),CNST,
     *        SENPAR,G_RES, A_RES, ADI)
         IF (IRES .LT. 0) GO TO 575
         GO TO 400
C
C     The iteration for solution has converged.
C
 470     M = 0
C
C     The Newton iteration for state variable has converged.  
C     If nonnegativity of
C     solution is required, set the solution nonnegative, if the
C     perturbation to do it is small enough.  If the change is too
C     large, then consider the corrector iteration to have failed.
C
         IF(NONNEG .EQ. 0) GO TO 475
         DO I = 0, NYMNQ-1
            DELTA(IPS+I) = MIN(Y(IPS+I),0.0D0)
         END DO
         DELNRM = DDWNRM(NYMNQ,DELTA(IPS),WT(IPS),RPAR,IPAR)
         IF(DELNRM .GT. EPCON) THEN
            IERNEW = -1
            RETURN
         END IF
         DO I = 0,NYMNQ-1
            E(IPS+I) = E(IPS+I) - DELTA(IPS+I)
         END DO
 475     CONTINUE
C
C     Evaluate the residual ONLY for the quadrature.
         IF (ISENFO(10) .GT. 0) THEN
            IRES = 3
            ISENFO(1) = 0
            CALL DDSEN(X,Y(IPS),YPRIME(IPS),CJ,
     *           DELTA(IPS),IRES,RPAR,IPAR,RES,NY,
     *           ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),
     *           WT(IPS),CNST,SENPAR,G_RES, A_RES, ADI)
            ISENFO(1) = NP
            IF (IRES .LT. 0) GOTO 575
            CJINV = 1.0D0/CJ
            DO IQ = NYMNQ, NY-1
               IIQ = IPS + IQ
               YPRIME(IIQ) = YPRIME(IIQ) - DELTA(IIQ)
               Y(IIQ) = Y(IIQ) - CJINV*DELTA(IIQ)
               E(IIQ) = -CJINV*DELTA(IIQ)
            END DO
         END IF
C     Error test for the state variable
C
         ENORM1 = DDWNRM(NY,E(IPS),VT(IPS),RPAR,IPAR)
         IF (NPHI .EQ. 2*NEQ) THEN
            ENORM2 = DDWNRM(NY,E(IPS+NEQ),VT(IPS+NEQ),RPAR,IPAR)
            IF (ENORM2 .GT. ENORM1) ENORM1 = ENORM2
         END IF
         IF (ENORM1 .GT. ENORM) ENORM = ENORM1
         IF (CK*ENORM1 .GT. 1.0D0) THEN
            IERNEW = -2
            IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
            IWM(LNNI) =IWM(LNNI) + LOCNNI/NSTATE
            RETURN
         END IF
         IPS = IPS + NY
      END DO
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
      IWM(LNNI) =IWM(LNNI) + LOCNNI/NSTATE
C>>>>>>>>>>>>>>>>>>>>>>state variables done<<<<<<<<<<<<<<<<<<<<<<<<<<
C
      IF (NY .EQ. NEQ .OR. ISENFO(6) .LT. 0) RETURN ! no sensitivity
      M = 0
C     
C     Then compute the sensitivites
C     
C     
C     Reevaluate the residual again for the sensitivites and solution
      IF (ISENFO(2) .LT. 2) THEN
C
C     For the finite difference method, DELTA is required in GMRES
C     linear solver and must be updated.
         IRES = 0
         CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *        SENPAR,G_RES, A_RES, ADI)
      END IF
      IF (ISENFO(2) .EQ. 5) THEN
C....................matrix times vector methods.....................
C     
C     evaluate the Jacobian and DF/DP
C     
         LJACI = 1
         LJACJ = LJACI + NY + 1
         LJAC = 1 + ISENFO(4)*NY
         CALL JRADFSP (
     1        NY, X, Y, YPRIME, DELTA, G_RES, CJ,  
     2        SENWRK(LJAC), ISENWK(LJACJ), ISENWK(LJACI), 
     3        IPAR, RPAR, SENPAR, IRES, SENWRK, ISENFO(4), 
     4        IWM(LIADF))
      END IF
      IRES = 1
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      DO I = 1, NY
         SAVR(I) = DELTA(I)
      END DO
      IF (ISENFO(7) .EQ. 2) THEN
C     
C>>>>>>>>>>>>>>>>>>>>Staggered direct method<<<<<<<<<<<<<<<<<<<<<<<<<<
C     
C     
         IF (MULDEL .EQ. 1) THEN
            DO J = NY+1, NEQ
               DELTA(J) = DELTA(J) * CONFAC
            END DO
         END IF
         NLI4STAT = IWM(LNLI)
         DO I = 1, NP
            II = I*NY + 1
            IRES = 0
            CALL DSLVK (
     *           NYMNQ, Y, X, YPRIME, SAVR, DELTA(II), WT(II), WM, IWM,
     *           RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, 
     *           RSQRTN, RHOK, RPAR, IPAR, SENPAR, 0, ID, 
     *           ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
            IF (IRES .NE. 0 .OR. IERSL .NE. 0) THEN
               IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-
     *              NLI4STAT+NP-1)/NP
               IWM(LNLI) = NLI4STAT
               GO TO 580
            END IF
C     
C     Update Y, E, and YPRIME.
C     
            DO J = 0, NYMNQ-1
               JII = II + J
               Y(JII)=Y(JII)-DELTA(JII)
               E(JII)=E(JII)-DELTA(JII)
               YPRIME(JII)=YPRIME(JII)-CJ*DELTA(JII)
            END DO
         END DO
         IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-NLI4STAT+NP-1)/NP
         IWM(LNLI) = NLI4STAT
         RETURN
      END IF
C     
C>>>>>>>>>>>>>>>>>>>>>>>>>staggered corrector method<<<<<<<<<<<<<<<<<<<
C 
      SS = S
 500  CONTINUE
C     
C     If necessary, multiply residual by convergence factor.
C     
      IF (MULDEL .EQ. 1) THEN
         DO J = NY+1, NEQ
            DELTA(J) = DELTA(J) * CONFAC
         END DO
      END IF
      IRES = 0         
      NLI4STAT = IWM(LNLI)
      DO I = 1, NP
         II = I*NY + 1
         CALL DSLVK (
     *        NYMNQ, Y, X, YPRIME, SAVR, DELTA(II), WT(II), WM, IWM,
     *        RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, 
     *        RSQRTN, RHOK, RPAR, IPAR, SENPAR, 0, ID, 
     *        ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
         IF (IRES .NE. 0 .OR. IERSL .NE. 0) THEN
            IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-NLI4STAT+NP-1)/NP
            IWM(LNLI) = NLI4STAT
            GO TO 580
         END IF
C     
C     Update Y, E, and YPRIME.
C     
         DO J = 0, NYMNQ-1
            JII = II + J
            Y(JII)=Y(JII)-DELTA(JII)
            E(JII)=E(JII)-DELTA(JII)
            YPRIME(JII)=YPRIME(JII)-CJ*DELTA(JII)
         END DO
      END DO
      IWM(LNLIS) = IWM(LNLIS)+(IWM(LNLI)-NLI4STAT+NP-1)/NP
      IWM(LNLI) = NLI4STAT
      DELNRM = 0.0d0
      DO I = 1, NP
         II = I*NY + 1
         ENORMTMP = DDWNRM(NYMNQ,DELTA(II),WT(II),RPAR,IPAR)
         IF (ENORMTMP .GT. DELNRM) DELNRM = ENORMTMP 
      END DO
C     
C     Test for convergence of the iteration.
C     
      IF (M .EQ. 0) THEN
         OLDNRM = DELNRM
         IF (DELNRM .LE. TOLNEW) GO TO 570
      ELSE
         RATE = (DELNRM/OLDNRM)**(1.0D0/M)
         IF (RATE .GT. 0.9D0) GO TO 580
         SS = RATE/(1.0D0 - RATE)
      ENDIF
      IF (SS*DELNRM .LE. EPCON) GO TO 570
C     
C     The corrector has not yet converged.
C     Update M and test whether the
C     maximum number of iterations have
C     been tried.
C     
      M=M+1
      IF(M.GE.MAXIT) GO TO 580
C     
C     Evaluate the residual,
C     and go back to do another iteration.
C     
      IRES = 1
      CALL DDSEN(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,IWM(LNRE),IWM(LNSE),WT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      IF (IRES .LT. 0) GO TO 580
      GO TO 500
C     
C     The iteration has converged.
 570  CONTINUE
      RETURN     
C     
C     The iteration has not converged.  Set IERNEW appropriately.
C     
 575  CONTINUE
      IWM(LNNI) =IWM(LNNI) + LOCNNI/NSTATE                 
      IF (NSTATE .GT. 1) ISENFO(6) = -NSTATE
 580  CONTINUE
      IF (IRES .LE. -2 ) THEN
         IERNEW = -1
      ELSE
         IERNEW = 1
      ENDIF
      RETURN
C
C------END OF SUBROUTINE DNSK-------------------------------------------
      END
      SUBROUTINE DSLVK (NEQ, Y, TN, YPRIME, SAVR, X, EWT, WM, IWM,
     *   RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK,
     *   RPAR, IPAR, SENPAR, ICOPT, ID, 
     *   ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, ADI)
C
C***BEGIN PROLOGUE  DSLVK
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940928   Removed MNEWT and added RHOK in call list.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C DSLVK uses a restart algorithm and interfaces to DSPIGM for
C the solution of the linear system arising from a Newton iteration.
C
C In addition to variables described elsewhere,
C communication with DSLVK uses the following variables..
C WM    = Real work space containing data for the algorithm
C         (Krylov basis vectors, Hessenberg matrix, etc.).
C IWM   = Integer work space containing data for the algorithm.
C X     = The right-hand side vector on input, and the solution vector
C         on output, of length NEQ.
C IRES  = Error flag from RES.
C IERSL = Output flag ..
C         IERSL =  0 means no trouble occurred (or user RES routine
C                    returned IRES < 0)
C         IERSL =  1 means the iterative method failed to converge
C                    (DSPIGM returned IFLAG > 0.)
C         IERSL = -1 means there was a nonrecoverable error in the
C                    iterative solver, and an error exit will occur.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DSCAL, DCOPY, DSPIGM
C
C***END PROLOGUE  DSLVK
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER NEQ, IWM, IRES, IERSL, IPAR, ICOPT, ID
      DOUBLE PRECISION Y, TN, YPRIME, SAVR, X, EWT, WM, CJ, EPLIN,
     1   SQRTN, RSQRTN, RHOK, RPAR, SENPAR, ADI
      DIMENSION Y(*), YPRIME(*), SAVR(*), X(*), EWT(*), ADI(*), 
     1  WM(*), IWM(*), RPAR(*), IPAR(*), SENPAR(*), ID(*)
C
      INTEGER IFLAG, IRST, NRSTS, NRMAX, LR, LDL, LHES, LGMR, LQ, LV,
     1        LWK, LZ, MAXLP1, NPSL
      INTEGER NLI, NPS, NCFL, NRE, MAXL, KMP, MITER
      INTEGER ISENFO(*), ISENWK(*)
      DOUBLE PRECISION SENWRK(*), CNST
      EXTERNAL G_RES, A_RES, K_RES
      EXTERNAL  RES, PSOL
C    
      PARAMETER (LNRE=12, LNCFL=16, LNLI=20, LNPS=21) 
      PARAMETER (LLOCWP=29, LLCIWP=30)
      PARAMETER (LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26)
      PARAMETER (LMTYPE=4)
C
C-----------------------------------------------------------------------
C IRST is set to 1, to indicate restarting is in effect.
C NRMAX is the maximum number of restarts.
C-----------------------------------------------------------------------
      DATA IRST/1/
C
      LIWP = IWM(LLCIWP)
      NLI = IWM(LNLI)
      NPS = IWM(LNPS)
      NCFL = IWM(LNCFL)
      NRE = IWM(LNRE)
      LWP = IWM(LLOCWP)
      MAXL = IWM(LMAXL) 
      KMP = IWM(LKMP)
      NRMAX = IWM(LNRMAX) 
      MITER = IWM(LMITER)
      IERSL = 0
      IRES = 0
C-----------------------------------------------------------------------
C Use a restarting strategy to solve the linear system
C P*X = -F.  Parse the work vector, and perform initializations.
C Note that zero is the initial guess for X.
C-----------------------------------------------------------------------
      NY = NEQ + ISENFO(10)
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + NY*MAXL
      LHES = LR + NY + 1
      LQ = LHES + MAXL*MAXLP1
      LWK = LQ + 2*MAXL
      LDL = LWK + MIN0(1,MAXL-KMP)*NY
      LZ = LDL + NY
      CALL DSCAL (NEQ, RSQRTN, EWT, 1)
      CALL DCOPY (NEQ, X, 1, WM(LR), 1)
      DO 110 I = 1,NEQ
 110     X(I) = 0.D0
C-----------------------------------------------------------------------
C Top of loop for the restart algorithm.  Initial pass approximates
C X and sets up a transformed system to perform subsequent restarts
C to update X.  NRSTS is initialized to -1, because restarting
C does not occur until after the first pass.
C Update NRSTS; conditionally copy DL to R; call the DSPIGM
C algorithm to solve A*Z = R;  updated counters;  update X with
C the residual solution.
C Note:  if convergence is not achieved after NRMAX restarts,
C then the linear solver is considered to have failed.
C-----------------------------------------------------------------------
      NRSTS = -1
 115  CONTINUE
      NRSTS = NRSTS + 1
      IF (NRSTS .GT. 0) CALL DCOPY (NEQ, WM(LDL), 1, WM(LR),1)
      CALL DSPIGM (NEQ, TN, Y, YPRIME, SAVR, WM(LR), EWT, MAXL, MAXLP1,
     1   KMP, EPLIN, CJ, RES, IRES, NRES, PSOL, NPSL, WM(LZ), WM(LV),
     2   WM(LHES), WM(LQ), LGMR, WM(LWP), IWM(LIWP), WM(LWK),
     3   WM(LDL), RHOK, IFLAG, IRST, NRSTS, RPAR, IPAR, SENPAR,
     4   ICOPT,ID, ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, 
     5   IWM(LMTYPE), ADI,NY)
      NLI = NLI + LGMR
      NPS = NPS + NPSL
      NRE = NRE + NRES
      DO 120 I = 1,NEQ
 120     X(I) = X(I) + WM(LZ+I-1) 
      IF ((IFLAG .EQ. 1) .AND. (NRSTS .LT. NRMAX) .AND. (IRES .EQ. 0))
     1   GO TO 115
C-----------------------------------------------------------------------
C The restart scheme is finished.  Test IRES and IFLAG to see if
C convergence was not achieved, and set flags accordingly.
C-----------------------------------------------------------------------
      IF (IRES .LT. 0) THEN
         NCFL = NCFL + 1
      ELSE IF (IFLAG .NE. 0) THEN
         NCFL = NCFL + 1
         IF (IFLAG .GT. 0) IERSL = 1 
         IF (IFLAG .LT. 0) IERSL = -1 
      ENDIF
C-----------------------------------------------------------------------
C Update IWM with counters, rescale EWT, and return.
C-----------------------------------------------------------------------
      IWM(LNLI)  = NLI
      IWM(LNPS)  = NPS
      IWM(LNCFL) = NCFL
      IWM(LNRE)  = NRE
      CALL DSCAL (NEQ, SQRTN, EWT, 1)
      RETURN
C
C------END OF SUBROUTINE DSLVK------------------------------------------
      END
      SUBROUTINE DSPIGM (NEQ, TN, Y, YPRIME, SAVR, R, WGHT, MAXL,
     *   MAXLP1, KMP, EPLIN, CJ, RES, IRES, NRE, PSOL, NPSL, Z, V,
     *   HES, Q, LGMR, WP, IWP, WK, DL, RHOK, IFLAG, IRST, NRSTS,
     *   RPAR, IPAR, SENPAR,ICOPT,ID, 
     *   ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES, MTYPE,
     *   ADI,NY)
C
C***BEGIN PROLOGUE  DSPIGM
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940927   Removed MNEWT and added RHOK in call list.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine solves the linear system A * Z = R using a scaled
C preconditioned version of the generalized minimum residual method.
C An initial guess of Z = 0 is assumed.
C
C      On entry
C
C          NEQ = Problem size, passed to PSOL.
C
C           TN = Current Value of T.
C
C            Y = Array Containing current dependent variable vector.
C
C       YPRIME = Array Containing current first derivative of Y.
C
C         SAVR = Array containing current value of G(T,Y,YPRIME).
C
C            R = The right hand side of the system A*Z = R.
C                R is also used as work space when computing
C                the final approximation and will therefore be
C                destroyed.
C                (R is the same as V(*,MAXL+1) in the call to DSPIGM.)
C
C         WGHT = The vector of length NEQ containing the nonzero
C                elements of the diagonal scaling matrix.
C
C         MAXL = The maximum allowable order of the matrix H.
C
C       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C
C          KMP = The number of previous vectors the new vector, VNEW,
C                must be made orthogonal to.  (KMP .LE. MAXL.)
C
C        EPLIN = Tolerance on residuals R-A*Z in weighted rms norm.
C
C           CJ = Scalar proportional to current value of 
C                1/(step size H).
C
C           WK = Real work array used by routine DATV and PSOL.
C
C           DL = Real work array used for calculation of the residual
C                norm RHO when the method is incomplete (KMP.LT.MAXL)
C                and/or when using restarting.
C
C           WP = Real work array used by preconditioner PSOL.
C
C          IWP = Integer work array used by preconditioner PSOL.
C
C         IRST = Method flag indicating if restarting is being
C                performed.  IRST .GT. 0 means restarting is active,
C                while IRST = 0 means restarting is not being used.
C
C        NRSTS = Counter for the number of restarts on the current
C                call to DSPIGM.  If NRSTS .GT. 0, then the residual
C                R is already scaled, and so scaling of R is not
C                necessary.
C
C
C      On Return
C
C         Z    = The final computed approximation to the solution
C                of the system A*Z = R.
C
C         LGMR = The number of iterations performed and
C                the current order of the upper Hessenberg
C                matrix HES.
C
C         NRE  = The number of calls to RES (i.e. DATV)
C
C         NPSL = The number of calls to PSOL.
C
C         V    = The neq by (LGMR+1) array containing the LGMR
C                orthogonal vectors V(*,1) to V(*,LGMR).
C
C         HES  = The upper triangular factor of the QR decomposition
C                of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C                entries are the scaled inner-products of A*V(*,I)
C                and V(*,K).
C
C         Q    = Real array of length 2*MAXL containing the components
C                of the givens rotations used in the QR decomposition
C                of HES.  It is loaded in DHEQR and used in DHELS.
C
C         IRES = Error flag from RES.
C
C           DL = Scaled preconditioned residual, 
C                (D-inverse)*(P-inverse)*(R-A*Z). Only loaded when
C                performing restarts of the Krylov iteration.
C
C         RHOK = Weighted norm of final preconditioned residual.
C
C        IFLAG = Integer error flag..
C                0 Means convergence in LGMR iterations, LGMR.LE.MAXL.
C                1 Means the convergence test did not pass in MAXL
C                  iterations, but the new residual norm (RHO) is
C                  .LT. the old residual norm (RNRM), and so Z is
C                  computed.
C                2 Means the convergence test did not pass in MAXL
C                  iterations, new residual norm (RHO) .GE. old residual
C                  norm (RNRM), and the initial guess, Z = 0, is
C                  returned.
C                3 Means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 Means there was an unrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   PSOL, DNRM2, DSCAL, DATV, DORTH, DHEQR, DCOPY, DHELS, DAXPY
C
C***END PROLOGUE  DSPIGM
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER NEQ,MAXL,MAXLP1,KMP,IRES,NRE,NPSL,LGMR,IWP,
     1   IFLAG,IRST,NRSTS,IPAR, ICOPT, ID, MTYPE
      DOUBLE PRECISION TN,Y,YPRIME,SAVR,R,WGHT,EPLIN,CJ,Z,V,HES,Q,WP,WK,
     1   DL,RHOK,RPAR, SENPAR, ADI
      DIMENSION Y(*), YPRIME(*), SAVR(*), R(*), WGHT(*), Z(*),
     1   V(NY,*), HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*),
     2   RPAR(*), IPAR(*), SENPAR(*), ID(*), ADI(*)
      INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
      DOUBLE PRECISION RNRM,C,DLNRM,PROD,RHO,S,SNORMW,DNRM2,TEM
      INTEGER ISENFO(*), ISENWK(*)
      DOUBLE PRECISION SENWRK(*), CNST
      EXTERNAL G_RES, A_RES, K_RES
      EXTERNAL  RES, PSOL
C
      IER = 0
      IFLAG = 0
      LGMR = 0
      NPSL = 0
      NRE = 0
C-----------------------------------------------------------------------
C The initial guess for Z is 0.  The initial residual is therefore
C the vector R.  Initialize Z to 0.
C-----------------------------------------------------------------------
      DO 10 I = 1,NEQ
 10     Z(I) = 0.0D0
C-----------------------------------------------------------------------
C Apply inverse of left preconditioner to vector R if NRSTS .EQ. 0.
C Form V(*,1), the scaled preconditioned right hand side.
C-----------------------------------------------------------------------
      IF (NRSTS .EQ. 0) THEN
         CALL PSOL (NEQ, TN, Y, YPRIME, SAVR, WK, CJ, WGHT, WP, IWP,
     1      R, EPLIN, IER, RPAR, IPAR, SENPAR)
         NPSL = 1
         IF (IER .NE. 0) GO TO 300
         DO 30 I = 1,NEQ
 30         V(I,1) = R(I)*WGHT(I)
      ELSE
         DO 35 I = 1,NEQ
 35         V(I,1) = R(I)
      ENDIF
C-----------------------------------------------------------------------
C Calculate norm of scaled vector V(*,1) and normalize it
C If, however, the norm of V(*,1) (i.e. the norm of the preconditioned
C residual) is .le. EPLIN, then return with Z=0.
C-----------------------------------------------------------------------
      RNRM = DNRM2 (NEQ, V, 1)
      IF (RNRM .LE. EPLIN) THEN
        RHOK = RNRM
        RETURN
      ENDIF
      TEM = 1.0D0/RNRM
      CALL DSCAL (NEQ, TEM, V(1,1), 1)
C-----------------------------------------------------------------------
C Zero out the HES array.
C-----------------------------------------------------------------------
      DO 65 J = 1,MAXL
        DO 60 I = 1,MAXLP1
 60       HES(I,J) = 0.0D0
 65     CONTINUE
C-----------------------------------------------------------------------
C Main loop to compute the vectors V(*,2) to V(*,MAXL).
C The running product PROD is needed for the convergence test.
C-----------------------------------------------------------------------
      PROD = 1.0D0
      
      DO 90 LL = 1,MAXL
         LGMR = LL
C-----------------------------------------------------------------------
C Call routine DATV to compute VNEW = ABAR*V(LL), where ABAR is
C the matrix A with scaling and inverse preconditioner factors applied.
C Call routine DORTH to orthogonalize the new vector VNEW = V(*,LL+1).
C call routine DHEQR to update the factors of HES.
C-----------------------------------------------------------------------
         CALL DATV (NEQ, Y, TN, YPRIME, SAVR, V(1,LL), WGHT, Z,
     1        RES, IRES, PSOL, V(1,LL+1), WK, WP, IWP, CJ, EPLIN,
     1        IER, NRE, NPSL, RPAR, IPAR, SENPAR, ICOPT, ID, 
     1        ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES,
     1        MTYPE,ADI)
         IF (IRES .LT. 0) RETURN
         IF (IER .NE. 0) GO TO 300
         CALL DORTH (V(1,LL+1), V, HES, NEQ, LL, MAXLP1, KMP,SNORMW,NY)
         HES(LL+1,LL) = SNORMW
         CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL)
         IF (INFO .EQ. LL) GO TO 120
C-----------------------------------------------------------------------
C Update RHO, the estimate of the norm of the residual R - A*ZL.
C If KMP .LT. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C necessarily orthogonal for LL .GT. KMP.  The vector DL must then
C be computed, and its norm used in the calculation of RHO.
C-----------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*RNRM)
        IF ((LL.GT.KMP) .AND. (KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              CALL DCOPY (NEQ, V(1,1), 1, DL, 1)
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,NEQ
 70                 DL(K) = S*DL(K) + C*V(K,IP1)
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,NEQ
 80        DL(K) = S*DL(K) + C*V(K,LLP1)
           DLNRM = DNRM2 (NEQ, DL, 1)
           RHO = RHO*DLNRM
        ENDIF
C-----------------------------------------------------------------------
C Test for convergence.  If passed, compute approximation ZL.
C If failed and LL .LT. MAXL, then continue iterating.
C-----------------------------------------------------------------------
        IF (RHO .LE. EPLIN) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C-----------------------------------------------------------------------
C Rescale so that the norm of V(1,LL+1) is one.
C-----------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL (NEQ, TEM, V(1,LL+1), 1)
 90     CONTINUE
 100  CONTINUE
      IF (RHO .LT. RNRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
      DO 130 I = 1,NEQ
 130     Z(I) = 0.D0
      RETURN
 150  IFLAG = 1
C-----------------------------------------------------------------------
C The tolerance was not met, but the residual norm was reduced.
C If performing restarting (IRST .gt. 0) calculate the residual vector
C RL and store it in the DL array.  If the incomplete version is 
C being used (KMP .lt. MAXL) then DL has already been calculated.
C-----------------------------------------------------------------------
      IF (IRST .GT. 0) THEN
         IF (KMP .EQ. MAXL) THEN
C
C           Calculate DL from the V(I)'s.
C
            CALL DCOPY (NEQ, V(1,1), 1, DL, 1)
            MAXLM1 = MAXL - 1
            DO 175 I = 1,MAXLM1
               IP1 = I + 1
               I2 = I*2
               S = Q(I2)
               C = Q(I2-1)
               DO 170 K = 1,NEQ
 170              DL(K) = S*DL(K) + C*V(K,IP1)
 175        CONTINUE
            S = Q(2*MAXL)
            C = Q(2*MAXL-1)/SNORMW
            DO 180 K = 1,NEQ
 180           DL(K) = S*DL(K) + C*V(K,MAXLP1)
         ENDIF
C
C        Scale DL by RNRM*PROD to obtain the residual RL.
C
         TEM = RNRM*PROD
         CALL DSCAL(NEQ, TEM, DL, 1)
      ENDIF
C-----------------------------------------------------------------------
C Compute the approximation ZL to the solution.
C Since the vector Z was used as work space, and the initial guess
C of the Newton correction is zero, Z must be reset to zero.
C-----------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
 210    R(K) = 0.0D0
      R(1) = RNRM
      CALL DHELS (HES, MAXLP1, LL, Q, R)
      DO 220 K = 1,NEQ
 220    Z(K) = 0.0D0
      DO 230 I = 1,LL
        CALL DAXPY (NEQ, R(I), V(1,I), 1, Z, 1)
 230    CONTINUE
      DO 240 I = 1,NEQ
 240    Z(I) = Z(I)/WGHT(I)
C Load RHO into RHOK.
      RHOK = RHO
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns forced by routine PSOL.
C-----------------------------------------------------------------------
 300  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
C
      RETURN
C
C------END OF SUBROUTINE DSPIGM-----------------------------------------
      END
      SUBROUTINE DATV (NEQ, Y, TN, YPRIME, SAVR, V, WGHT, YPTEM, RES,
     *   IRES, PSOL, Z, VTEM, WP, IWP, CJ, EPLIN, IER, NRE, NPSL,
     *   RPAR, IPAR, SENPAR, ICOPT, ID, 
     *   ISENFO,SENWRK,ISENWK, CNST, G_RES, A_RES, K_RES,
     *   MTYPE,ADI)
C
C***BEGIN PROLOGUE  DATV
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine computes the product
C
C   Z = (D-inverse)*(P-inverse)*(dF/dY)*(D*V),
C
C where F(Y) = G(T, Y, CJ*(Y-A)), CJ is a scalar proportional to 1/H,
C and A involves the past history of Y.  The quantity CJ*(Y-A) is
C an approximation to the first derivative of Y and is stored
C in the array YPRIME.  Note that dF/dY = dG/dY + CJ*dG/dYPRIME.
C
C D is a diagonal scaling matrix, and P is the left preconditioning
C matrix.  V is assumed to have L2 norm equal to 1.
C The product is stored in Z and is computed by means of a
C difference quotient, a call to RES, and one call to PSOL.
C
C      On entry
C
C          NEQ = Problem size, passed to RES and PSOL.
C
C            Y = Array containing current dependent variable vector.
C
C       YPRIME = Array containing current first derivative of y.
C
C         SAVR = Array containing current value of G(T,Y,YPRIME).
C
C            V = Real array of length NEQ (can be the same array as Z).
C
C         WGHT = Array of length NEQ containing scale factors.
C                1/WGHT(I) are the diagonal elements of the matrix D.
C
C        YPTEM = Work array of length NEQ.
C
C         VTEM = Work array of length NEQ used to store the
C                unscaled version of V.
C
C         WP = Real work array used by preconditioner PSOL.
C
C         IWP = Integer work array used by preconditioner PSOL.
C
C           CJ = Scalar proportional to current value of 
C                1/(step size H).
C
C
C      On return
C
C            Z = Array of length NEQ containing desired scaled
C                matrix-vector product.
C
C         IRES = Error flag from RES.
C
C          IER = Error flag from PSOL.
C
C         NRE  = The number of calls to RES.
C
C         NPSL = The number of calls to PSOL.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   AMUVSP, RES, PSOL
C
C***END PROLOGUE  DATV
C
      IMPLICIT NONE
      INTEGER NEQ, IRES, IWP, IER, NRE, NPSL, IPAR, ICOPT, ID, MTYPE
      DOUBLE PRECISION Y, TN, YPRIME, SAVR, V, WGHT, YPTEM, Z, VTEM,
     1   WP, CJ, RPAR, SENPAR, ADI
      DIMENSION Y(*), YPRIME(*), SAVR(*), V(*), WGHT(*), YPTEM(*),
     1   Z(*), VTEM(*), WP(*), IWP(*), RPAR(*), IPAR(*), SENPAR(*),
     2   ID(*), ADI(*)
      INTEGER ISENFO(*), ISENWK(*)
      DOUBLE PRECISION SENWRK(*), CNST
      EXTERNAL G_RES, A_RES, K_RES
      INTEGER I, NSE, NY
      DOUBLE PRECISION EPLIN
      
      EXTERNAL  RES, PSOL
C
      IRES = 0
C-----------------------------------------------------------------------
C Set VTEM = D * V.
C-----------------------------------------------------------------------
      DO 10 I = 1,NEQ
 10     VTEM(I) = V(I)/WGHT(I)
      IER = 0
C======================Matrix vector product option for MV=============
c$$$      METHOD_MV = 0
c$$$      IF (METHOD_MV .EQ. 1) THEN
c$$$C---------------set up pointer for IWP---------------------------------
c$$$         N2EQP3 = 2*NEQ + 3
c$$$         CALL AMUVSP(NEQ,VTEM,Z,WP,IWP(N2EQP3),IWP)
c$$$         GOTO 80
c$$$      END IF
      IF (MTYPE .GT. 0) THEN
         IF (MTYPE .EQ. 1) THEN
            CALL K_RES(TN,Y,YPRIME,CJ,IRES,RPAR,IPAR,SENPAR,VTEM,Z)
         ELSE IF (MTYPE .EQ. 2) THEN
            NY = NEQ + ISENFO(10)
            CALL KMVBYADF(TN,Y,YPRIME,CJ,SAVR,VTEM,Z,IRES,RPAR,IPAR,
     *           SENPAR,K_RES,NY,SENWRK)
         ELSE IF (MTYPE .EQ. 3) THEN
C
C     Matrix-not-free matrix-vector product. The matrix information is
C     saved in WP and IWP
            CALL K_RES(NEQ, WP, IWP, VTEM, Z)      
         END IF
         GOTO 80
      END IF
C=======================for adjoint DAE ===============================
      IF (ISENFO(6).LT.0.AND.ICOPT.EQ.0) THEN
         DO I = 1, NEQ
            YPTEM(I) = VTEM(I)*CJ
         END DO
         NY = NEQ + ISENFO(10)
         CALL KMVAD(TN,VTEM,YPTEM,CJ,Z,IRES,RPAR,IPAR,RES,NY,
     *        ISENFO,SENWRK,ISENWK,SENPAR,G_RES,ADI)
         GOTO 80
      END IF      
C=======================Finite difference method for AV================
C
C-----------------------------------------------------------------------
C Store Y in Z and increment Z by VTEM.
C Store YPRIME in YPTEM and increment YPTEM by VTEM*CJ.
C-----------------------------------------------------------------------
      IF (ICOPT .EQ. 0) THEN
         DO I = 1,NEQ
            YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
            Z(I) = Y(I) + VTEM(I)
         END DO
      ELSE IF (ICOPT.EQ.1) THEN
         DO I = 1,NEQ
            IF (ID(I) .GT. 0) THEN
               Z(I) = Y(I)
               YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
            ELSE 
               Z(I) = Y(I) + VTEM(I)
               YPTEM(I) = YPRIME(I)
            END IF
         END DO
      ELSE IF (ICOPT .EQ. 5) THEN
         DO I = 1,NEQ
            IF (ID(I) .GT. -2) THEN
               Z(I) = Y(I)
               YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
            ELSE 
               Z(I) = Y(I) + VTEM(I)
               YPTEM(I) = YPRIME(I)
            END IF
         END DO         
      ELSE IF (ICOPT.EQ.2) THEN
         DO I = 1,NEQ
            Z(I) = Y(I) + VTEM(I)
            YPTEM(I) = YPRIME(I)
         END DO
      ELSE IF (ICOPT.EQ.3 .OR. ICOPT.EQ.4) THEN
         DO I = 1,NEQ
            IF (ID(I) .EQ. 1) THEN
               Z(I) = Y(I) + VTEM(I) 
               YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
            ELSE IF (ID(I) .EQ. 2) THEN
               Z(I) = Y(I)
               YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
            ELSE IF (ID(I) .EQ. 3 .OR. ID(I) .LT. 0) THEN
               Z(I) = Y(I) + VTEM(I) 
               YPTEM(I) = YPRIME(I)              
            END IF
         END DO
      ELSE IF (ICOPT .EQ. 6) THEN
         DO I = 1,NEQ
            Z(I) = Y(I) + VTEM(I)
            YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
         END DO         
      END IF
C-----------------------------------------------------------------------
C Call RES with incremented Y, YPRIME arguments
C stored in Z, YPTEM.  VTEM is overwritten with new residual.
C-----------------------------------------------------------------------
      CONTINUE
C      CALL RES(TN,Z,YPTEM,CJ,VTEM,IRES,RPAR,IPAR,SENPAR)
C      NRE = NRE + 1
      NY = NEQ + ISENFO(10)
      IF (ISENFO(6).LT.0 .AND. ISENFO(10)-ISENFO(4).GT.0) THEN
         DO I = 1, ISENFO(10)-ISENFO(4)
            Z(NEQ+I) = Y(NEQ+I)
            YPTEM(NEQ+I) = YPRIME(NEQ+I)
         END DO
      END IF
      CALL DDSEN(TN,Z,YPTEM,CJ,VTEM,IRES,RPAR,IPAR,RES,NY,
     *     ISENFO,SENWRK,ISENWK,NRE,NSE,WGHT,CNST,
     *     SENPAR,G_RES, A_RES, ADI)
      IF (IRES .LT. 0) RETURN
C-----------------------------------------------------------------------
C Set Z = (dF/dY) * VBAR using difference quotient.
C (VBAR is old value of VTEM before calling RES)
C-----------------------------------------------------------------------
      DO 70 I = 1,NEQ
 70     Z(I) = VTEM(I) - SAVR(I)
C
 80   CONTINUE
C-----------------------------------------------------------------------
C Apply inverse of left preconditioner to Z.
C-----------------------------------------------------------------------
      CALL PSOL (NEQ, TN, Y, YPRIME, SAVR, YPTEM, CJ, WGHT, WP, IWP,
     1   Z, EPLIN, IER, RPAR, IPAR, SENPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) RETURN
C-----------------------------------------------------------------------
C Apply D-inverse to Z and return.
C-----------------------------------------------------------------------
      DO 90 I = 1,NEQ
 90     Z(I) = Z(I)*WGHT(I)
      RETURN
C
C------END OF SUBROUTINE DATV-------------------------------------------
      END
      SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW, LDV)
C
C***BEGIN PROLOGUE  DORTH
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine orthogonalizes the vector VNEW against the previous
C KMP vectors in the V array.  It uses a modified Gram-Schmidt
C orthogonalization procedure with conditional reorthogonalization.
C
C      On entry
C
C         VNEW = The vector of length N containing a scaled product
C                OF The Jacobian and the vector V(*,LL).
C
C         V    = The LDV x LL array containing the previous LL
C                orthogonal vectors V(*,1) to V(*,LL).
C
C         HES  = An LL x LL upper Hessenberg matrix containing,
C                in HES(I,K), K.LT.LL, scaled inner products of
C                A*V(*,K) and V(*,I).
C
C        LDHES = The leading dimension of the HES array.
C
C         N    = The order of the matrix A, and the length of VNEW.
C
C         LL   = The current order of the matrix HES.
C
C          KMP = The number of previous vectors the new vector VNEW
C                must be made orthogonal to (KMP .LE. MAXL).
C
C
C      On return
C
C         VNEW = The new vector orthogonal to V(*,I0),
C                where I0 = MAX(1, LL-KMP+1).
C
C         HES  = Upper Hessenberg matrix with column LL filled in with
C                scaled inner products of A*V(*,LL) and V(*,I).
C
C       SNORMW = L-2 norm of VNEW.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DDOT, DNRM2, DAXPY 
C
C***END PROLOGUE  DORTH
C
      INTEGER N, LL, LDHES, KMP
      DOUBLE PRECISION VNEW, V, HES, SNORMW
      DIMENSION VNEW(*), V(LDV,*), HES(LDHES,*)
      INTEGER I, I0
      DOUBLE PRECISION ARG, DDOT, DNRM2, SUMDSQ, TEM, VNRM
C
C-----------------------------------------------------------------------
C Get norm of unaltered VNEW for later use.
C-----------------------------------------------------------------------
      VNRM = DNRM2 (N, VNEW, 1)
C-----------------------------------------------------------------------
C Do Modified Gram-Schmidt on VNEW = A*V(LL).
C Scaled inner products give new column of HES.
C Projections of earlier vectors are subtracted from VNEW.
C-----------------------------------------------------------------------
      I0 = MAX0(1,LL-KMP+1)
      DO 10 I = I0,LL
        HES(I,LL) = DDOT (N, V(1,I), 1, VNEW, 1)
        TEM = -HES(I,LL)
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
 10     CONTINUE
C-----------------------------------------------------------------------
C Compute SNORMW = norm of VNEW.
C If VNEW is small compared to its input value (in norm), then
C Reorthogonalize VNEW to V(*,1) through V(*,LL).
C Correct if relative correction exceeds 1000*(unit roundoff).
C Finally, correct SNORMW using the dot products involved.
C-----------------------------------------------------------------------
      SNORMW = DNRM2 (N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0D0
      DO 30 I = I0,LL
        TEM = -DDOT (N, V(1,I), 1, VNEW, 1)
        IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
        HES(I,LL) = HES(I,LL) - TEM
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
        SUMDSQ = SUMDSQ + TEM**2
 30     CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
      RETURN
C
C------END OF SUBROUTINE DORTH------------------------------------------
      END
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
C
C***BEGIN PROLOGUE  DHEQR
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine performs a QR decomposition of an upper
C     Hessenberg matrix A.  There are two options available:
C
C          (1)  performing a fresh decomposition
C          (2)  updating the QR factors by adding a row and A
C               column to the matrix A.
C
C     DHEQR decomposes an upper Hessenberg matrix by using Givens
C     rotations.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                The matrix to be decomposed.
C
C        LDA     INTEGER
C                The leading dimension of the array A.
C
C        N       INTEGER
C                A is an (N+1) by N Hessenberg matrix.
C
C        IJOB    INTEGER
C                = 1     Means that a fresh decomposition of the
C                        matrix A is desired.
C                .GE. 2  Means that the current decomposition of A
C                        will be updated by the addition of a row
C                        and a column.
C     On return
C
C        A       The upper triangular matrix R.
C                The factorization can be written Q*A = R, where
C                Q is a product of Givens rotations and R is upper
C                triangular.
C
C        Q       DOUBLE PRECISION(2*N)
C                The factors C and S of each Givens rotation used
C                in decomposing A.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  If  A(K,K) .EQ. 0.0.  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHELS will divide by zero
C                     if called.
C
C     Modification of LINPACK.
C     Peter Brown, Lawrence Livermore Natl. Lab.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED (NONE)
C
C***END PROLOGUE  DHEQR
C
      INTEGER LDA, N, INFO, IJOB
      DOUBLE PRECISION A(LDA,*), Q(*)
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      DOUBLE PRECISION C, S, T, T1, T2
C
      IF (IJOB .GT. 1) GO TO 70
C-----------------------------------------------------------------------
C A new factorization is desired.
C-----------------------------------------------------------------------
C
C     QR decomposition without pivoting.
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute Kth column of R.
C           First, multiply the Kth column of A by the previous
C           K-1 Givens rotations.
C
            IF (KM1 .LT. 1) GO TO 20
            DO 10 J = 1, KM1
              I = 2*(J-1) + 1
              T1 = A(J,K)
              T2 = A(J+1,K)
              C = Q(I)
              S = Q(I+1)
              A(J,K) = C*T1 - S*T2
              A(J+1,K) = S*T1 + C*T2
   10         CONTINUE
C
C           Compute Givens components C and S.
C
   20       CONTINUE
            IQ = 2*KM1 + 1
            T1 = A(K,K)
            T2 = A(KP1,K)
            IF (T2 .NE. 0.0D0) GO TO 30
              C = 1.0D0
              S = 0.0D0
              GO TO 50
   30       CONTINUE
            IF (ABS(T2) .LT. ABS(T1)) GO TO 40
              T = T1/T2
              S = -1.0D0/SQRT(1.0D0+T*T)
              C = -S*T
              GO TO 50
   40       CONTINUE
              T = T2/T1
              C = 1.0D0/SQRT(1.0D0+T*T)
              S = -C*T
   50       CONTINUE
            Q(IQ) = C
            Q(IQ+1) = S
            A(K,K) = C*T1 - S*T2
            IF (A(K,K) .EQ. 0.0D0) INFO = K
   60 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C The old factorization of A will be updated.  A row and a column
C has been added to the matrix A.
C N by N-1 is now the old size of the matrix.
C-----------------------------------------------------------------------
  70  CONTINUE
      NM1 = N - 1
C-----------------------------------------------------------------------
C Multiply the new column by the N previous Givens rotations.
C-----------------------------------------------------------------------
      DO 100 K = 1,NM1
        I = 2*(K-1) + 1
        T1 = A(K,N)
        T2 = A(K+1,N)
        C = Q(I)
        S = Q(I+1)
        A(K,N) = C*T1 - S*T2
        A(K+1,N) = S*T1 + C*T2
 100    CONTINUE
C-----------------------------------------------------------------------
C Complete update of decomposition by forming last Givens rotation,
C and multiplying it times the column vector (A(N,N),A(NP1,N)).
C-----------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF (T2 .NE. 0.0D0) GO TO 110
        C = 1.0D0
        S = 0.0D0
        GO TO 130
 110  CONTINUE
      IF (ABS(T2) .LT. ABS(T1)) GO TO 120
        T = T1/T2
        S = -1.0D0/SQRT(1.0D0+T*T)
        C = -S*T
        GO TO 130
 120  CONTINUE
        T = T2/T1
        C = 1.0D0/SQRT(1.0D0+T*T)
        S = -C*T
 130  CONTINUE
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C
C------END OF SUBROUTINE DHEQR------------------------------------------
      END
      SUBROUTINE DHELS (A, LDA, N, Q, B)
C
C***BEGIN PROLOGUE  DHELS
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This is similar to the LINPACK routine DGESL except that
C A is an upper Hessenberg matrix.
C
C     DHELS solves the least squares problem
C
C           MIN (B-A*X,B-A*X)
C
C     using the factors computed by DHEQR.
C
C     On entry
C
C        A       DOUBLE PRECISION (LDA, N)
C                The output from DHEQR which contains the upper
C                triangular factor R in the QR decomposition of A.
C
C        LDA     INTEGER
C                The leading dimension of the array  A .
C
C        N       INTEGER
C                A is originally an (N+1) by N matrix.
C
C        Q       DOUBLE PRECISION(2*N)
C                The coefficients of the N givens rotations
C                used in the QR factorization of A.
C
C        B       DOUBLE PRECISION(N+1)
C                The right hand side vector.
C
C
C     On return
C
C        B       The solution vector X.
C
C
C     Modification of LINPACK.
C     Peter Brown, Lawrence Livermore Natl. Lab.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DAXPY 
C
C***END PROLOGUE  DHELS
C
      INTEGER LDA, N
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
      INTEGER IQ, K, KB, KP1
      DOUBLE PRECISION C, S, T, T1, T2
C
C        Minimize (B-A*X,B-A*X).
C        First form Q*B.
C
         DO 20 K = 1, N
            KP1 = K + 1
            IQ = 2*(K-1) + 1
            C = Q(IQ)
            S = Q(IQ+1)
            T1 = B(K)
            T2 = B(KP1)
            B(K)   = C*T1 - S*T2
            B(KP1) = S*T1 + C*T2
   20    CONTINUE
C
C        Now solve R*X = Q*B.
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C
C------END OF SUBROUTINE DHELS------------------------------------------
      END
      SUBROUTINE DDSEN(T,Y,YP,CJ,DELTA,IRESI,RPAR,IPAR,RES,NY,ISENFO,
     *     SENWRK,ISENWK,NRE,NSE,WT,CNST,SENPAR,G_RES, A_RES, ADI)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YP(*),DELTA(*),ISENWK(*),SENWRK(*),RPAR(*),IPAR(*),
     *          ISENFO(*),WT(*),SENPAR(*), ADI(*)
C-------------------------------------------------------------------*
C     This routine approximates the solution to the sensitivity equations
C     for NP parameters using either a forward or centered finite
C     difference method (specified via ISENFO(2)). The number of calls to
C     the RES routine with IRES = 0 or 1 (number of evaluations of the
C     state system) is counted via the constant NRE.  The number of
C     evaluations of the sensitivity residuals, if they are provided
C     by the user, is counted via the constant NSE.  The vector
C     SENWRK is used as work space for the specified finite difference
C     scheme.
C-------------------------------------------------------------------*
      EXTERNAL RES, G_RES, A_RES
      IF (ISENFO(6) .LT. 0) THEN ! Adjoint equation evaluations
         CALL DDRESAD(
     *        T,Y,YP,CJ,DELTA,IRESI,RPAR,IPAR,SENPAR,
     *        RES, NY, ISENFO, SENWRK, ISENWK, G_RES, A_RES, ADI)
         NRE = NRE + 1
         RETURN
      END IF
C     
C     If we are not doing sensitivity, just call RES and return
C
      IRES = 0
      IF(ISENFO(1).EQ.0 .OR. IRESI.EQ.0) THEN
         CALL RES(T,Y,YP,CJ,DELTA,IRESI,RPAR,IPAR,SENPAR)
         IF (IRESI .NE. 3) NRE = NRE + 1
         RETURN
      ELSE IF(ISENFO(2).EQ.2) THEN
C
C     If we are doing sensitivity with analytic residuals, call
C     RES and return
C
         CALL RES(T,Y,YP,CJ,DELTA,IRESI,RPAR,IPAR,SENPAR) ! IRES is input
         IF (IRESI .NE. 3) THEN
            NSE = NSE + 1
            NRE = NRE + 1
         END IF
         RETURN
      ELSE IF(ISENFO(2).EQ.3) THEN
C
C     evaluate the residuals by the ADIFOR with seed matrix options
C
         IF (IRESI .EQ. 3) IRES = 3
         CALL rAdfSM(
     *        T,Y,YP,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,G_RES,NY,ISENFO,
     *        SENWRK)
         IF (IRESI .NE. 3) THEN
            NSE = NSE + 1
            NRE = NRE + 1
         END IF
         RETURN
      ELSE IF(ISENFO(2).EQ.4) THEN
C
C     evaluate the residuals by the ADIFOR with matrix-vector product options
C
         IF (IRESI .EQ. 3) IRES = 3
         CALL rAdfMV(
     *        T,Y,YP,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,G_RES,NY,ISENFO,
     *        SENWRK)
         IF (IRESI .NE. 3) THEN
            NSE = NSE + 1
            NRE = NRE + 1
         END IF
         RETURN
      ELSE IF (ISENFO(2) .EQ. 5) THEN
C
C     evaluate the residuals by matrix times vector options
C     df/dp saved in senwrk(1:np*ny), jac saved in senwrk(np*ny+1..)
C     ia(*) saved in isenwk(1:ny+1), ja(*) saved in isenwk(ny+2:..)
C
         CALL RMATV(NY,T,Y,YP,CJ,DELTA,RPAR,IPAR,SENPAR,ISENFO,
     *        SENWRK, ISENWK)
         NSE = NSE + 1
         RETURN
      END IF
C
C     If ISENFO(2) = 0 or 1, compute sensitivity residuals via
C     forward or centered finite differences, respectively
C
C     Initialize work space markers.
C   
      MYID = ISENFO(8)
      NUMPROCS = ISENFO(9)
      M1 = NY + NY
      M2 = M1 + NY
      NSE = NSE + 1
C
      NP = ISENFO(1)
C
C     Branch to different finite-difference schemes.
C
      IF (IRESI .EQ. 3) IRES = 3
      IF (IRES .EQ. 3) THEN
         IF (ISENFO(2).EQ.1) 
     *        CALL RES(T,Y,YP,CJ,DELTA,IRES,RPAR,IPAR,SENPAR)
         IF(IRES.LT.0) RETURN
         IBEG = NY - ISENFO(10) + 1
         IEND = NY
      ELSE
         IBEG = 1
         IEND = NY - ISENFO(10)
      END IF
C         
C
      IF(ISENFO(2).EQ.0) GO TO 200
C
C     Single difference O(1) finite difference approximation.
C
 100  CONTINUE
C
C     Iterate on the number of parameters, NP.
C
      DO IP = 1, NP
         MYIP = MYID+1 + (IP-1)*NUMPROCS
         ISROW = IP*NY
C
C       Determine the perturbation.
C
         VNORM = 0.0D0
         DO J = IBEG, IEND
            VNORM = VNORM + (WT(J)/WT(ISROW + J))**2
         ENDDO
         VNORM = DSQRT(VNORM)
         IF (MYIP .GT. ISENFO(4)) THEN
            DEL = CNST/VNORM
         ELSE
            DEL = CNST*MAX(ABS(SENPAR(MYIP)),1.0D0/VNORM )
         END IF
C
C       Save and perturb Y, YP and the parameter.
C
         IF (MYIP .LE. ISENFO(4)) THEN
            PSAVE = SENPAR(MYIP)
            SENPAR(MYIP) = PSAVE + DEL
         END IF
C
         DO I = 1, NY
            MRKR = ISROW + I
            SENWRK(I)      = Y(I)  + DEL*Y(MRKR)
            SENWRK(NY + I) = YP(I) + DEL*YP(MRKR)
         ENDDO
C
C     Call RES with perturbed values.
C
         CALL RES(T,SENWRK,SENWRK(NY+1),CJ,
     *        SENWRK(M1+1),IRES,RPAR,IPAR,SENPAR)
         IF(IRES.LT.0) RETURN
C
C     Approx. first order sensitivity, and restore Y, YP, parameter values.
C
         DELINV = 1.0D0/DEL
         DO I = IBEG, IEND
            DELTA(ISROW+I)=DELINV*(SENWRK(M1+I)-DELTA(I))
         ENDDO
C
         IF (MYIP .LE. ISENFO(4)) SENPAR(MYIP) = PSAVE
      ENDDO
C
C     First order forward finite difference approximation complete.
C
      RETURN
C
C     Second order centered finite difference scheme.
C
 200  CONTINUE
C
C     Iterate on the number of parameters, NP.
C
      DO IP = 1, NP
         MYIP = MYID+1 + (IP-1)*NUMPROCS
         ISROW = IP*NY
C 
C       Determine the perturbation.  
C
         VNORM = 0.0D0
         DO J = IBEG, IEND
            VNORM = VNORM + (WT(J)/WT(ISROW + J))**2
         ENDDO
         VNORM = DSQRT(VNORM)
         IF (MYIP .GT. ISENFO(4)) THEN
            DEL = CNST/VNORM
         ELSE
            DEL = CNST*MAX( ABS(SENPAR(MYIP)),1.0D0/VNORM )
         END IF
C
C       Save and perturb Y, YP and the parameter in the forward (+)      
C       direction.
C
         IF (MYIP .LE. ISENFO(4)) THEN
            PSAVE = SENPAR(MYIP)
            SENPAR(MYIP) = PSAVE + DEL
         END IF
C
         DO I = 1, NY
            MRKR = ISROW + I
            SENWRK(I)      = Y(I)  + DEL*Y(MRKR)
            SENWRK(NY + I) = YP(I) + DEL*YP(MRKR)
         ENDDO
C
C     Call RES with perturbed values.
C
         CALL RES(
     *        T,SENWRK,SENWRK(NY+1),CJ,SENWRK(M1 + 1),
     *        IRES,RPAR,IPAR,SENPAR)
C
         IF(IRES.LT.0)RETURN
C
C     Save and perturb Y, YP and the parameter in the backward (-)      
C     direction.
C 
         IF (MYIP .LE. ISENFO(4)) SENPAR(MYIP) = PSAVE - DEL
         DO I = 1, NY
            MRKR = ISROW + I
            SENWRK(I)      = Y(I)  - DEL*Y(MRKR)
            SENWRK(NY + I) = YP(I) - DEL*YP(MRKR)
         ENDDO
C
C     Call RES with perturbed values.
C
         CALL RES(T,SENWRK,SENWRK(NY+1),CJ,
     *        SENWRK(M2 + 1),IRES,RPAR,IPAR,SENPAR)
C
         IF(IRES.LT.0)RETURN
C
C     Approx. second order sensitivity, and restore Y, YP, parameter values.
C
         DELINV = 1.0D0/(2.0D0*DEL)
         DO I = IBEG, IEND
            DELTA(ISROW+I)=DELINV*(SENWRK(M1+I)-SENWRK(M2+I))
         ENDDO
C
         IF (MYIP .LE. ISENFO(4)) SENPAR(MYIP) = PSAVE
      ENDDO
C
C     Second order centered finite difference approximation complete.
C
      RETURN
C-------------------------END OF DDSEN ROUTINE----------------------
      END
      SUBROUTINE RMATV(
     *     NY,T,Y,YPRIME,CJ,DELTA,RPAR,IPAR,SENPAR,ISENFO,
     *     SENWRK,ISENWK)
C====================================================================
C  
C     matrix times vector method to evaluate the residuals
C
C====================================================================
      IMPLICIT NONE
      DOUBLE PRECISION T, Y(*), YPRIME(*), CJ(*), DELTA(*), RPAR(*), 
     *     SENPAR(*), SENWRK(*)
      INTEGER NY, IPAR(*), ISENFO(*), ISENWK(*), I, IP, IPOS, II, 
     *     LJAC, LJACI, LJACJ, myid, numprocs, myip
      LJAC = 1 + ISENFO(4)*NY
      LJACI = 1
      LJACJ = LJACI + NY + 1
      myid = isenfo(8)
      numprocs = isenfo(9)
      do ip = 1, isenfo(1)
         myip = myid+1 + (ip-1)*numprocs
         IPOS = myip*NY+1
         CALL AMUV2SP(NY,Y(IPOS),YPRIME(IPOS),SENWRK(LJAC),
     *        ISENWK(LJACJ), ISENWK(LJACI), DELTA(IPOS))
         if (myip .le. isenfo(4)) then
            ipos = ipos - 1
            DO I = 1, NY 
               II = IPOS + I
               DELTA(II) = DELTA(II) + SENWRK(II-NY)  
            END DO
         end if
      END DO
      RETURN
      END         
C
      SUBROUTINE RESIDX2(
     1     T,Y,YPRIME,CJ,DELTA,IRESI,RPAR,IPAR,SENPAR,
     2     WT,NY,IDE,RES,G_RES,CNST,RWK,IWK,ISENFO,A_RES,ADI)
C================================================================
C    This routine caculates 
C          g_uu' + g_t =0 
C    for algebraic constraints in a DAEs
C          0   = g(u, t)
C     If the sensitivity is also considered, the equation is
C          0   = g_uu_p + g_p = h(u, u_p, t).
C     This routine will calculate
c          0   = h_u u' + h_(u_p) u_p' + h_t
C    
C================================================================
      IMPLICIT NONE
      REAL*8 T, Y(*), YPRIME(*), CJ, DELTA(*), RPAR(*),
     *     SENPAR(*), WT(*), CNST, RWK(*), ADI(*)
      INTEGER IRESI, IPAR(*), NY, IDE(*),  IWK(*), ISENFO(*)
      EXTERNAL RES, G_RES, A_RES
C===  Local variables
      INTEGER I, NP, II, IRES, J
      INTEGER NEQ, NRE, NSE
      REAL*8  VNORM, DEL, DELINV, G_T
C
      IF (ISENFO(6) .LT. 0) THEN ! Adjoint equation evaluations
         CALL RESADIDX2(
     *        T,Y,YPRIME,CJ,DELTA,IRESI,RPAR,IPAR,SENPAR,
     *        RES, NY, IDE, ISENFO, RWK, IWK, G_RES, ADI)
         RETURN
      END IF
      NP = ISENFO(4)
      IRES = 0
      IF (ISENFO(2) .EQ. 2) THEN
         IRES = 2
         CALL RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,ADI)
         RETURN
      ELSE IF (ISENFO(2) .LT. 2) THEN
C
C     result of finite difference, require 2*NY rwork space.
C     Determine the perturbation.
C
         IF (IRESI .EQ. 0) THEN
            VNORM = 0.0D0
            DO I = 1,NY
               VNORM = VNORM + (WT(I)*YPRIME(I))**2
            END DO
            VNORM = DSQRT(VNORM)
            DEL = CNST/MAX(VNORM, 1.0D0)
C     
C     Perturb Y
C     
            DO I = 1,NY
               RWK(I) =  Y(I) + DEL*YPRIME(I)
            ENDDO
C     
C     Call RES with perturbed values.
C     
            CALL RES(T+DEL,RWK,YPRIME,CJ,RWK(NY+1),IRES,RPAR,
     *           IPAR,SENPAR,ADI)
            IF(IRES.LT.0)RETURN
C     
C     Approx. first order sensitivity, 
C     
            DELINV = 1.0D0/DEL
            DO I = 1,NY
               IF (IDE (I) .EQ. 1) THEN
                  DELTA(I)=DELINV*(RWK(NY+I)-DELTA(I))
               END IF
            END DO
         ELSE IF (IRESI .EQ. 1) THEN
C
C     If sensitivity is considered,
            VNORM = 0.0D0
            NEQ = NY*(ISENFO(1)+1)
            DO I = 1,NEQ
               VNORM = VNORM + (WT(I)*YPRIME(I))**2
            END DO
            VNORM = DSQRT(VNORM)
            DEL = CNST/MAX(VNORM, 1.0D0)
C     
C     Perturb Y
C     
            DO I = 1,NEQ
               RWK(I) =  Y(I) + DEL*YPRIME(I)
            ENDDO
C     
C     Call RES with perturbed values.
C     
            CALL DDSEN(T+DEL,RWK,YPRIME,CJ,RWK(NEQ+1),IRESI,RPAR,
     *           IPAR,RES,NY,ISENFO,
     *           RWK(2*NEQ+1),IWK,NRE,NSE,WT,CNST,SENPAR,G_RES,
     *           A_RES, ADI)
            IF(IRESI.LT.0)RETURN
C     
C     Approx. first order sensitivity, 
C     
            DELINV = 1.0D0/DEL
            DO I = 1,NY
               IF (IDE (I) .EQ. 1) THEN
                  DO J = 1, ISENFO(1)
                     II = J*NY + I
                     DELTA(II)=DELINV*(RWK(NEQ+II)-DELTA(II))
                  END DO
               END IF
            END DO
         END IF
      ELSE IF (ISENFO(1) .EQ. 0) THEN
C
C     If no sensitivity is considered
         G_T = 1.0D0
         CALL G_RES(
     *        T, G_T, Y, YPRIME, YPRIME, CJ, DELTA,
     *        RWK, IRES, RPAR, IPAR, SENPAR)
         DO I = 1, NY
            IF (IDE (I) .EQ. 1) THEN
               DELTA(I) = RWK(I)
            END IF
         END DO
      ELSE 
C
C     If sensitivity is considered
c
         call ri2AdfMV(
     1        t, y, yprime, cj, delta, iresi, 
     1        rpar, ipar, senpar, ny, np, 
     1        ide, g_res, rwk, isenfo, ADI)
      END IF 
C
C     Rescale the residuals to the original form
C
      DO J = 1, ISENFO(1)+1
         II = (J-1)*NY
         DO I = 1, NY
            IF (IDE(I) .EQ. 1) DELTA(II+I) = DELTA(II+I)/CJ
         END DO
      END DO
      RETURN
C--------------End of RESIDX2----------------------------------------
      END
C
      SUBROUTINE ri2AdfMV(
     1     t, y, yp, cj, delta, ires, 
     1     rpar, ipar, senpar, ny, np, 
     1     ide, g_res, wrk, isenfo, ADI)
C==============================================================
C
C    This routine calculates the residual of 
C         g_u u' + g_t = 0
C    for algebraic constraints in the sensitivity analysis
C         g(u) = 0
C    via ADIFOR with matrix-vector product only method
C
C==============================================================
C
      implicit none
      real*8   t                ! time
      real*8   y(*)             ! solutions
      real*8   yp(*)            ! derivatives
      real*8   delta(*)         ! residuals
      real*8   rpar(*), senpar(*) ! parameters
      real*8   cj               ! coefficient of Jacobian
      integer  ide(*)           ! indicator on which is index-2 constraint
      integer  ipar(*)          ! integer parameters
      integer  ny               ! number of variables 
      integer  np               ! number of sensitivities
      external g_res            ! Adifor-generated routine
      integer  ires             ! other parameters
      integer  isenfo(*)
      real*8   wrk(*), adi(*)   ! work array at least 2*ny + np
c === local variables
      integer  i, j, ig_delta, itotal, ig_senpar,
     *         igg_delta, iresloc
      real*8   g_t
      data     g_t/1.0d0/
C           
      ig_senpar = 1
      ig_delta  = ig_senpar + np
      igg_delta = ig_delta + ny
      itotal = igg_delta + ny
C      
      iresloc = 0
C     
      do i = 1, isenfo(1)
         if (i .le. isenfo(4)) wrk(i) = 1.0d0
         call g_res(t,g_t,y,yp,y(1+i*Ny),yp(1+i*Ny),yp,yp(1+i*Ny),
     1        cj,delta,wrk(ig_delta),delta(1+i*Ny),wrk(igg_delta),
     2        iresloc,rpar,ipar,senpar,wrk)
         if (i .le. isenfo(4)) wrk(i) = 0.0d0
         if (ires .eq. 0) then
            do j = 1, ny
               if (ide(j) .eq. 1) then
                  delta(j) = wrk(j+ig_delta-1)
               end if
            end do
            return
         end if
         do j = 1, ny
            if (ide(j) .eq. 1) delta(i*ny + j) = wrk(j+igg_delta-1)
         end do
      end do
      return
      end
C      
      SUBROUTINE jdAdfsp(
     1     T,y, yp, PD,CJ,RPAR,IPAR, SENPAR, NY, idense, MU, ML,
     2     i_res,ires, wrk, indvec, iwk, nq, iad, index2ad, iwkad)
c=======================================================================
c
C    Jacobian computation via sparse-linC Adifor. 
C
C           A step by step procedure to generate I_RES via ADIFOR with
C           SparsLinC option is 
C            1. Put all the codes related to
C               subroutine RES(
C              *      T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C               in a file called "res.f"
C            2. Create a file "res.cmp" with one line:
C                res.f
C            3. Create a file "resjacsp.adf" with the following lines:
C               AD_PROG = res.cmp
C               AD_TOP = res
C               AD_IVARS = y, yprime
C               AD_OVARS = delta
C               AD_PREFIX = i
C               AD_OUTPUT_DIR = jacsp
C               AD_FLAVOR = sparse
C            4. Run Adifor to generate I_RES(...) with command:
C               % Adifor AD_SCRIPT=resjacsp.adf
C
C            I_RES is stored in jacsp/i_res.f. 
C            The generated routine I_RES may have different argument list
C            depending on whether or not RPAR(*) have data dependence
C            with Y(*) and YPRIME(*).   
C             
C**** Copyright (C) 1998, Shengtai Li
C         
C=======================================================================   
      IMPLICIT NONE
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  pd(*)             ! Jacobian matrix
      real*8  cj                ! const coming from the DDASPK
      integer ires              ! error indicator 
      integer idense            ! indicator on whether the Jacobian is dense
                                !    1----banded
                                !    0----dense
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! parameter array for sensitivity
      integer ipar(*)           ! integer parameter of the problems
      external i_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables only
      integer MU                ! upper half bandwidth
      integer ML                ! lower half bandwidth
      real*8  wrk(*)            ! work array, size = ny
      integer indvec(*)         ! temp iwork array
      integer iwk(*)            ! integer work array, 
                                ! size = 3*ny
      integer nq                ! number of quadratures
      integer iad               ! for adjoint equation or not
      integer index2ad
      integer iwkad(*)          ! iwork space for adjoint method             
*      save iwk
c === local variable
      integer i, j, ig_y, ig_yp, ig_delta, width, 
     *     mband, outlen, info, ii,  itotal, nymnq
      character msg*80
      logical first
      save first
      data first/.true./
c
      ig_y     = 0
      ig_yp    = ig_y + ny
      ig_delta = ig_yp + ny
      itotal  = ig_delta + ny
      nymnq = ny - nq
      if (idense .eq. 1) then
         width = ML + mu + 1
         mband = 2*ML + mu
         if (iad .lt. 0) mband = 2*MU + ML
      else if (idense .eq. 0) then
         width = nymnq
      end if
      if (first) then
c
c     === initialize the sparselink       
         call xspini
         first = .false.
      end if
      do i = 1, nymnq
         call dspsd(iwk(ig_y  + i), i, 1.0d0, 1)
         call dspsd(iwk(ig_yp + i), i, cj, 1)
      end do
      do i = nymnq+1, ny
         call dspsd(iwk(ig_y  + i), i, 0.0d0, 1)
         call dspsd(iwk(ig_yp + i), i, 0.0d0, 1)
      end do
c
c === call Adifor generated routine
c
      call i_res(
     $     t, y, iwk(ig_y+1), yp,  iwk(ig_yp+1), cj, wrk, 
     $     iwk(ig_delta+1), ires, rpar, ipar, senpar)
c
c      call ehrpt
c
      if (idense .eq. 0) then   ! dense matrix
         do i = 1, nymnq
            call dspxsq(indvec,wrk,width,iwk(ig_delta+i),
     $           outlen,info)
            if (iad .ge. 0) then
               do j = 1, outlen
                  PD(i+(indvec(j)-1)*Nymnq)=wrk(j)
               end do              
            else if (index2ad .eq. 1) then
               do ii = 1, outlen
                  j = indvec(ii)
                  if (iwkad(j).eq.-2) then
                     PD(j+(i-1)*Nymnq)=wrk(ii)*cj
                  else if(iwkad(i+ny) .eq. -2) then 
                     PD(j+(i-1)*Nymnq)=wrk(ii)/cj
                  else
                     PD(j+(i-1)*Nymnq)=wrk(ii)
                  end if
               end do
            else
               do ii = 1, outlen
                  j = indvec(ii)
                  PD(j+(i-1)*Nymnq)=wrk(ii)
               end do
            end if
         end do
      else if (idense.eq. 1) then ! banded matrix
         do i = 1, nymnq
            call dspxsq(indvec,wrk,width,iwk(ig_delta+i),
     $           outlen,info)
            if (info .lt. 0) then
               msg = 'DASPK-- ERROR IN JdADFSP'
               CALL XERRWD(MSG,24,40,0,0,0,0,1,0.D0,0.0D0)
               ires = info
               RETURN
            end if
            if (iad .ge. 0) then
               do ii = 1, outlen
                  j = indvec(ii)
                  pd(i+j*mband - ML) = wrk(ii)
               end do
            else if (index2ad .eq. 1) then
               do ii = 1, outlen
                  j = indvec(ii)
                  if (iwkad(j).eq.-2) then                  
                     pd(j+i*mband - MU) = wrk(ii)*cj
                  else if(iwkad(i+ny) .eq. -2) then
                     pd(j+i*mband - MU) = wrk(ii)/cj
                  else
                     pd(j+i*mband - MU) = wrk(ii)
                  end if                    
               end do
            else
               do ii = 1, outlen
                  j = indvec(ii)
                  pd(j+i*mband - MU) = wrk(ii)
               end do          
            end if
         end do
      end if
      return
C----------------End of jdAdfsp--------------------------------------
      end
c
      SUBROUTINE jidAdfsp(
     1     T,y, yp, PD,CJ,RPAR,IPAR, SENPAR, NY, idense, MU, ML,
     2     id, icopt, i_res,ires, wrk, indvec, iwk, nq, iad)
c=======================================================================
c
C    Jacobian computation for initial conditions via sparsLinC Adifor.
C
C           A step by step procedure to generate I_RES via ADIFOR with
C           SparsLinC option is 
C            1. Put all the codes related to
C               subroutine RES(
C              *      T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C               in a file called "res.f"
C            2. Create a file "res.cmp" with one line:
C                res.f
C            3. Create a file "resjacsp.adf" with the following lines:
C               AD_PROG = res.cmp
C               AD_TOP = res
C               AD_IVARS = y, yprime
C               AD_OVARS = delta
C               AD_PREFIX = i
C               AD_OUTPUT_DIR = jacsp
C               AD_FLAVOR = sparse
C            4. Run Adifor to generate I_RES(...) with command:
C               % Adifor AD_SCRIPT=resjacsp.adf
C
C            I_RES is stored in jacsp/i_res.f. 
C            The generated routine I_RES may have different argument list
C            depending on whether or not RPAR(*) have data dependence
C            with Y(*) and YPRIME(*).   
C             
C**** Copyright (C) 1998, Shengtai Li
C         
C=======================================================================   
      IMPLICIT NONE
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  pd(*)             ! Jacobian matrix
      real*8  cj                ! const coming from the DDASPK
      integer ires              ! error indicator 
      integer idense            ! indicator on whether the Jacobian is dense
                                !    1----banded
                                !    0----dense
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! parameter array for sensitivity
      integer ipar(*)           ! integer parameter of the problems
      integer id(*)             ! integer array for indicator of variables
      integer icopt             ! option for initialization
      external i_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables only
      integer MU                ! upper half bandwidth
      integer ML                ! lower half bandwidth
      real*8  wrk(*)            ! work array, size = ny
      integer indvec(*)
      integer iwk(*)            ! integer work array, 
                                ! size = 3*ny
      integer nq                ! number of quadratures
      integer iad               ! monitor for adjoint method
c      save iwk
c === local variable
      integer i, j, ig_y, ig_yp, ig_delta, width, 
     *     mband, outlen, info, ii,  itotal, nymnq
      character msg*80
      logical first
      save first
      data first/.true./
c
      ig_y     = 0
      ig_yp    = ig_y + ny
      ig_delta = ig_yp + ny
      itotal  = ig_delta + ny
      nymnq = ny - nq
      if (idense .eq. 1) then
         width = ML + mu + 1
         mband = 2*ML + mu
         if (iad .lt. 0) mband = 2*MU + ML 
      else if (idense .eq. 0) then
         width = nymnq
      end if
      if (first) then      
         call xspini
         first = .false.
      end if
c
c...  initialize the seed matrix first
      do i = 1, ny
         call dspsd(iwk(ig_y  + i), i, 0.0d0, 1)
         call dspsd(iwk(ig_yp + i), i, 0.0d0, 1)
      end do
         
      do i = 1, nymnq
         if (id(i) .lt. 0) then
c
c...  algebraic variables
            call dspsd(iwk(ig_y + i), i, 1.0d0, 1)
         end if
      end do
      if (icopt .EQ. 3 .OR. icopt .EQ. 4) then
c.... icopt = 3,4
         do i = 1, nymnq
            if (id(i) .eq. 1 .or. id(i) .eq. 2) then
c
c...  differential variables, fixed Y, free YPRIME
               call dspsd(iwk(ig_yp + i), i, cj, 1)
            end if
         end do
         do i = 1, nymnq
            if (id(i) .eq.1 .or. id(i).eq. 3) then
c
c...  differential variables, fixed YPRIME, free Y
               call dspsd(iwk(ig_y + i), i, 1.0d0, 1)
            end if
         end do
      else if (icopt .eq. 1) then
c.....icopt = 1
         do i = 1, nymnq
            if (id(i) .gt. 0) then
c
c...  differential variables, fixed Y, free YPRIME
               call dspsd(iwk(ig_yp + i), i, cj, 1)
            end if
         end do
      else if (icopt .EQ. 5) then
c.....icopt = 5
         do i = 1, nymnq
            call dspsd(iwk(ig_y + i), i, 1.0d0, 1)
         end do         
      end if 
c
c === call Adifor generated routine
c
      call i_res(
     $     t, y, iwk(ig_y+1), yp,  iwk(ig_yp+1), cj, wrk, 
     $     iwk(ig_delta+1), ires, rpar, ipar, senpar)
      if (ires .lt. 0) return
c
c       call ehrpt
c
      if (icopt .lt. 5) then
         if (idense .eq. 0) then ! dense matrix
            do i = 1, nymnq
               call dspxsq(indvec,wrk,width,iwk(ig_delta+i),
     $              outlen,info)
               if (info .lt. 0) then
                  msg = 'DASPK-- ERROR IN JidADFSP'
                  CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                  ires = info
                  RETURN
               end if
               if (iad .lt. 0) then
                  do ii = 1, outlen
                     j = indvec(ii)
                     if (id(j)*id(i+ny) .gt. 0) then
                        PD(j+(i-1)*Nymnq)=wrk(ii)
                     end if
                  end do
               else
                  do j = 1, outlen
                     PD(i+(indvec(j)-1)*Nymnq)=wrk(j)
                  end do
               end if
            end do
         else if (idense.eq. 1) then ! banded matrix
            do i = 1, nymnq
               call dspxsq(indvec,wrk,width,iwk(ig_delta+i),
     $              outlen,info)  
               if (info .lt. 0) then
                  msg = 'DASPK-- ERROR IN JidADFSP'
                  CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                  ires = info
                  RETURN
               end if
               if (iad .lt. 0) then
                  do ii = 1, outlen
                     j = indvec(ii)
                     if (id(j)*id(i+ny) .gt. 0) then
                        pd(j+i*mband - MU) = wrk(ii)
                     end if
                  end do
               else
                  do ii = 1, outlen
                     j = indvec(ii)
                     pd(i+j*mband - ML) = wrk(ii)
                  end do
               end if
            end do
         end if
         if (iad .ge. 0) return ! forward initialization
c
c...  second step for the adjoint Jacobian computation
         do i = 1, ny
            call dspsd(iwk(ig_y  + i), i, 1.0d0, 1)
            call dspsd(iwk(ig_yp + i), i, 0.0d0, 1)
         end do
         call i_res(
     $        t, y, iwk(ig_y+1), yp,  iwk(ig_yp+1), cj, wrk, 
     $        iwk(ig_delta+1), ires, rpar, ipar, senpar)
         if (ires .lt. 0) return
         if (idense .eq. 0) then ! dense matrix
            do i = 1, nymnq
               if (id(i+ny) .lt. 0) then
                  call dspxsq(indvec,wrk,width,iwk(ig_delta+i),
     $                 outlen,info)  
                  do ii = 1, outlen
                     j = indvec(ii)
                     if (id(j) .gt. 0) then
                        PD(j+(i-1)*Nymnq)=wrk(ii)
                     end if
                  end do
               end if
            end do
         else if (idense.eq. 1) then ! banded matrix
            do i = 1, nymnq
               if (id(i+ny) .lt. 0) then
                  call dspxsq(indvec,wrk,width,iwk(ig_delta+i),
     $                 outlen,info)  
                  if (info .lt. 0) then
                     msg = 'DASPK-- ERROR IN JidADFSP'
                     CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                     ires = info
                     RETURN
                  end if
                  do ii = 1, outlen
                     j = indvec(ii)
                     if (id(j) .gt. 0) then
                        pd(j+i*mband - MU) = wrk(ii)
                     end if
                  end do
               end if
            end do
         end if
         return
      else                      ! icopt == 5
c
c     Second stage for index-2 DAEs
C
         if (idense .eq. 0) then ! dense matrix
            if (iad .lt. 0) then
               do i = 1, nymnq
                  if (id(ny+i) .lt. 0) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)
                     if (id(ny+i) .eq. -2) then
                        do j = 1, outlen
                           PD(indvec(j) + (i-1)*Nymnq)=wrk(j)/cj
                        end do
                     else
                        do j = 1, outlen
                           PD(indvec(j) + (i-1)*Nymnq)=wrk(j)
                        end do
                     end if
                  end if
               end do
            else
               do i = 1, nymnq
                  if (id(ny+i) .eq. 1) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)  
                     do j = 1, outlen
                        PD(i+(indvec(j)-1)*Nymnq)=wrk(j)
                    end do
                  end if
               end do
            end if
         else if (idense.eq. 1) then ! banded matrix
            if (iad .lt. 0) then
               do i = 1, nymnq
                  if (id(ny+i) .lt. 0) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)  
                     if (info .lt. 0) then
                        msg = 'DASPK-- ERROR IN JidADFSP'
                        CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                        ires = info
                        RETURN
                     end if
                     if (id(ny+i) .eq. -2) then
                        do ii = 1, outlen
                           j = indvec(ii)
                           pd(j+i*mband - MU) = wrk(ii)/cj
                        end do
                     else
                        do ii = 1, outlen
                           j = indvec(ii)
                           pd(j+i*mband - MU) = wrk(ii)
                        end do
                     end if  
                  end if
               end do
            else
               do i = 1, nymnq
                  if (id(ny+i) .eq. 1) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)  
                     if (info .lt. 0) then
                        msg = 'DASPK-- ERROR IN JidADFSP'
                        CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                        ires = info
                        RETURN
                     end if
                     do ii = 1, outlen
                        j = indvec(ii)
                        pd(i+j*mband - ML) = wrk(ii)
                     end do
                  end if
               end do 
            end if
         end if
c
         do i = 1, ny
            if (id(i) .gt. -2) then
c
c...  fixed Y, compute YP
               call dspsd(iwk(ig_y + i),  i, 0.0d0, 1)
               call dspsd(iwk(ig_yp + i), i, cj, 1)
            end if
         end do         
c
c === call Adifor generated routine again
c
         call i_res(
     $        t, y, iwk(ig_y+1), yp,  iwk(ig_yp+1), cj, wrk, 
     $        iwk(ig_delta+1), ires, rpar, ipar, senpar)
c
c     compute the part corresponding to the index-0 and 1 equations
c
         if (idense .eq. 0) then ! dense matrix
            if (iad .lt. 0) then
               do i = 1, nymnq
                  if (id(ny+i) .gt. 0) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)
                     do ii = 1, outlen
                        j = indvec(ii)
                        if (id(j) .gt. 0) then
                           PD(j+(i-1)*Nymnq)=wrk(ii)
                        else if (id(j) .eq. -2) then
                           PD(j+(i-1)*Nymnq)=wrk(ii)*cj
                        end if
                     end do
                  end if
               end do
            else
               do i = 1, nymnq
                  if (id(ny+i) .eq. 0) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)  
                     do j = 1, outlen
                        PD(i+(indvec(j)-1)*Nymnq)=wrk(j)
                     end do
                  end if
               end do
            end if
         else if (idense.eq. 1) then ! banded matrix
            if (iad .lt. 0) then
               do i = 1, nymnq
                  if (id(ny+i) .gt. 0) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)
                     if (info .lt. 0) then
                        msg = 'DASPK-- ERROR IN JidADFSP'
                        CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                        ires = info
                        RETURN
                     end if
                     do ii = 1, outlen
                        j = indvec(ii)
                        if (id(j) .gt. 0) then
                           pd(j+i*mband - MU) = wrk(ii)
                        else if (id(j) .eq. -2) then
                           pd(j+i*mband - MU) = wrk(ii)*cj
                        end if
                     end do
                  end if
               end do
            else
               do i = 1, nymnq
                  if (id(ny+i) .eq. 0) then
                     call dspxsq(
     $                    indvec,wrk,width,iwk(ig_delta+i),
     $                    outlen,info)  
                     if (info .lt. 0) then
                        msg = 'DASPK-- ERROR IN JidADFSP'
                        CALL XERRWD(MSG,25,40,0,0,0,0,1,0.D0,0.0D0)
                        ires = info
                        RETURN
                     end if
                     do ii = 1, outlen
                        j = indvec(ii)
                        pd(i+j*mband - ML) = wrk(ii)
                     end do
                  end if
               end do
            end if
         end if
      end if
      return
      end
c
      SUBROUTINE jidAdfSM(
     1     T, y, yp, PD,CJ,RPAR,IPAR, senpar, NY, idense, MU, ML,  
     2     id, icopt, j_res, ires, delta, wrk, nq, iad)
c=======================================================================
c
C  Jacobian for initialization by ADIFOR with seed matrix option
C
C  This version is for using the work array
C
c   How to use this routine:
C    1. Put all the routines related to 
C       subroutine RES(
C      *         T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C       in a file called "res.f"
C    2. Create file "res.cmp" with one line:
C        res.f
C    3. Create file "resjac.adf":
C        AD_PROG = res.cmp
C        AD_TOP = res
C        AD_IVARS = y, yprime
C        AD_OVARS = delta
C        AD_PREFIX = j
C        AD_SUPPRESS_LDG = true
C        AD_PMAX = bandwidth = (ML + MU + 1)
C        AD_OUTPUT_DIR = jac
C
C    4. run Adifor generate j_res(...) with
C       % Adifor AD_SCRIPT=resjac.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1998, Shengtai Li
C         
C=======================================================================   
      IMPLICIT NONE
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  pd(*)             ! Jacobian matrix
      real*8  cj                ! const coming from the DDASPK
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
      external j_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables only
      integer idense            ! indicator on whether the Jacobian is dense
                                !    1----banded
                                !    0----dense
      integer MU                ! upper half bandwidth
      integer ML                ! lower half bandwidth
      integer id(*)             ! integer array for indicator of variables
      integer icopt             ! initialization options
      real*8  delta(*)          ! work array for residuals, size = NY
      real*8  wrk(*)            ! work array, size = width*(3*NY)
      integer nq                ! number of quadratures
      integer iad
c === local variable
      integer i, j, ig_y, ig_yp, ig_delta, ipos, width, ll,ul, 
     *     length, mebm1, ii,jj, itotal, nymnq
c
c === define the width of the seed matrix
c
      nymnq = ny - nq
      if (idense .eq. 1) then
         WIDTH=min(nymnq,1+ML+MU)
      else 
         width = nymnq
      end if
c
c === set index for the work array
c
      length = ny*width
      ig_y = 0
      ig_yp = ig_y + length
      ig_delta = ig_yp + length
      itotal = ig_delta + length
c
c === initialize the seed matrix
c
      do i = 1, itotal
         wrk(i) = 0.0d0
      end do
      mebm1 = 2*ML + MU
      if (iad .lt. 0) mebm1 = 2*MU + ML
      DO I=1,Nymnq
         if (idense .eq. 1) then
            ipos = (i-1)*width + mod(i-1,width) + 1
         else
            ipos = I + (I-1)*width
         end if
         If (icopt .eq. 1) then
            if (id(i) .lt. 0) then
               wrk(ig_y + ipos ) = 1.0d0
            else
               wrk(ig_yp+ ipos ) = cj
            end if
         else if (icopt .EQ. 3 .OR. icopt .EQ. 4) then
            if (id(i) .lt. 2) then
               wrk(ig_y + ipos )=1.0D0
               wrk(ig_yp+ ipos )=cj
            else if (id(i) .eq. 2) then
               wrk(ig_yp+ ipos )=cj
            else if (id(i) .eq. 3) then
               wrk(ig_y + ipos )=1.0D0
            end if
         else if (icopt .EQ. 5) then
            wrk(ig_y + ipos ) = 1.0d0
         end if               
      end do
c
c === call Adifor generated routine
c
      call j_res(
     *     width, t, y, wrk(ig_y+1), yp, wrk(ig_yp+1), 
     *     cj, delta,  wrk(ig_delta+1), 
     *     ires, rpar, ipar, senpar)
c
      if (icopt .lt. 5) then
         if (idense .eq. 1) then
            do i = 1, nymnq
               LL=max(1,I-MU)
               UL=min(Nymnq,I+ML)
               jj = mod(i-1, width) + 1
               if (iad .lt. 0) then
                  DO j=LL,UL
                     if (id(j+ny)*id(i) .gt. 0) then
                        pd(j*mebm1-MU+i)=wrk(ig_delta + (j-1)*width +jj)
                     end if
                  end do
               else
                  II=i*MEBM1-ML
                  DO j=LL,UL
                     pd(II+j)=wrk(ig_delta + (j-1)*width + jj)
                  end do
               end if
            end do
         else if (idense .eq. 0) then ! dense case
            do i=1, Nymnq
               ipos = (I-1)*width
               if (iad .lt. 0) then
                  do j=1,Nymnq
                     if (id(i)*id(j+ny) .gt. 0) then
                        PD(j + ipos)=wrk(ig_delta+ j + ipos)
                     end if
                  end do
               else
                  do j=1,Nymnq
                     PD(I + (j-1)*width)=wrk(ig_delta+ j + ipos)
                  end do
               end if
            end do
         end if
         if (iad .ge. 0) return ! forward initialization
c
c...  second step for the adjoint Jacobian computation
         do i = 1, nymnq
            if (idense .eq. 1) then
               ipos = (i-1)*width + mod(i-1,width) + 1
            else
               ipos = I + (I-1)*width
            end if
            wrk(ig_y + ipos )=1.0D0
            wrk(ig_yp + ipos )=0.0D0
         end do         
c
c === call Adifor generated routine again
c
         call j_res(
     *        width, t, y, wrk(ig_y+1), yp, wrk(ig_yp+1), 
     *        cj, delta,  wrk(ig_delta+1), 
     *        ires, rpar, ipar, senpar)
         if (ires .lt. 0) return
         if (idense .eq. 1) then
            do i = 1, nymnq
               if (id(i) .gt. 0) then
                  LL=max(1,I-MU)
                  UL=min(Nymnq,I+ML)
                  jj = mod(i-1, width) + 1
                  DO j=LL,UL
                     if (id(j+ny) .lt. 0) then
                        pd(j*mebm1-MU+i)=wrk(ig_delta + (j-1)*width +jj)
                     end if
                  end do
               end if
            end do
         else if (idense .eq. 0) then ! dense case
            do i=1, Nymnq
               if (id(i) .gt. 0) then
                  ipos = (I-1)*width
                  do j=1,Nymnq
                     if (id(j+ny) .lt. 0) then
                        PD(j + ipos)=wrk(ig_delta+ j + ipos)
                     end if
                  end do
               end if
            end do
         end if
         return
      else                      ! icopt == 5
c
C     two-stage initialization for index-2 and general index-1
C
         if (idense .eq. 1) then ! banded case
            if (iad .lt. 0) then
               do i = 1, nymnq                  
                  LL=max(1,I-MU)
                  UL=min(Nymnq,I+ML)
                  jj = mod(i-1, width) + 1
                  ii = i - MU
                  DO j=LL,UL
                     if (id(j+ny) .eq. -2) then
                        pd(j*mebm1+ii)=wrk(ig_delta+(j-1)*width + jj)/cj
                     else if (id(j+ny) .lt. 0) then
                        pd(j*mebm1+ii)=wrk(ig_delta+(j-1)*width + jj)                       
                     end if
                  end do
               end do                    
            else
               do i = 1, nymnq
                  LL=max(1,I-MU)
                  UL=min(Nymnq,I+ML)
                  jj = mod(i-1, width) + 1
                  II=i*MEBM1-ML
                  DO j=LL,UL
                     if (id(j+ny) .eq. 1) then
                        pd(II+j)=wrk(ig_delta + (j-1)*width + jj)
                     end if
                  end do
               end do
            end if
         else if (idense .eq. 0) then
            if (iad .lt. 0) then
               do i=1, Nymnq
                  if (id(i+ny) .eq. -2) then
                     ipos = (I-1)*Nymnq
                     do j=1,Nymnq
                        PD(j + ipos)=wrk(ig_delta+ j+ipos)/cj
                     end do
                  else if (id(i+ny) .lt. 0) then
                     ipos = (I-1)*Nymnq
                     do j=1,Nymnq
                        PD(j + ipos)=wrk(ig_delta+ j+ipos)
                     end do
                  end if
               end do
            else
               do i=1, Nymnq
                  if (id(i+ny) .eq. 1) then
                     ipos = (I-1)*Nymnq
                     do j=1,Nymnq
                        PD(I + (J-1)*Nymnq )=wrk(ig_delta+ j + ipos)
                     end do
                  end if
               end do
            end if
         end if
c
         do i = 1, nymnq
            if (idense .eq. 1) then
               ipos = (i-1)*width + mod(i-1,width) + 1
            else
               ipos = I + (I-1)*width
            end if
            if (id(i) .gt. -2) then
               wrk(ig_y + ipos )=0.0D0
               wrk(ig_yp + ipos )=cj
            end if
         end do         
c
c === call Adifor generated routine again
c
         call j_res(
     *        width, t, y, wrk(ig_y+1), yp, wrk(ig_yp+1), 
     *        cj, delta,  wrk(ig_delta+1), 
     *        ires, rpar, ipar, senpar)
         if (ires .lt. 0) return
c
         if (idense .eq. 0) then ! dense matrix
            if (iad .lt. 0) then
               do i = 1, Nymnq
                  if (id(i+ny) .gt. 0) then
c
c...  differential equations
                     ipos = (I-1)*Nymnq
                     do j = 1, Nymnq
                        if (id(j) .gt. -2) then
                           PD(j + ipos)=wrk(ig_delta+ j + ipos)
                        else 
                           PD(j + ipos)=wrk(ig_delta+ j+ipos)*cj
                        end if
                     end do
                  end if
               end do
            else
               do i=1, Nymnq
                  if (id(i+ny) .eq. 0) then
                     ipos = (I-1)*Ny
                     do j=1,Nymnq
                        PD(I + (J-1)*Nymnq)=wrk(ig_delta+ j + ipos)
                     end do
                  end if
               end do
            end if
         else if (idense .eq. 1) then ! BANDED CASE
            if (iad .lt. 0) then
               do i = 1, nymnq
                  LL=max(1,I-MU)
                  UL=min(Nymnq,I+ML)
                  jj = mod(i-1, width) + 1
                  ii = i - MU
                  DO j=LL,UL
                     if (id(j+ny) .gt. 0) then
                        if (id(i) .gt. -2) then
                           pd(j*mebm1+ii)=wrk(
     *                          ig_delta + (j-1)*width + jj)
                        else
                           pd(j*mebm1+ii)=wrk(
     *                          ig_delta + (j-1)*width + jj)*cj
                        end if
                     end if
                  end do
               end do
            else
               do i = 1, nymnq
                  LL=max(1,I-MU)
                  UL=min(Nymnq,I+ML)
                  jj = mod(i-1, width) + 1
                  II=i*MEBM1-ML
                  DO j=LL,UL
                     if (id(j+ny) .eq. 0) then
                        pd(II+j)=wrk(ig_delta + (j-1)*width + jj)
                     end if
                  end do
               end do
            end if
         end if
      end if
c
      RETURN
c--------------------End of jidAdfsm-------------------------------
      END
C
      SUBROUTINE jdAdfSM(
     1     T, y, yp, PD,CJ,RPAR,IPAR, senpar, NY, idense, MU, ML,  
     2     j_res, ires, delta, wrk, nq, iad, index2ad, iwkad)
c=======================================================================
c
C  This version is for using the work array
C
c   How to use this routine:
C    1. Put all the routines related to 
C       subroutine RES(T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR)
C       in a file called "res.f"
C    2. Create file "res.cmp" with one line:
C        res.f
C    3. Create file "resjac.adf":
C        AD_PROG = res.cmp
C        AD_TOP = res
C        AD_IVARS = y, yprime
C        AD_OVARS = delta
C        AD_PREFIX = j
C        AD_SUPPRESS_LDG = true
C        AD_PMAX = bandwidth = (ML + MU + 1)
C        AD_OUTPUT_DIR = jac
C
C    4. run Adifor generate j_res(...) with
C       % Adifor AD_SCRIPT=resjac.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1998, Shengtai Li
C         
C=======================================================================   
      IMPLICIT NONE
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  pd(*)             ! Jacobian matrix
      real*8  cj                ! const coming from the DDASPK
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
      external j_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables only
      integer idense            ! indicator on whether the Jacobian is dense
                                !    1----banded
                                !    0----dense
      integer MU                ! upper half bandwidth
      integer ML                ! lower half bandwidth
      real*8  delta(*)          ! work array for residuals, size = NY
      real*8  wrk(*)            ! work array, size = width*(3*NY)
      integer nq
      integer iad
      integer index2ad          ! index-2 adjoint system
      integer iwkad(*)          ! integer work space for adjoint method
c === local variable
      integer i, j, ig_y, ig_yp, ig_delta, ipos, width, ll,ul, 
     *     length, mebm1, ii,jj, itotal, nymnq
c
c === define the width of the seed matrix
c
      nymnq = ny - nq
      if (idense .eq. 1) then
         WIDTH=min(nymnq,1+ML+MU)
      else 
         width = nymnq
      end if
c
c === set index for the work array
c
      length = ny*width
      ig_y = 0
      ig_yp = ig_y + length
      ig_delta = ig_yp + length
      itotal = ig_delta + length
c
c === initialize the seed matrix
c
      do i = 1, itotal
         wrk(i) = 0.0d0
      end do
      if (idense .eq. 1) then
         mebm1 = 2*ML + MU
         if (iad .lt. 0) mebm1 = 2*MU + ML 
         DO I=1,Nymnq
            ipos = (i-1)*width + mod(i-1,width) + 1
            wrk(ig_y + ipos )=1.0D0
            wrk(ig_yp+ ipos )=cj
         end do
      else if (idense .eq. 0) then
         DO I=1, Nymnq
            ipos = I + (I-1)*Nymnq
            wrk(ig_y + ipos )=1.0D0
            wrk(ig_yp + ipos )=CJ
         end do
      end if
c
c === call Adifor generated routine
c
      call j_res(
     *     width, t, y, wrk(ig_y+1), yp, wrk(ig_yp+1), 
     *     cj, delta,  wrk(ig_delta+1), 
     *     ires, rpar, ipar, senpar)
c
      if (idense .eq. 1) then
         if (iad .ge. 0) then
c
c...  forward options
            do i = 1, nymnq
               LL=max(1,I-MU)
               UL=min(Nymnq,I+ML)
               II=i*MEBM1-ML
               jj = mod(i-1, width) + 1
               DO j=LL,UL
                  pd(II+j)=wrk(ig_delta + (j-1)*width + jj)
               end do
            end do
         else if (index2ad .eq. 1) then
c
c...  index-2 adjoint options
            do i = 1, nymnq
               LL=max(1,I-MU)
               UL=min(Nymnq,I+ML)
               jj = mod(i-1, width) + 1
               ii = i - MU
               DO j=LL,UL                  
                  if (iwkad(i).eq.-2) then
                     pd(j*mebm1+ii)=wrk(ig_delta+(j-1)*width + jj)*cj
                  else if (iwkad(j+ny) .eq. -2) then
                     pd(j*mebm1+ii)=wrk(ig_delta+(j-1)*width + jj)/cj
                  else
                     pd(j*mebm1+ii)=wrk(ig_delta+(j-1)*width + jj)
                  end if
               end do
            end do
         else
c
c...  non-index-2 adjoint options
            do i = 1, nymnq
               LL=max(1,I-MU)
               UL=min(Nymnq,I+ML)
               jj = mod(i-1, width) + 1
               ii = i - MU
               DO j=LL,UL                  
                  pd(j*mebm1+ii)=wrk(ig_delta+(j-1)*width + jj)
               end do
            end do            
         end if
      else if (idense .eq. 0) then
c
c...  dense matrix
         if (iad .ge. 0) then
            do i=1, Nymnq
               ipos = (I-1)*Nymnq
               do j=1,Nymnq
                  PD(I + (J-1)*Nymnq )=wrk(ig_delta+ j + ipos)
               end do
            end do
         else if (index2ad .eq. 1) then
            do i=1, Nymnq
               ipos = (I-1)*Nymnq
               if (iwkad(i+ny) .eq. -2) then
                  do j = 1, Nymnq
                     PD(j + ipos)=wrk(ig_delta + j + ipos)/cj
                  end do
               else
                  do j = 1, Nymnq
                     if (iwkad(j).eq.-2) then
                        PD(j + ipos)=wrk(ig_delta + j + ipos)*cj
                     else
                        PD(j + ipos)=wrk(ig_delta+ j + ipos)
                     end if
                  end do               
               end if
            end do
         else
            do i=1, Nymnq
               ipos = (I-1)*Nymnq
               do j=1,Nymnq
                  PD(j + ipos)=wrk(ig_delta+ j + ipos)
               end do
            end do
         end if
      end if
c
      RETURN
      END
      SUBROUTINE kmvbyAdf(
     *     T,Y,YP,CJ,DELTA,V,AV,IRES,RPAR,IPAR,
     *     SENPAR,K_RES,NY,WRK)
c===================================================================
c  This routine computes the matrix-vector product in Krylov iteration.
C     AD_SCALAR_GRADIENTs = true
c
c   How to use this routine:
C    1. put all the staffs related with 
C           subroutine RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,SENPAR)
C       in a file called "res.f"
C    2. create file "res.cmp" with one line as
C        res.f
C    3. create file "resmv.adf" as
C     AD_PROG = res.cmp
C     AD_TOP = res
C     AD_IVARS = y, yprime
C     AD_OVARS = delta
C     AD_PMAX = 1
C     AD_PREFIX = g
C     AD_SCALAR_GRADIENTS = true
C     AD_OUTPUT_DIR = resmv
C
C    4. run Adifor generate k_res(...) with
C       % Adifor AD_SCRIPT=resmv.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 2000, Shengtai Li
C         
C
C===================================================================
      implicit none
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  cj                ! const coming from the DDASPK
      real*8  delta(*)          ! residual 
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
      real*8  v(*)              ! vector in matrix-vector product
      real*8  av(*)             ! result for matrix-vector product
      external k_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables
      real*8  wrk(*)            ! work array, size >= ny
c === local variable
      integer i

      do i = 1, ny
         wrk(i) = cj*v(i)
      end do
      ires = 0
      call k_res(t,y,v,yp,wrk,cj,delta,av,ires,rpar,ipar,senpar)
      return
      end
c
      SUBROUTINE rAdfMV(
     *     T,Y,YP,CJ,DELTA,IRES,RPAR,IPAR,senpar,G_RES,NY,isenfo,
     *     wrk)
c===================================================================
c  This routine computes the state variable and sensitivity residuals 
C  through the Adifor with Matrix-vector product form.
C     AD_SCALAR_GRADIENTs = true
c
c   How to use this routine:
C    1. put all the staffs related with 
C         subroutine RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,SENPAR)
C       in a file called "res.f"
C    2. create file "res.cmp" with one line as
C        res.f
C    3. create file "res4senmv.adf" as
C     AD_PROG = res.cmp
C     AD_TOP = res
C     AD_IVARS = y, yprime, senpar
C     AD_OVARS = delta
C     AD_PMAX = 1
C     AD_PREFIX = g
C     AD_SCALAR_GRADIENTS = true
C     AD_OUTPUT_DIR = resmv
C
C    4. run Adifor generate g_res(...) with
C       % Adifor AD_SCRIPT=resmv.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1998, Shengtai Li
C         
C
C===================================================================
      implicit none
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  cj                ! const coming from the DDASPK
      real*8  delta(*)          ! residual 
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
      external g_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables
      integer isenfo(*)         ! information for sensitivity size>9
                                ! see routine ddsen(*) in DDASPK for details 
      real*8  wrk(*)            ! work array, size >= isenfo(1)
c === local variable
      integer i, myid, numprocs, myi
      do i = 1, isenfo(4)
         wrk(i) = 0.0d0
      end do
      myid = isenfo(8)
      numprocs = isenfo(9)
      do i = 1, isenfo(1)
         myi = myid+1 + (i-1)*numprocs
         if (myi .le. isenfo(4)) wrk(myi) = 1.0d0
         call g_res(t,y,y(1+i*Ny),yp,yp(1+i*Ny),
     1        cj,delta,delta(1+i*Ny),ires,rpar,ipar,senpar,wrk)
         if (myi .le. isenfo(4)) wrk(myi) = 0.0d0
      end do
      return
      end
C
C      
      SUBROUTINE rAdfSM(
     *     T,Y,YP,CJ,DELTA,IRES,RPAR,IPAR,senpar,G_RES,NY,isenfo,
     *     wrk)
c===================================================================
c  This routine computes the state variable and sensitivity residuals 
C  through the Adifor with seed matrix input.
C     AD_SCALAR_GRADIENTS = false
c
c   How to use this routine:
C    1. put all the staffs related with 
C           subroutine RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR,SENPAR)
C       in a file called "res.f"
C    2. create file "res.cmp" with one line as
C        res.f
C    3. create file "res4sensm.adf" as
C     AD_PROG = res.cmp
C     AD_TOP = res
C     AD_IVARS = y, yprime, senpar
C     AD_OVARS = delta
C     AD_PREFIX = h
C     AD_SUPPRESS_LDG = true
C     AD_PMAX   = # of sensitivities
C     AD_OUTPUT_DIR = res4sen
C
C    4. run Adifor generate h_res(...) with
C       % Adifor AD_SCRIPT=res4sensm.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1998, Shengtai Li
C         
C
C===================================================================
      implicit none
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  cj                ! const coming from the DDASPK
      real*8  delta(*)          ! residual 
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)         ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
      external g_res            ! residual generated by the Adifor
      integer ny                ! number of equations totally
      integer isenfo(*)         ! information for sensitivity size>9
                                ! see routine ddsen(*) in DDASPK for details 
      real*8  wrk(*)            ! work array, 
                                !size = isenfo(1)*(3*ny+isenfo(4))
c === local variable
      integer i, np, ig_y, ig_yp,ig_delta, 
     *     ipos1, ip, ipos2, length, ig_senpar, itotal
      np = isenfo(1)
      length = np*ny
C === set Pointer for the work array
      ig_senpar = 0
      ig_y  = ig_senpar + isenfo(4)*np
      ig_yp = ig_y + length
      ig_delta = ig_yp + length
      itotal = ig_delta + length
c === Initialize the work array for the residual array
      do i = ig_delta, itotal
         wrk(i) = 0.0d0
      end do
C
C === Set up the seed matrices for the call to RES
C
      do ip = 1, np
         ipos1 = ip*Ny
         do i=1,Ny
            ipos2 = ip + (i-1)*np
            wrk(ig_y + ipos2 ) = y(i+ ipos1 )
            wrk(ig_yp + ipos2 ) = yp(i+ ipos1 )
         end do
         do i = 1, isenfo(4)
            wrk(ip + (i-1)*np) = 0.0d0
         end do
C
C === Initialize the the vector that picks out the correct dF/dp(i)
C
         if (ip .le. isenfo(4)) wrk(ip+(ip-1)*np)=1.0d0
      end do
C
C Now call the ADIFOR generated  code 
c
      call g_res(np, t, y, wrk(ig_y+1), yp, wrk(ig_yp+1),
     *     cj, delta, wrk(ig_delta+1), ires,
     *     rpar, ipar, senpar, wrk)
c      
      do i=1,ny
         ipos1 = ig_delta + (i-1)*np
         do ip=1,np
            delta(i+ip*Ny)=wrk(ipos1 + ip)
         end do
      end do
c
      return
      end
C
      SUBROUTINE jrAdfsp (
     1     NY, T, Y, YPRIME, DELTA, I_RES, CJ, JAC, JA, IA, 
     2     IPAR, RPAR, SENPAR, IRES, DFDP, NP, IWK)
c==================================================================
c     Compute the Jacobian dF/dY, dF/dY', and dF/dP for residual
c     function F(Y, Y', P) = 0,
c     used when ISENFO(2) = 5.
c
C     Estimates of Jacobian elements are computed by Adifor2.1.
C     The Jacobian is stored in compressed sparse row format in 
C     JAC(ny, 1:2*ny) and DFDP(ny, np)
c==================================================================
      IMPLICIT NONE
C ... Input arguments:
      INTEGER NY                ! number of equations for state variables
      REAL*8 T                  ! independent variable t
      REAL*8 Y(*)               ! most recent iterate of solution vector y
      REAL*8 YPRIME(*)          ! most recent iterate of solution vector y'
      EXTERNAL I_RES            ! function generated by the Adifor
      INTEGER IRES              ! error flag for RES routine
      REAL*8 DFDP(NY,*)         ! DF/DP
      INTEGER NP                ! number of the sensitivities
      REAL*8 CJ                 ! scalar proportional to 1/H
      REAL*8 RPAR(*)            ! user real workspace
      real*8 SENPAR(*)          ! sensitivity parameters that appear in RES
      INTEGER IPAR(*)           ! user integer workspace
C ... Work-array argument:
      REAL*8 DELTA(*)           ! residual for the equations.
C ... Output arguments:
      REAL*8 JAC(*)             ! nonzero Jacobian elements
      INTEGER JA(*)             ! col indices of nonzero Jacobian elements
      INTEGER IA(NY+1)          ! pointers to beginning of each row in JAC,JA
c ... Workspace for  adifor with sparslink
      integer iwk(*)            ! integer work array,size>=3*ny+np
c === local variable
      integer i, j, ig_y, ig_yp, ig_delta, width, jj, itotal
      integer ig_senpar, nnz, nnzpj, j1, ny2, info
      character msg*80
      logical first
      save first
      data first/.true./
c
      ig_y      = 0
      ig_yp     = ig_y + ny
      ig_delta  = ig_yp + ny
      ig_senpar   = ig_delta + ny
      itotal    = ig_senpar + np
c
*      ipiwk = malloc(4*itotal)
      ny2 = 2*ny
      if (first) then
         call xspini        
         first = .false.
      end if
      do i = 1, ny
         call dspsd(iwk(ig_y + i), i, 1.0d0, 1)
         call dspsd(iwk(ig_yp + i), i+ny, 1.0d0, 1)
      end do
      do i = 1, np
         call dspsd(iwk(ig_senpar + i), i+ny2, 1.0d0, 1)
      end do
c
c     === call Adifor generated routine
      call i_res(
     $     t, y, iwk(ig_y+1), yprime, iwk(ig_yp+1), cj, delta, 
     $     iwk(ig_delta+1), ires, rpar, ipar, 
     $     SENPAR, iwk(ig_senpar+1))
c
c ... call error-handling routine
c      call ehrpt
c
      nnz = 0
      ia(1) = 1
      width = ny2 + np
      do i = 1, ny
         call dspxsq(ja(nnz+1),jac(nnz+1),width,
     *        iwk(ig_delta+i),jj,info)
         if (info .lt. 0) then
            msg = 'DASPK-- ERROR IN JRADFSP'
            CALL XERRWD(MSG,24,40,0,0,0,0,1,0.D0,0.0D0)
            ires = info
            RETURN
         end if
         j1 = 0
         do j = 1, jj
            nnzpj = nnz+j
            if (ja(nnzpj) .gt. ny2) then
               dfdp(i, ja(nnzpj)-ny2) = jac(nnzpj)
            else
               j1 = j
            end if
         end do
         nnz = nnz + j1
         ia(i+1) = nnz+1
      end do
      RETURN
C------------  End of Subroutine JrAdfSP  -------------------------------
      END
c
      SUBROUTINE amuvsp (n, x, y, a,ja,ia) 
      real*8  x(*), y(*), a(*) 
      integer n, ja(*), ia(*)
c-----------------------------------------------------------------------
c         A times a vector
c----------------------------------------------------------------------- 
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c----------
c n     = row dimension of A
c x     = real array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c-----------
c y     = real array of length n, containing the product y=Ax
c
c-----------------------------------------------------------------------
c local variables
c
      real*8 t
      integer i, k
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c 
         t = 0.0D0
         do 99 k=ia(i), ia(i+1)-1 
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux---------------------------------------------------
c-----------------------------------------------------------------------
      end
c
      SUBROUTINE amuv1sp (n, s, cj, a, ja, ia, y)
      real*8  s(*), cj, y(*), a(*) 
      integer n, ja(*), ia(*)
      real*8 t
      integer i, k, j
c-----------------------------------------------------------------------
c         y = cj* a(n+1:2n)*s + a(1:n)*s
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector sprime and s
c 
         t = 0.0d0
         do 99 k=ia(i), ia(i+1)-1
            j = ja(k)
            if (j .gt. n) then
               t = t + cj*a(k)*s(j-n)
            else 
               t = t + a(k)*s(j)
            end if
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux1---------------------------------------------------
c-----------------------------------------------------------------------
      end
      SUBROUTINE amuv2sp (n, s, sprime, a, ja, ia, y) 
      real*8  s(*), sprime(*), y(*), a(*) 
      integer n, ja(*), ia(*)
      real*8 t
      integer i, k, j
c-----------------------------------------------------------------------
      do 100 i = 1,n
c
c     compute the inner product of row i with vector sprime and s
c 
         t = 0.0D0
         do 99 k=ia(i), ia(i+1)-1
            j = ja(k)
            if (j .gt. n) then
               t = t + a(k)*sprime(j-n)
            else 
               t = t + a(k)*s(j)
            end if
 99      continue
c
c     store result in y(i) 
c
         y(i) = t
 100  continue
c
      return
c---------end-of-amux2---------------------------------------------------
c-----------------------------------------------------------------------
      end
      SUBROUTINE SD_JAC_ONLY(NY, A, JA, IA, CJ, WM, IWM, ier)
      IMPLICIT NONE
      INTEGER NY, JA(*), IA(*), IWM(*)
      REAL*8  CJ,a(*),wm(*)
      integer lml, lmu, lmtype, lnre, lnpd, llciwp
      PARAMETER (LML=1,LMU=2,LMTYPE=4,LNRE=12,LNPD=36,LLCIWP=30)

      integer i, j, k, jj, mu, ml, mband, meband, lipvt,ier,
     *     idense

      LIPVT = IWM(LLCIWP)
      ier = 0
      idense = 0
      if (iwm(lmtype) .gt. 4) idense = 1

      if (idense .eq. 0) then   ! dense matrix
         do i = 1, ny
            do j = 1, ny
               wm(i+(j-1)*ny) = 0.0d0
            end do
            do k = ia(i), ia(i+1)-1
               j = ja(k)
               if (j .gt. ny) then 
                  jj = j - ny
                  wm(i+(jj-1)*ny) = wm(i+(jj-1)*ny) + cj*a(k)
               else
                  wm(i+(j-1)*ny) = wm(i+(j-1)*ny) + a(k)
               end if
            end do
         end do
c
c     === dense factorization
c
         CALL DGEFA(WM,NY,NY,IWM(LIPVT),IER)
      else if (idense.eq. 1) then ! banded matrix
         MU = iwm(LMU)
         ML = iwm(LML)
         mband = mu + ml + 1
         do i = 1, ny            
            do j = 1, mband
               wm(i+j*mband-ml) = 0.0d0
            end do
            do k = ia(i), ia(i+1)-1
               j = ja(k)
               if (j .gt. ny) then 
                  jj = j - ny
                  wm(i+jj*mband-ML) = wm(i+jj*mband-ML) + cj*a(k)
               else
                  wm(i+j*mband-ML) = wm(i+j*mband-ML) + a(k)
               end if
            end do
         end do
         meband = mband + ML
c
c     === sparse factorization
c      
         CALL DGBFA (WM,MEBAND,NY,ML,MU,IWM(LIPVT),IER)
      end if
c
      return
      end
      SUBROUTINE adfsp_init
      logical first
      save first
      data first/.true./
      
      if (first) then  
c
c === initialize the sparselink 
c      
c====routine4adf
         call xspini
         first = .false.
      end if

      return
      end
c
      SUBROUTINE MPI_STEPSIZE(H)
      IMPLICIT NONE
      DOUBLE PRECISION H, HNEW
      INTEGER MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD
      PARAMETER (MPI_DOUBLE_PRECISION=27, MPI_MIN=101)
      PARAMETER (MPI_COMM_WORLD=91)
      INTEGER IERR
      HNEW = H
      CALL MPI_ALLREDUCE(H, HNEW, 1, MPI_DOUBLE_PRECISION,
     *     MPI_MIN, MPI_COMM_WORLD, IERR)
      H = HNEW                   ! new stepsize for MPI
      RETURN
      END
c
      SUBROUTINE FXPV(T,DY,CJ,DELTA,IRES,RPAR,IPAR,ADRES,NEQ,
     *     ISENFO,ISENWK,SENPAR,ADI,imeth)
c------------------------------------------------------------------
c   This routine calculates the matrix-vector product
c       delta = dy * F_y';
c
c   where  dy * F_y' is calculated by
c           ADRES
c
c   (c) Shengtai Li, 2000
c------------------------------------------------------------------
      implicit none
      integer ires, neq, isenfo(*), ipar(*), isenwk(*), imeth
      double precision t, dy(*), cj, delta(*), rpar(*),senpar(*),adi(*)
      external adres
c
c...  local variable
      integer ly, lyp, lad1, j, ny, ldelta
c
c...  routine body
c
c...  NY here is the number of state variables excluding the 
c     Quadrature variables from adjoint method
c     NQ here is the number of quadrature variables in original
c     forward equations.
      ny = neq - isenfo(4)
c
c...  pointers in ADI(*)
      ly   = 1
      lyp  = ly   + neq 
      ldelta = lyp + neq
      lad1 = ldelta + neq
      ires = 0
c     
c...  input F_y'
      if (imeth .eq. 0) then
c
c...  for the F_y'
         do j = 1, ny
            delta(j) = 0.0d0
            adi(lad1+j-1) = dy(j)
         end do 
      else if (imeth .eq. 1) then
c
c...  for F_y
c     
c...  AD-generated routine g_res for Y, dy*F_y
         if (isenfo(7) .eq. 0) then
            do j = 1, ny
               delta(j) = 0.0d0
               adi(lad1+j-1) = dy(j)
            end do
         else              
c     
c...  index-2 adjoint equation
            do j = 1, ny
               delta(j) = 0.0d0
               if (isenwk(2+j+neq) .eq. -2) then
                  adi(lad1+j-1) = dy(j)/cj
               else
                  adi(lad1+j-1) = dy(j)
               end if                    
            end do
         end if   
      end if
c     
c...  adres = resUSER_adyp is generated by TAMC, adyp*F_y'
      call adres(t, adi, adi(lyp), delta, cj, adi(ldelta),adi(lad1),
     &     ires, rpar, ipar, senpar)
      do j = ny+1, neq
         delta(j) = 0.0d0
      end do
c
      return
      end 
C
      SUBROUTINE KMVAD(
     *     T,DY,DYP,CJ,DELTA,IRES,RPAR,IPAR,ADRES,NEQ,
     *     ISENFO,SENWRK,ISENWK,SENPAR,G_RES,ADI)
c------------------------------------------------------------------
c   This routine calculates the matrix-vector product
c       delta = dy * F_y + dyp * F_y;
c
c   where dy * F_y is calculted by
c           G_RES
c   and dyp * F_y' is calculated by
c           ADRES
c
c   (c) Shengtai Li, 2000
c------------------------------------------------------------------
      implicit none
      integer ires, neq, isenfo(*), isenwk(*), ipar(*)
      double precision t, dy(*), dyp(*), cj, delta(*), 
     *     rpar(*), senpar(*), senwrk(*), adi(*)
      external adres, g_res
c
c...  local variable
      integer ly, lyp, Meva, lad1, lad2, nq, j, ny, ldelta
c
c...  routine body
      meva = isenfo(12)
c
c...  NY here is the number of state variables excluding the 
c     Quadrature variables from adjoint method
c     NQ here is the number of quadrature variables in original
c     forward equations.
      ny = neq - isenfo(4)
      nq = isenfo(10) - isenfo(4) ! original quadratures in the RES
c
c...  pointers in ADI(*)
      ly     = 1
      lyp    = ly  + neq 
      ldelta = lyp + neq
      lad1   = ldelta + neq
      lad2   = lad1 + neq
C     
C...  begin the evaluation of RES
      if (Meva .eq. 3 .or. Meva .eq. 7) then
c     
c...  Input by the user
         call adres(t,dy,dyp,cj,delta,ires,rpar,ipar,senpar,adi)
      else if (ires .eq. 0 .or. isenfo(10) .eq. 0) then
         ires = 0
         if (Meva .eq. 0 .or. Meva .eq. 4) then
c     
c...  AD-generated routine for Y only, ady*F_y
            do j = 1, ny
               delta(j) = 0.0d0
               adi(lad1+j-1) = dy(j)
            end do
            call g_res(t, adi, delta, adi(lyp), cj, 
     &           adi(ldelta), adi(lad1), ires, rpar, ipar, senpar)
            do j = 1, ny - nq
               delta(j) = dyp(j) + delta(j)
            end do
         else
c     
c...  AD-generated routine g_res for Y, dy*F_y
            if (isenfo(7) .eq. 0) then
               do j = 1, ny
                  delta(j) = 0.0d0
                  adi(lad1+j-1) = dy(j)
               end do
            else              
c     
c...  index-2 adjoint equation
               do j = 1, ny
                  delta(j) = 0.0d0
                  if (isenwk(2+j+neq) .eq. -2) then
                     adi(lad1+j-1) = dy(j)/cj
                  else
                     adi(lad1+j-1) = dy(j)
                  end if                    
               end do
            end if   
            call g_res(t, adi, delta, adi(lyp), cj, 
     &           adi(ldelta), adi(lad1), ires, rpar, ipar, senpar)
c     
c...  input F_y' 
            do j = 1, ny
               adi(lad1+j-1) = 0.0d0
               adi(lad2+j-1) = dyp(j)
            end do            
c     
c...  adres = resUSER_adyp is generated by TAMC, adyp*F_y'
            call adres(t, adi, adi(lyp), adi(lad1), cj, adi(ldelta), 
     &           adi(lad2), ires, rpar, ipar, senpar)
            do j = 1, ny - nq
               delta(j) =  delta(j) + adi(lad1+j-1)
            end do
         end if                 ! if ires == 0
      end if
c
      return
      end 
C
      SUBROUTINE ddresad(
     *     T,adY,adYP,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,
     *     adres, neq, isenfo, SENWRK, ISENWK, g_res, A_RES, adi)
      implicit none
      integer ires, neq, isenfo(*), isenwk(*), ipar(*)
      double precision t, ady(*), adyp(*), cj, delta(*), 
     *     rpar(*), senpar(*), senwrk(*), adi(*)
      external adres, g_res, A_RES
c
c...  local variable
      double precision g_t
      integer ltpre, lt, ly, lyp, Nt, i, Meva
      integer lad1, lad2, lg_adi, lg_x, ladyp, nq, j, ny
      integer neqt, ldelta
c
c...  routine body
      meva = isenfo(12)
c
c...  NY here is the number of state variables excluding the 
c     Quadrature variables from adjoint method
c     NQ here is the number of quadrature variables in original
c     forward equations.
      ny = neq - isenfo(4)
      nq = isenfo(10) - isenfo(4) ! original quadratures in the RES
      if (Meva .lt. 4 .and. isenfo(11).eq.0) then 
         ltpre = isenwk(2) + 1
c     
c...  Check and see if the interpolation is needed or not
c     If time has changed from last interpolation, redo it
c     Otherwise, using the old value
         if (t .ne. senwrk(ltpre)) then
            Nt      = isenwk(1)
            lt      = ltpre + 1
            ly      = lt    + Nt
            lyp     = ly    + Nt*Ny
c     
c...  do the interpolation at tf, compute new Y, YP at tf
            call dchi(
     *           t, NY, Nt, senwrk(lt), senwrk(ly), senwrk(lyp),
     *           adi, adi(neq+1))
            senwrk(ltpre) = t
         end if
      end if
      ly     = 1
      lyp    = ly   + neq
      ldelta = lyp + neq
      lad1   = ldelta  + neq
      lad2   = lad1 + neq
      lg_adi = lad2 + neq
C     
C...  begin the evaluation of RES
      if (Meva .eq. 3 .or. Meva .eq. 7) then
c     
c...  Input by the user
         call adres(T,adY,adYP,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,adi)
      else if (ires .eq. 0 .or. isenfo(10) .eq. 0) then
         ires = 0
         if (Meva .eq. 0 .or. Meva .eq. 4) then
c     
c...  AD-generated routine for Y only, ady*F_y
            do j = 1, ny
               delta(j) = 0.0d0
               adi(lad1+j-1) = ady(j)
            end do
            call g_res(t, adi, delta, adi(lyp), cj, 
     &           adi(ldelta), adi(lad1), ires, rpar, ipar, senpar)
            if (nq.gt.0) then
               ires = 3
               call g_res(t, adi, delta, adi(lyp), cj, 
     &              adi(ldelta), adi(lad1), ires,rpar,ipar,senpar)
               ires = 0
            end if
            do j = 1, ny - nq
               delta(j) = adyp(j) + delta(j)
            end do
         else
c     
c...  AD-generated routine g_res for Y, ady*F_y
            if (isenfo(7) .eq. 0) then
               do j = 1, ny
                  delta(j) = 0.0d0
                  adi(lad1+j-1) = ady(j)
               end do
            else              
c
c...  index-2 adjoint equation
               do j = 1, ny
                  delta(j) = 0.0d0
                  if (isenwk(2+j+neq) .eq. -2) then
                     adi(lad1+j-1) = ady(j)/cj
                  else
                     adi(lad1+j-1) = ady(j)
                  end if                    
               end do
            end if   
            call g_res(t, adi, delta, adi(lyp), cj, 
     *           adi(ldelta), adi(lad1), ires, rpar, ipar, senpar)
            if (nq.gt.0) then
               ires = 3
               call g_res(t, adi, delta, adi(lyp), cj, 
     *              adi(ldelta), adi(lad1),ires,rpar,ipar,senpar)
               ires = 0
            end if               
            if (Meva .eq. 1 .or. Meva .eq. 5) then
c
c...  F_y' is a constant matrix
               do j = 1, ny
                  adi(lad1+j-1) = 0.0d0
                  adi(lad2+j-1) = adyp(j)
               end do            
c
c...  adres = resUSER_adyp is generated by TAMC, adyp*F_y'
               call adres(t, adi, adi(lyp), adi(lad1), cj, 
     *              adi(ldelta), adi(lad2), ires, rpar, ipar,
     *              senpar)
               do j = 1, ny - nq
                  delta(j) =  delta(j) + adi(lad1+j-1)
               end do
            else if (Meva .eq. 2 .or. Meva .eq. 6) then
c     
c...  F_y' depends on either T or Y
               if (isenfo(11) .eq. 0) then
c
c...  during the integration
                  NEQT = NEQ*ISENFO(1)
c     
c...  first mu' + F_x^* \lambda = 0
                  do j = 1, ny-nq
                     delta(j) = adyp(j+neqt) + delta(j)
                  end do
c     
c...  second mu - F_{x'}^*\lambda = 0
                  do j = 1, ny
                     adi(lad1+j-1) = 0.0d0
                     adi(lad2+j-1) = ady(j)
                  end do            
c     
c...  adres = resUSER_adyp is generated by TAMC, ady*F_y'
                  call adres(t, adi, adi(lyp), adi(lad1), cj, 
     *                 adi(ldelta), adi(lad2), ires, rpar, ipar,
     *                 senpar)
                  do j = 1, ny - nq
                     delta(j+neqt) = ady(j+neqt) - adi(lad1+j-1)
                  end do
               else
c
c...  during the initialization
                  g_t = 1.0d0
                  ladyp = lg_adi + neq
                  do j = 1, ny
                     adi(lad1 + j - 1) = 0.0d0
                     adi(lad2 + j - 1) = ady(j)
                     adi(lg_adi+j - 1) = 0.0d0
                     adi(ladyp +j - 1) = - adyp(j)
                  end do
c     
c...  adres = g_resUSER_adyp is generated by ADIFOR and TAMC,
c     adyp*F_y' - ady*(dF_y'/dt)
                  call adres(t, g_t, adi, adi(lyp), adi(lyp), 
     *                 adi(lad1), adi(lg_adi), cj, adi(ldelta),
     *                 adi(lad2), adi(ladyp),ires, rpar, ipar, senpar)  
                  do j = 1, ny - nq
                     delta(j) =  delta(j) - adi(lg_adi+j-1)
                  end do 
               end if
            end if            
         end if                 ! if ires == 0
      else if (ires .eq. 3) then
c
c...  residual for the original quadratures
         do i = 0, nq-1
            delta(ny-i) = adyp(ny-i)
         end do
c     
c...   compute \dot c = ady*F_p
         if (isenfo(4) .gt. 0) then
            do i = 1, ny
               adi(lad1+i-1) = ady(i)
            enddo
            do i = 1, isenfo(4)
               delta(ny+i) = 0.0d0
            end do
            ires = 0
            call A_RES(t, adi, adi(lyp), cj, adi(ldelta), adi(lad1),
     *           ires, rpar, ipar, senpar, delta(ny+1))
            ires = 3
            do i = 1, isenfo(4)
               delta(ny +i) = adyp(ny+i) - delta(ny+i)
            end do
         end if
      end if
c     
c...  for the initialization where g_x is required.
      if (isenfo(11) .eq. 1) then
c     
c...  -isenfo(6) save the index of the derived function
         lg_x = 6*neq + (-isenfo(6)-1)*ny
         if (meva .eq. 2) lg_x = lg_x + neq
         do i = 1, ny
            delta(i) = delta(i) - adi(lg_x + i)
         end do
      end if
      if (isenfo(7) .eq. 1 .AND. (Meva.ne.3.and.Meva.ne.6)) then
c
c...  multiply CJ to index-2 constraints
         do j = 1, ny
            if (isenwk(2 + j) .eq. -2) then
               delta(j) = delta(j) * cj
            end if
         end do
      end if
      return
      end
c
c...  Cubic Hermit Interpolation
      SUBROUTINE dchi(t, NY, Nt, ta, ya, ypa, y, yp)
      implicit none
      integer NY, Nt, i, j
      double precision t, ta(Nt), ya(NY,Nt), ypa(NY,Nt), y(NY), yp(*)
      double precision dy, y1tl, y1tr, c1, c2, c3, c4, qdt, dif, dt
      i = 1
      do while (t .gt. ta(i) .and. i .lt. Nt)
         i = i + 1
      end do
      i = i - 1
      if (i .eq. 0) i = 1
c
      dt = ta(i+1) - ta(i)
      qdt = 1.0d0/dt
      dif = t - ta(i)
      do j = 1, NY
         dy = (ya(j,i+1) - ya(j,i))*qdt
         y1tl = ypa(j,i)
         y1tr = ypa(j,i+1)
         c1 = ya(j, i)
         c2 = y1tl
         c3 = (3.0d0*dy - y1tr - 2.0d0*y1tl)*qdt
         c4 = (y1tr + y1tl - 2.0d0 * dy) * qdt * qdt
         y(j) = c1 + dif * (c2 + dif*(c3 + dif*c4))
         yp(j)= c2 + dif * (2.0d0 * c3 + dif * 3.0d0 * c4)
      end do
      return
      end
      SUBROUTINE resadidx2(
     *     T,adY,adYP,CJ,DELTA,IRESI,RPAR,IPAR,SENPAR,
     *     adres, neq, ide, isenfo, SENWRK, ISENWK, g_res, adi)
      implicit none
      integer iresi, neq, isenfo(*),isenwk(*),ipar(*),ide(*)
      double precision t, ady(*), adyp(*), cj, delta(*), 
     *     rpar(*), senpar(*), senwrk(*), adi(*)
      external adres, g_res
c
c...  local variable
      double precision g_t
      integer ly, lyp, i, ires, Meva
      integer lad1, lad2, lg_adi, lg_x, nq, j, ny, ladyp, ldelta
      meva = isenfo(12)
c
c...  NY here is the number of state variables excluding the 
c     Quadrature variables from adjoint method
c     NQ here is the number of quadrature variables in original
c     forward equations.
      ny = neq - isenfo(4)
      nq = isenfo(10) - isenfo(4) ! original quadratures in the RES
      ly     = 1
      lyp    = ly   + neq 
      ldelta = lyp + neq
      lad1   = ldelta + neq
      lad2   = lad1 + neq
      lg_adi = lad2 + neq
      ladyp  = lg_adi + neq
C     
C...  begin the evaluation of RES
      if (Meva .eq. 3 .or. Meva .eq. 7) then
         ires = 2
c     
c...  Input by the user
         call adres(T,adY,adYP,CJ,DELTA,IRES,RPAR,IPAR,SENPAR,adi)
      else 
         ires = 0
c
c...  F_y' depends on either T or Y
         g_t = 1.0d0
         do j = 1, ny
            adi(lad1 + j - 1) = 0.0d0
            adi(lg_adi+j - 1) = 0.0d0
            adi(lad2 + j - 1) = ady(j)
            adi(ladyp + j - 1) = -adyp(j)
         end do
c
c...  adres = g_resUSER_ady is generated by ADIFOR and TAMC
         call g_res(t, g_t, adi, adi(lyp), adi(lad1), adi(lg_adi), 
     *        adi(lyp), cj, adi(ldelta), adi(lad2), adi(ladyp), ires, 
     *        rpar, ipar, senpar)
         do j = 1, ny
            if (ide(j) .eq. 1) then
               delta(j) =  -adi(lg_adi+j-1)
            end if
         end do 
      end if
c     
c...  for the initialization where g_x is required.
      if (isenfo(11) .eq. 3) then
c     
c...  -isenfo(6) save the index of the derived function
         lg_x = 7*neq + (-isenfo(6)-1)*ny
         do i = 1, ny
            if (ide(i) .eq. 1) then
               delta(i) = delta(i) - adi(lg_x + i)
            end if
         end do
      else if (isenfo(11) .eq. 4) then
         lg_x = 7*neq + (-isenfo(6)-1)*ny
         do i = 1, ny
            if (ide(i) .eq. 0) then
               delta(i) = delta(i) - adi(lg_x + i)
            end if
         end do
      end if
      return
      end
c
c...  Preconditioner for the Krylov method and adjoint mode
      subroutine jackad(
     *     jack, adres, ires, neq, t, y, yprime, rewt, savr,
     *     wk, h, cj, wp, iwp, ierr, rpar, ipar, senpar, adi,
     *     isenfo, senwrk, isenwk, g_res, p_res, icopt, id)
      IMPLICIT NONE

C ... Input arguments:
      INTEGER NEQ               ! total number of equations
      REAL*8 T                  ! independent variable t
      REAL*8 Y(NEQ)             ! most recent iterate of solution vector y
      REAL*8 YPRIME(NEQ)        ! most recent iterate of solution vector y'
      REAL*8 SAVR(NEQ)          ! current residual evaluated at (T,Y,YPRIME)
      REAL*8 REWT(NEQ)          ! scale factors for Y and YPRIME
      INTEGER IRES              ! error flag for RES routine
      REAL*8 WK(NEQ)            ! work space available to this subroutine
      REAL*8 H                  ! current step size
      REAL*8 CJ                 ! scalar proportional to 1/H
      REAL*8 RPAR(*)            ! user real workspace
      INTEGER IPAR(*)           ! user integer workspace
      real*8 SENPAR(*)          ! sensitivity parameter array

C ... Output arguments:
      REAL*8 WP(*)              ! matrix elements of ILU
      INTEGER IWP(*)            ! array indices for elements of ILU
      REAL*8 ADI(*)             ! forward state information 
      INTEGER IERR              ! error flag (0 means success, else failure)
      REAL*8  SENWRK(*)         ! work space 
      INTEGER ISENFO(*), ISENWK(*)
      INTEGER ICOPT, ID(*)
      EXTERNAL JACK, ADRES, G_RES, P_RES

      call JACK(ADRES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
     *     WK, H, CJ, WP, IWP, IERR, RPAR, IPAR, SENPAR, ADI, 
     *     ISENFO,SENWRK,ISENWK,G_RES, P_RES, ICOPT, ID)

      return
      end



