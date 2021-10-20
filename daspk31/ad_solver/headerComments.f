c======================================================================
C***BEGIN PROLOGUE  DDASPKADJOINT
C***DATE WRITTEN   991207   (YYMMDD)
C***REVISION DATE  000620   (Release version 1.0)
C***REVISION DATE  000723   (Bug fixed in the output)
C***REVISION DATE  000726   (Bug fixed in the work space)
C***REVISION DATE  010501   (Added one pass option)
C***REVISION DATE  020506   (Bug fixed for info(16)=0 and index>1 problem)
C***REVISION DATE  020517   (Bug fixed for work space problem)
C***REVISION DATE  020801   (changed argument list for AD-generated routine)
C***REVISION DATE  020801   (Added INFOBI(8) for Adifor3.0 adjoint mode)
C***REVISION DATE  020802   (Added script file for Adifor3.0 adjoint mode)
C***REVISION DATE  020813   (Bug fixed for SENPAR(*)=0)
C***REVISION DATE  060701   (Added support for intermediate values of adjoint variable)
C***REVISION DATE  061401   (Bug fixed by setting LNLIS=38 instead of 39)
C***REVISION DATE  062505   (Bug fixed by using LTPRE instead of LTPRE+LRWORKB)
C***REVISION DATE  070131   (Bug fixed for initial interpolation time)
C***REVISION DATE  070306   (Updated comments to reflect additional features)
C   
C***AUTHORS  Shengtai Li and Linda Petzold
C            with additional bug-fixes and features by Stephanie R Taylor
C
C            Department of Computer Science
C            University of California 
C            Santa Barbara, CA  93106
C            Email: petzold@engineering.ucsb.edu
C                   shengtai@engineering.ucsb.edu
C                   staylor@cs.ucsb.edu
C
C***PURPOSE  For a system of differential/algebraic equations (DAEs) 
C            with sensitivity parameters p
C                F(t, y, y', p) = 0,      (*)
C            and initial conditions 
C                Y_0 = Y0(p),    
C            this code solves for the sensitivities of the derived functions 
C                g = g(T, y, p) 
C            with respect to p at time T by solving the adjoint system
C            backwardly with a combination of Backward Differentiation 
C            Formula (BDF) methods and a choice of two linear system 
C            solution methods: direct (dense or banded) or Krylov 
C            (iterative). The class of problems that DASPKADJOINT can 
C            solve includes index-0, index-1 and index-2 DAE system with
C            nonsingular mass matrix. This version is in double precision.
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C *Usage:
C
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR(*)
C      DOUBLE PRECISION T, Y(*), YPRIME(*), TOUT, TFINAL, RTOL(*), ATOL(*),
C         RWORK(LRW), RPAR(*), SENPAR(*)
C      EXTERNAL  RES, JAC, PSOL, G_RES, ADRES, ADJAC, RES_ADP, 
C     *          RES_ADY, G_RES_ADY, ADINIT, K_RES, T_RES
C
C      CALL DDASPKADJOINT (RES, NEQ, T, Y, YPRIME, TOUT, TFINAL, INFO, RTOL, ATOL,
C     *  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL, SENPAR,
C     *  ADRES, NQ, QRES, QSEN, INFOB, RTOLB, ATOLB, NBUF, ADJAC, IEOPT,
C     *  RES_ADP, RES_ADY, G_RES_ADY, ADINIT, K_RES, T_RES) 
C
C    Quantities which may be altered by the code are:
C       T, Y(*), YPRIME(*), RTOL, ATOL, IDID, RWORK(*), IWORK(*), QSEN(*)
C
C *Arguments:
C   
C  RES:EXT          This is the name of a routine which you
C                   provide to define the residual function F(t,y,y',p)
C                   of the differential/algebraic system.
C
C  NEQ:IN           This is the number of state equations in the system.
C
C  T:INOUT          This is the current value of the independent 
C                   variable.
C
C  Y(*):INOUT       This array contains the solution components at T.
C
C  YPRIME(*):INOUT  This array contains the derivatives of the solution
C                   components at T.
C
C  TOUT:IN          This is a point at which a solution is desired (but 
C                   this cannot include the sensitivity values).
C
C  TFINAL:IN        This is the point at which the solution is desired and
C                   (possibly) sensitivities are desired.
C
C  INFOI(N):IN      This is an integer array used to communicate details
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
C                   program and RES, JAC, and PSOL routines.
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
C  SENPAR:IN        This is a real array of sensitivity parameters that 
C                   appear in the RES routine. If you wish to compute 
C                   sensitivities and RES depends on the problem parameters, 
C                   SENPAR must be used to store and pass the problem 
C                   parameters to RES.
C  
C  ADRES:EXT        This is the name of a routine which you must provide 
C                   to define the residual function for the adjoint 
C                   differential/algebraic system if you select the 
C                   option to construct the adjoint DAE by yourself, 
C                   or a routine which you must provide for the vector-
C                   matrix product vF_y' if you select the option to 
C                   construct the adjoint system by DASPKADJOINT and F_y'
C                   is not the identity.
C
C  NQ:IN            Total number of derived funtions.
C
C  QRES:EXT         This is the name of a subroutine which you must
C                   provide to define the derived quantities g(T,y,p)
C                   and/or their derivatives with respect to y and p.
C
C  QSEN(*):INOUT    This array contains the values of the derived 
C                   quantities and their derivatives with respect to y and p.
C                   It also contains the sensitivity information
C                   at output.
C
C  INFOBI(N):IN     This is an integer array used to communicate details
C                   of how the adjoint solution is to be carried out.
C                   N must be at least 10.
C                   
C  RTOLB,ATOLB:INOUT  These quantities represent absolute and relative
C                   error tolerances (on local error) which you provide
C                   to indicate how accurately you wish the adjoint 
C                   solution to be computed.  You may choose them to be 
C                   both scalars or else both arrays of length NEQ.
C
C  NBUF:IN          This is the size of the buffers (double precision) 
C                   used to communicate between the forward integration and 
C                   the backward (adjoint) integration.
C
C  ADJAC:EXT        This is the name of a routine which you may
C                   provide (optionally) for calculating Jacobian 
C                   (partial derivative) data involved in solving linear
C                   systems within DDASPKADJOINT. 
C                   When you use the ADIFOR option to calculate the Jacobian,
C                   this is the name of the ADIFOR-generated routine.
C
C  IEOPT(*):IN      This integer array indicates the properties of each
C                   equation and each variable.
C
C  RES_ADP:EXT      This is the name of a routine which you may provide
C                   (optionally) for calculating the vector-matrix
C                   product vF_p.
C
C  RES_ADY:EXT      This is the name of a routine which you may provide
C                   (optionally) for calculating the vector-matrix
C                   product vF_y.
C
C  NEQAD:IN         An integer specifying the size of the adjoint variables.
C
C  ADY(NEQAD):OUT   Contains the adjoint variable. If INFO(3) = 0, then ADY may
C                   be of length 0 (but NEQAD must still be specified correctly).
C
C  ADYP(NEQAD):OUT  Contains the adjoint variable associated with the derivative.
C                   If INFO(3) = 0, then ADY may be of length 0 (but NEQAD must
C                   still be specified correctly).
C
C  G_RES_ADY:EXT    This is the name of a routine which you may provide 
C                   (optionally) for initialization of an index-2 adjoint
C                   system.
C
C  G_RES_ADYP:EXT   This is the name of a routine which you may provide 
C                   (optionally) for initialization when F_y' is time-variant.
C
C  ADINIT:EXT        This is the name of a routine which you may provide 
C                   (optionally) for the initialization of the adjoint system.
C 
C  K_RES:EXT        This is the name of the ADIFOR-generated routine which 
C                   you may provide (optionally) for the matrix-vector product 
C                   in the Krylov iteration. 
C
C  T_RES:EXT        This is the name of the ADIFOR-generated routine which 
C                   you may provide (optionally) for the initialization of 
C                   an index-2 system.
C
C  G_RES_ADP:EXT    This is the name of a routine which you may provide 
C                   (optionally) for the sensitivity evaluation of an index-2
C                   system.
C
C *Overview 
C
C     For a DAE system (DAEs) with sensitivity parameters p
C          F(t, y, y', p) = 0,      (*)
C     and initial conditions 
C          Y_0 = Y0(p),    
C     DASPKADJOINT is used to compute the sensitivity of the derived 
C     functions (which are defined by QRES)
C          g = g(T, y, p) 
C     with respect to p at time T.
C
C     DASPKADJOINT is designed to be used in conjunction with DASPK3.0 as 
C     described below. Many of the arguments are the same and they will not 
C     be described here. The user is advised to read carefully the documentation
C     of DASPK3.0.
C
C     DASPKADJOINT does not support directly those cases where the derived  
C     function is of integral form, i.e.,
C          G = \int^{T}_{0} g(tau, y, p) d (tau),
C     where \int^T_0 denotes the integration operation from 0 to T.
C     However, G can be represented as a quadrature variable in DASPK3.0
C     and then G = y_{N} is of the required form.
C
C     DASPKADJOINT requires that a sensitivity parameter appear either in the 
C     initial conditions or in the DAE(RES subroutine). Sensitivity parameters
C     appearing in the DAE must be stored in the SENPAR array.
C
C     The adjoint system for a given DAE F(t, y, y', p) = 0 is of the form
C
C         (F_y')^*ADYPRIME - ((F_y)^* - ((F_y')_t)^*) ADY = 0  (**)
C
C     where ADY and ADYPRIME are the adjoint variables and their time
C     derivatives, (A)^* denotes the adjoint matrix of A, and (F_y')_t denotes
C     the total derivative of F_y' with respect to time t. The sensitivity
C     of g = g(T, y, p) with respect to p is calculated by
C
C         (g_p - ADY*F_p)|t=T - \int^T_0 (ADY*F_p) + (ADY*F_y')|t=0 * (Y_0)_p.
C                                                              (***)
C
C     DASPKADJOINT solves the adjoint DAEs (**) backwards for the variables
C     ADY. After a time reversing transformation, DASPKADJOINT actually solves
C
C         (F_y')^*ADYPRIME + ((F_y)^* - ((F_y')_t)^*) ADY = 0  (****)
C
C     forward by DASPK3.0. Note that F_y' and F_y depend on Y for a nonlinear
C     problem. Thus we must compute Y first before constructing the residual 
C     routine for the DAE (****). DASPKADJOINT does this by calling DASPK3.0. 
C     A two level checkpointing with cubic Hermite interpolation is used inside 
C     DASPKADJOINT to reduce the required storage in cases where the user defined 
C     BUFFER is not large enough. A checkpoint in DASPKADJOINT is defined as a 
C     specific time at which ENOUGH information is saved for later computation. 
C     The 'ENOUGH' means that if we start DASPK at the checkpoint with the saved 
C     information, it should produce the same results.
C      
C     Inside DASPKADJOINT, we first calculate two constants according to the size
C     of the BUFFER specified in argument NBUF by the user: NCHKSTP, which is the
C     number of time steps allowed between two consecutive checkpoints, and NCHK,
C     which is the number of checkpoints the BUFFER could allow.
C     Next we solve the original DAE (*) forward by DASPK3.0, and save enough 
C     information at each checkpoint in the BUFFER. After 'NCHK' checkpoints, 
C     the BUFFER is full. Then we write the whole BUFFER to a DISK file and 
C     reuse the BUFFER for later checkpoints. 
C
C     After finishing the forward integration, we compute the consistent 
C     initial conditions for the adjoint system at T = TFINAL. Then we solve 
C     the adjoint system backwards by DASPK3.0. If Y is needed in constructing
C     the adjoint system, we integrate the original system forward first 
C     between two checkpoints by retrieving the data from either the BUFFER 
C     or the DISK file, and save the values of Y and YPRIME at each time step 
C     in a work array. The required values of Y during the backward 
C     integration are obtained by cubic Hermite interpolation in time. 
C     
C     This version handles only those DAE problems where the index-1 and/or
C     index-2 variables can be identified explicitly, for example semi-explicit
C     Hessenberg form DAEs.
C
C *Note*
C    For any EXTERNAL routines that you do not need, do not declare them
C    as EXTERNAL when calling DASPKADJOINT. Instead, declare them as 
C    integer or double precision variables, or just omit them if they are
C    at the end of the calling list. In this way, you do not need to 
C    construct dummy routines for these unused external routines.
C
C    Constructing the adjoint equation inside DASPKADJOINT requires an
C    automatic differentiation (AD) tool with reverse mode. 
C    By the time this code is written, ADIFOR3.0 is not yet released.
C    We suggest you use Tapenade whenever you need to use a reverse
C    mode AD tool. A local installation of Tapenade can be obtained
C    freely from the website http://www-sop.inria.fr/tropics/tapenade.html.
C    See the README.TXT file in the examples/DM directory for
C    notes about installation. There is an alternative to Tapenade,
C    TAMC (Tangent linear and adjoint model compiler) by
C    Ralf Giering (http://puddle.mit.edu/~ralf/tamc). However, local
C    installations of TAMC are no longer freely available. It is
C    possible to use a program that sends email to a TAMC engine. 
C    We suggest using the `-pure' option when you invoke TAMC. 
C    The code generated by Tapenade and TAMC is usually not in the format that 
C    DASPKADJOINT requires. You need to change the argument list in the TAMC/Tapenade 
C    generated routine or construct another routine which is in our format and 
C    calls the TAMC/Tapenade generated routine. 
C
C    For problems where a forward mode AD tool is required, we recommend 
C    using ADIFOR 2.0 because ADIFOR-generated code is in the required format 
C    and does not need user modification.
C
C------INPUT - WHAT TO DO ON THE FIRST CALL TO DDASPKADJOINT-------------------
C
C *Arguments*: 
C
C    (Here we describe only those arguments which may have a different meaning
C     from DDASPK3.0. For those variables which have the same meaning as in 
C     DDASPK3.0, please refer to the documentation of DDASPK3.0 for more details.)
C
C  RES -- RES has the same meaning as in DDASPK3.0 except that the user
C         does not need to define the sensitivity equations.
C
C  NEQ -- Set it to the number of equations in the system. This is 
C         the number of equations for the state variables, excluding
C         the sensitivity equations.
C
C  T   -- Set it to the initial point of the forward integration.
C 
C  Y(*) -- Set this array to the initial values of the solution.
C          You must dimension Y of length at least NEQ in your calling program.
C
C          If some sensitivity parameters appear in the initial conditions 
C          and you want DASPKADJOINT to calculate the sensitivity (by 
C          specifying INFOBI(6) = 1), you must also set the CONSISTENT 
C          initial values--(Y_0)_p--of the sensitivity components at the 
C          initial point. In this case, you must dimension Y of length at 
C          least NEQ*(1+INFO(19)) in your calling program.
C
C  YPRIME(*) -- Set this array to the initial values of the first
C               derivatives of the solution. You must dimension YPRIME 
C               at least NEQ in your calling program. 
C
C  TOUT -- Set it to the first output point at which the solution is desired.
C          Note: This parameter is used only when INFOI(3) = 2.
C
C  TFINAL -- Set it to the output point at which the sensitivity is desired.
C
C  INFO(*) -- The INFO array has the same meaning as in DDASPK3.0. Please
C             refer to the documentation of DASPK3.0 on how to set it.
C
C             In the following we describe only those components of INFO that 
C             have changed from the DASPK3.0 documentation. 
C
C       INFO(2) - How much accuracy you want of your solution
C                 is specified by the error tolerances RTOL and ATOL (tolerance
C                 for the state variables), RTOLB and ATOLB (tolerance for
C                 the adjoint variables).
C 
C                 The simplest use is to take them to be scalars.
C                 To obtain more flexibility, they can both be arrays.
C                 The code must be told your choice.
C
C            **** Are both error tolerances RTOL, ATOL scalars ...
C                 yes - set INFO(2) = 0
C                       and input scalars for RTOL and ATOL, RTOLB and ATOLB.
C                  no - set INFO(2) = 1
C                       and input arrays of length NEQ for
C                       RTOL and ATOL. Dimension RTOLB and ATOLB
C                       of length at least NEQAD=NQ*NEQPNS where NQ is the
C                       number of derived functions and NEQPNS is the number
C                       of states in the original system plus the number of
C                       parameters on which sensitivity analysis is performed. 
C                       Each entry in RTOLB and ATOLB corresonds to an element 
C                       ADY, e.g. the element associated with dQI/dYJ is at index 
C                       (NEQPNS-1)*I+J. For more information, see below for the 
C                       discussion of the output value ADY. ****
C
C       INFO(3) - The overall goal of DASPKAdjoint is to compute the 
C              sensitivities of the system at time TFINAL.  To do this
C              it performs a forward integration of the system from
C              time T to time TFINAL.  Then, the adjoint system is
C              solved by integrating backwards from time TFINAL to time
C              T.  If you wish, you may retrieve the values of Y, YPRIME,
C              ADY, and ADYPRIME at intermediate steps during the _backwards_ 
C              integration.
C   
C
C          **** Do you want the solution only at
C               TFINAL (i.e. do you want DASPKAdjoint to run the entire
C                       forward and entire backward integration before it
C                       returns?) ...
C                yes - set INFO(3) = 0
C
C                 no - there are two options:
C                       1) Do you want DASPKAdjoint to return Y,YPRIME,ADY,ADYPRIME
C                          at the next intermediate step of the backward integration?
C                          (this is similar to intermediate-step mode in DASPK).
C
C                       set INFO(3) = 1
C
C                       In this case TOUT will be ignored.  At each intermediate step,
C                       values will be returned in T (the current time), Y, YPRIME, 
C                       ADY, and ADYPRIME.  ****
C
C                      2) Do you want DASPKAdjoint to return Y,YPRIME,ADY,ADYPRIME
C                         at a set output time TOUT?
C
C                       set INFO(3) = 2 
C
C                       On the first call to DASPKAdjoint, you should set T to be the
C                       starttime of the integration (e.g. T=0), TFINAL to the
C                       stoptime of the integration, and TOUT to the time at which you want
C                       the first output.  In this case, T < TOUT <= TFINAL.
C                   
C                       On any subsequent call to DASPKAdjoint, T will be the current time
C                       in the backward integration.  (i.e. the TOUT from the previous call).
C                       In this case you should set TOUT < T.
C
C       INFO(19) -- Number of sensitivity parameters.
C                 INFO(19) includes the number of parameters in the initial 
C                 conditions and the number of parameters in the equations 
C                 (specified in INFO(22)). 
C                 
C       INFO(20) -- options for obtaining the derivatives g_y and g_p for
C                 the derived function g = g(T, y, p). 
C
C         ****   Do you wish to use a finite difference method (FDM)
C                for the approximation of the derivatives of the derived 
C                functions?
C                Yes - Set INFO(20) = 0 or 1.
C                No  - There are two options. See the description of QRES in
C                      this documentation for more details.
C 
C       INFO(22) -- Number of parameters that appear in RES. If INFO(20) is 
C                   NOT equal to 2, INFO(22) must be equal to or greater than 
C                   the dimension of SENPAR(*).
C
C  RTOL, ATOL -- See the documentation of DDASPK3.0.
C 
C  RWORK(*) -- a real work array, which should be dimensioned in your
C              calling program with a length equal to the value of
C              LRW (or greater). 
C
C              The RWORK(*) in DASPKADJOINT is not only a real work array 
C              for the forward integration, but also a real work array
C              for the backward integration. It also serves as the BUFFER
C              to communicate between the forward and backward integration.
C
C  LRW -- Set it to the declared length of the RWORK array. The
C              maximum length depends on the options you have selected.
C 
C              If the system is linear or no checkpointing is needed 
C              (i.e. the system is nonlinear but INFOBI(7) is nonzero),
C              then
C
C                  BASE = max(LRWF,LRWB)
C
C              Otherwise,
C
C                  BASE = LRWF + LRWB + NBUF
C
C              where NBUF is an argument specified by the user (see the 
C              description of NBUF in this documentation for details), LRWF 
C              is the length of the real work array for the forward 
C              integration without considering sensitivity (i.e. there are 
C              NY equations), and LRWB is the length of the real work array 
C              for the backward integration (which uses NEQAD equations).
C
C              If the mass matrix F_y' is constant, the total length of 
C              RWORK(*) is 
C
C                 LRW = BASE + nsenpar*nq + 2*neqad + 1 + 2*neq + RCONTLEN
C
C              0therwise, if F_y' depends on t or y, then
C
C                 LRW = BASE + nsenpar*nq + 4*neqad + 1 + 2*neq + RCONTLEN
C
C              where NY=NEQ is the number of state variables, NQ is the 
C              number of derived functions, NSENPAR is the number of parameters
C              used in the sensitivity analysis, NEQAD is NQ*(NEQ+NSENPAR), 
C              and RCONTLEN is 10.
C
C              Please see the description of LRW in the 
C              documentation of DASPK3.0 for the length LRWF. 
C
C  IWORK(*) -- An integer work array, which should be dimensioned in
C              your calling program with a length equal to the value
C              of LIW (or greater).
C
C  LIW -- Set it to the declared length of the IWORK array.
C             
C             LIW must be at least LIWF + LIWB + NBUF/(3*NEQ) + 3*NEQ + 90, 
C             where LIWF is the length of the integer work space required by 
C             the forward integration without considering sensitivity and
C             LIWB is the length of the integer work space required by the
C             backward integration (which involved NEQAD equations).
C
C  RPAR, IPAR -- See the documentation of DDASPK3.0.
C
C  JAC -- This is  the name of a routine you may supply
C         (optionally) that relates to the Jacobian matrix of the
C         nonlinear system during the forward integration. 
C         See the documentation of DDASPK3.0 for more details.
C
C  PSOL -- This is the name of a routine you must supply if you have
C          selected the Krylov method (INFO(12) = 1) with preconditioning.
C          See the documentation of DDASPK3.0 for more details.
C                    
C          DDASPKADJOINT assumes that the PSOL for the backward integration
C          is the same as for the forward integration.
C
C  SENPAR -- This is an array of sensitivity parameters that appear in 
C            RES. They are not altered by DDASPKAJOINT. 
C  
C  ADRES:EXT -- This is the name of a subroutine which you
C              provide to define the residual function for the 
C              adjoint differential/algebraic system when
C              INFOBI(3) = 3 or 7, or the name of the AD-generated routine 
C              RES_ADYPRIME when INFOBI(3) = 1,2,5,6.
C              (see the description of the item INFOBI(3) in this
C               documentation).
C
C              If the initialization of an index-2 DAE system is required, 
C              ADRES should evaluate the time derivative of the algebraic 
C              constraints for the adjoint system if IRES=2 and 
C              INFOBI(3) = 3 or 7.
C
C              If INFOBI(3) = 3 or 7, ADRES must have the form 
C
C              SUBROUTINE ADRES(T, ADY, ADYPRIME, CJ, DELTA, IRES, RPAR, 
C            *                  IPAR, SENPAR, Y)
C
C              where Y is the vector of state variables at T.
C
C  NQ -- Total number of derived functions. 
C                
C  QRES:EXT -- This is the name of a subroutine which you must
C              provide to define the derived quantity g(T,y,p).
C              QRES must have the form
C
C              SUBROUTINE QRES(T, Y, YPRIME, QSEN, IRES, RPAR, IPAR, SENPAR)
C
C              to define the derived quantity QSEN = g(T, y, p).
C                 
C              If INFO(20) = 0 or 1 (finite difference is selected), QRES
C                 defines only QSEN;
C
C              If INFO(20) = 2 (user-input is selected), QRES defines
C                 both QSEN and the gradients of QSEN with respect
C                 to Y and p.
C                 The first NQ elements contain the values of g.
C                 The elements from NQ+1 to NQ*(NEQ+1) contain the gradient
C                    of g with respect to Y. The derivative of the ith function
C                    with respect to the jth element in Y is stored at
C                    NQ + (j-1)*NQ + i.
C                 The elements from NQ*(NEQ+1)+1 to NQ*(NEQ+1+INFO(22)) 
C                    contain the gradient of g with respect to p. The 
C                    derivative of the ith function with respect to the jth 
C                    element in p is stored at NQ*(NEQ+1) + (j-1)*NQ + i.
C                 For example, if there are NEQ=2, NQ=2, and NP=2, then QSEN 
C                    has the form
C                        Y1
C                        Y2
C                        dg1/dY1
C                        dg2/dY1
C                        dg1/dY2
C                        dg2/dY2
C                        dg1/dP1
C                        dg2/dP1
C                        dg1/dP2
C                        dg2/dP2
C
C              If INFO(20) = 3 or 4 (ADIFOR is selected), QRES is
C                 the name of an ADIFOR-generated routine for defining QSEN and
C                 the gradients of QSEN with respect to Y and p. It must have 
C                 the form
C
C                 SUBROUTINE QRES(NP,T,Y,G_Y,QSEN,G_QSEN,IRES,
C                *                RPAR,IPAR,SENPAR,G_SENPAR)        
C 
C                 SUBROUTINE QRES(NP,T,Y,G_Y,QSEN,G_QSEN,IRES,
C                *                RPAR,IPAR,SENPAR,G_SENPAR)        
C
C                 Steps on how to generate the routine by ADIFOR:
C                     1. Put all the routines related to 
C
C                   SUBROUTINE QRES(T, Y, YPRIME, Q, IRES, RPAR, IPAR, SENPAR)
C
C                        in a file called "qres.f"
C
C                     2. Create file "qres.cmp" with one line:
C                        qres.f
C
C                     3. Create file "qres.adf":
C                            AD_PROG = qres.cmp
C                            AD_TOP = qres
C                            AD_IVARS = y,senpar
C                            AD_OVARS = Q
C                            AD_PREFIX = g
C                            AD_OUTPUT_DIR = jac
C                            AD_SUPPRESS_LDG = true
C                            AD_EXCEPTION_FLAVOR=performance 
C
C                     4. Run Adifor to generate g_QRES(...) with
C                        % Adifor AD_SCRIPT=qres.adf
C
C                     5. Use g_QRES when calling DASPK_AJOINT instead of
C                        QRES.    
C                If ADIfor3.0 (or newer version ADIFOR) is used, the script
C                "spec.ad3" can be written as
C
C                   Function
C                     JM
C                   MaxColsM
C                     a number >= NEQ + INFO(22)
C                   Independent
C                     y
C                     senpar
C                   Dependent
C                     qsen
C                   Top
C                     qres
C                   ExceptionDefaultMode
C                     Performance
C                   Cache
C                     JM_Cache
C                   OutputDir
C                     JM_Output  
C
C                 and then run 
C                   %Adifor3.0 -spec spec.ad3 qres.f
C
C  QSEN(*) -- This array of dimension NQ*(NEQ+1+INFO(19))
C             will contain the values of the derived function g and/or 
C             the gradients of g with respect to Y and parameter p, 
C             depending on the values of INFO(20).
C                   
C  INFOBI(*) -- This is an integer array used to communicate details
C               of how the adjoint solution is to be carried out.
C               It must be dimensioned at least 10. You must respond to 
C               all of the following items, which are arranged as questions.  
C               The simplest use of DDASPKADJOINT corresponds to setting 
C               all entries of INFOBI to 0.
C
C   INFOBI(1) - The error tolerance for the adjoint variables may be input by
C               the user or given by default. The default tolerances for
C               the adjoint variables are 2*ATOL(i) and 2*RTOL(i) where ATOL 
C               and RTOL are the tolerances for the corresponding state 
C               variables.
C
C            **** Do you want to use the default RTOLB and ATOLB?
C                 Yes --- set INFOBI(1) = 0;
C                 No  --- set INFOBI(1) = 1 and input RTOLB and ATOLB.
C                         RTOLB and ATOLB can be scalars if INFO(2) = 0,
C                         and can be a vector of length NEQ+INFO(22) if 
C                         INFO(2) = 1. 
C                 
C   INFOBI(2) - Options for initialization of the adjoint system.
C               The initial conditions for the adjoint DAEs are very easy
C               to obtain for an index-0 DAE. However, for index-1 and
C               index-2 systems, the initial conditions need special 
C               treatment (see reference[1] at the end of this
C               documentation). Although it is less trouble for you
C               to have DASPKADJOINT compute the consistent initial
C               conditions for the adjoint DAEs, the initialization inside 
C               DASPKADJOINT is for the general case without any 
C               assumptions about the structure of the original DAEs, and 
C               it might not be efficient for some special problems.
C
C          **** Do you want DASPKADJOINT to do the initialization for the 
C               adjoint DAEs ?
C               No  -- set INFOBI(2) = 3 and provide routine ADINIT to 
C                      specify the consistent initial conditions.
C                     (see the description of the item ADINIT in the call list)
C               Yes -- There are three options.
C                   1. INFOBI(2) = 0 (default): 
C                      This option is for a standard ODE system y' = f(t, y).
C                   2. INFOBI(2) = 1: This option is for index-1 DAE system,
C                      or index-0 DAE (implicit ODE system). 
C                   3. INFOBI(2) = 2: This option is for Hessenberg index-2
c                      DAE systems, or mixed index-1 and index-2 DAE systems.
C                      If this option has been chosen, routines RES_ADY  
C                      and G_RES_ADY must be provided (see the descriptions of 
C                      RES_ADY and G_RES_ADY in the call list).
C                 
C      ***Note*** For INFOBI(2) = 1 or 2, IEOPT(*) MUST be provided by the 
C                 user. If INFO(16) > 0, IEOPT(*) is also required.
C                 (see the description of item IEOPT in this documentation).
C
C    INFOBI(3) - Options for obtaining the adjoint DAEs. The adjoint DAE system
C 
C                (F_y')^*ADYPRIME + ((F_y)^* - ((F_y')_t)^*) ADY = 0  (****)
C
C                is solved by DASPK3.0. The adjoint DAE system can be 
C                constructed either by automatic differentiation (AD) tools 
C                inside DASPK3.0, or for efficiency reasons by the user 
C                outside DASPK3.0.
C
C                If the original DAEs are linear, F_y and F_y' are
C                independent of Y and YPRIME. In this case, the forward 
C                information is not required in constructing the 
C                adjoint DAEs.
C                          
C          ****  Do you want to obtain the adjoint DAE system by automatic 
C                differentiation (AD)? 
C                No - set INFOBI(3) = 3 if the original DAE system is 
C                                          nonlinear, or
C                     set INFOBI(3) = 7 if the original DAE system is linear,
C                         and provide a user-supplied ADRES, 
C
C                     SUBROUTINE ADRES(
C                    *           T,ADY,ADYPRIME,CJ,DELTA,IRES,RPAR,
C                    *           IPAR,SENPAR,Y)
C 
C                     where DELTA is the output variable and the rest of the
C                     arguments are the input variables.
C
C               Yes - There are three options.
C                     1. If the system is a standard-form ODE y' = f(t,y), 
C                        set INFOBI(3) = 0 if the original ODE system is 
C                                             nonlinear or
C                        set INFOBI(3) = 4 if the original ODE system is linear,
C                        and provide the AD-generated routine 
C                        in argument RES_ADY (see the description of  the item  
C                                             RES_ADY in the call list).
C                     2. If the system is of general form F(t,y,y') = 0,
C                        and F_y' is a constant matrix, 
C                        set INFOBI(3) = 1 if the original DAE system is 
C                                          nonlinear or
C                        set INFOBI(3) = 5 if the original DAE system is linear,
C                        and provide RES_ADY and the other AD-generated 
C                        routine--RES_ADYPRIME--in argument ADRES.
C
C                        RES_ADYPRIME can be generated via a reverse mode 
C                        AD-tool or constructed by the user. If it is 
C                        generated by a reverse mode AD-tool, YPRIME is 
C                        the independent variable and DELTA is the dependent 
C                        variable. The form of RES_ADYPRIME must be
C
C                        SUBROUTINE RES_ADYPRIME(
C                      *            T, Y, YPRIME, ADYPRIME, CJ, DELTA, ADDELTA,
C                      *            IRES, RPAR, IPAR, SENPAR)
C
C                        where ADYPRIME is the output variable and the rest of
C                        the arguments are the input variables, ADDELTA is a 
C                        seed matrix.
C
C                        If RES_ADYPRIME is constructed by the user, it 
C                        should calculate the vector-matrix product and 
C                        return the result in ADYPRIME:
C
C                               ADYPRIME = ADDELTA*F_y'
C
C                        Example using Tapenade as AD tool:
C                        We want to perform automatic differentiation on the function RES
C                        (as described in the DAPSK3.1 docs), but we have to play 
C                        a game with filenames. If we ask Tapenade to create the code
C                        for RES differentiated with respect to SENPAR, it will create
C                        a file named RES_B. Then, later when we want Tapenade to
C                        create the code for RES differentiated with respect to Y, it
C                        will overwrite the file RES_B with the new code. The same 
C                        follows for differentiating RES backwards with respect to
C                        YPRIME. So, in order to keep all codes around, we employ a 
C                        work-around.  We make three copies of the same code - name one 
C                        RES the second RES_OFP, and the third RES_OFYP (and do the same 
C                        for any subroutines called by RES). Then, we differentiate RES 
C                        w.r.t. Y, RES_OFP w.r.t. SENPAR, and RESOFYP w.r.t. YPRIME. 
C                        So, if RES_OFYP resides in a file odes.f, then we will use 
C                        tapenade in gradient (reverse) mode, to produce ADDELTA*F_y'.
C                        We set "YPRIME" as the independent variable (-vars "YPRIME") and
C                        "DELTA" as the dependent variable (-outvars "DELTA"):
C
C                        java -mx256m -classpath \
C                          "C:\tapenade2.1.3\jars\tapenade.jar" \
C                          -Dali.tapenade_home="C:\TAPENA~1.3" \
C                          -Dali.browser="dummy" \
C                          topLevel.Differentiator \
C                          -backward -head RES_OFYP -vars "YPRIME" -outvars "DELTA" \
C                          "C:\DASPK3P1\examples\odes.f"
C                 
C                        The result will be a file res_ofyp_b.f containing the function
C
C                        SUBROUTINE RES_OFYP_B(t, y, yprime, cj, delta, deltab, 
C               +           ires, rpar, ipar, senpar, senparb)
C
C                        Even though DELTAB is the input, it is also treated as output.
C                        The implication is that if we pass ADDELTA into RES_OFYP_B, it 
C                        will be clobbered. Additionally, since RES_OFYP presumably 
C                        alters DELTA, we must assume RES_OFYP_B alters DELTA. This means 
C                        RES_ADYPRIME must be a wrapper for RES_OFYP_B that protects both 
C                        DELTA and ADDELTA.  For this particular example, we have a two-
C                        dimensional system, so we dimension DELTA and DELTAB to be 2.
C                        Additionally, we need to save the value of ADYPRIME, and
C                        add it to the result of RES_OFYP_B.
C
C                        SUBROUTINE RES_ADYPRIME(
C                      *            T, Y, YPRIME, ADYPRIME, CJ, DELTA, ADDELTA,
C                      *            IRES, RPAR, IPAR, SENPAR)
C                        IMPLICIT NONE
C                        DOUBLE PRECISION T, Y, YPRIME, DELTA, RPAR, ADDELTA
C                        DOUBLE PRECISION SENPAR,CJ
C                        INTEGER IPAR,IRES,II
C                        DIMENSION YPRIME(*),DELTA(*),Y(*),RPAR(*),ADDELTA(*)
C                        DIMENSION SENPAR(*)
C                        DIMENSION IPAR(*)
C                        DOUBLE PRECISION DELTAB(2), DELTA_COPY(2), ADYPRIME_SAV(2)
C
C                        DO II = 1, 2
C                          DELTAB(II) = ADDELTA(II)
C                          DELTA_COPY(II) = DELTA(II)
C                          ADYPRIME_SAV(II) = ADYPRIME(II)
C                        END DO
C                        CALL RESWEB_OFYP_B(T, Y, YPRIME, ADYPRIME, CJ, DELTA_COPY, 
C                       +                   DELTAB, IRES, RPAR, IPAR, SENPAR)
C                        DO II = 1, 2
C                          ADYPRIMEII) = ADYPRIME(II)+ADYPRIME_SAV(II)
C                        END DO
C                        RETURN
C                        END
C
C                     3. If the system is of general form F(t,y,y') and F_y' 
C                        depends on t and/or y, set INFOBI(3) = 2 if the 
C                        original DAE system is nonlinear or
C                        set INFOBI(3) = 6 if the original DAEs is linear. In
C                        addition to provide routines RES_ADY and RES_ADYPRIME
C                        as the case of INFOBI(3) = 1, the user Must provide 
C                        another routine G_RES_ADYP, which can be generated 
C                        by the combination of forward and adjoint mode of 
C                        an AD tool (see the descriptions of G_RES_ADYP in the 
C                        call list).
C               
C       ***Note*** If some sensitivity parameters appear in RES (specified 
C                  by INFO(22) > 0), the forward information is required in
C                  every time step in the sensitivity computation. Therefore,
C                  INFOBI(3) must be less than 4.
C
C    INFOBI(4) -- Used only when INFO(22) > 0.
C                 If some of the sensitivity parameters appear in RES, 
C                 the sensitivity computation requires the quadrature 
C                     \int_0^T (ADY*F_p) dt.
C                 The user can choose to calculate the quadrature outside of
C                 DASPK by the trapezoidal rule, or to treat the quadrature
C                 via quadrature variables inside DASPK3.0. The second option
C                 might take a little more computation time but is more accurate.
C
C           ****  Do you want to compute the quadrature of the parameter 
C                 sensitivities externally from DASPK?
C                   No  -- Set INFOBI(4) = 0. 
C                          The quadrature will be taken as a variable in DASPK.
C                   Yes -- Set INFOBI(4) = 1; then the quadrature will be 
C                          computed by trapzoidal rule.
C
C    INFOBI(5) -- Used only when the DAE system is index-2.
C                 An index-2 DAE requires special treatment in the Jacobian
C                 evaluation by ADIFOR.
C
C           ****  Is the system of index-2?
C                   No  -- Set INFOBI(5) = 0. 
C                   Yes -- Set INFOBI(5) = 1 and set up the IEOPT(*) array.
C
C    INFOBI(6) -- Used only if some of the sensitivity parameters appear in the
C                 initial conditions. The sensitivity computation requires  
C                 the value of ADY*F_y'*(Y_0)_p. Usually the initial values of 
C                 the sensitivities are very simple, and (Y_0)_p contains
C                 a few nonzero elements. Although DASPKADJOINT can compute
C                 ADY*F_y'*(Y_0)_p, it requires the user to input the
C                 dense form of (Y_0)_p, which might be huge if both NEQ and 
C                 INFO(19) are very large. 
C
C            **** Do you want DASPKADJOINT to calculate ADY*F_y'*(Y_0)_p?
C                    No  -- Set INFOBI(6) = 0. DASPKADJOINT will only 
C                           calculate ADY*F_y' and store it in QSEN.
C                           (see the description of QSEN in the OUTPUT part
C                           in this documentation). 
C                    Yes -- Set INFOBI(6) = 1, and provide the consistent
C                           initial values of (Y_0)_p in the Y(*) array and
C                           dimension Y at least NEQ*(1+INFO(19)).
C    
C    INFOBI(7) -- This is the number of time steps required for the forward 
C                 problem. If you do not know in advance, set INFOBI(7) = 0.
C                 This number can usually be obtained by solving the forward 
C                 problem first. Otherwise, this number must be set to no less
C                 than the actual number of time steps of the forward problem.
C                 
C                 The purpose of setting this number is to reduce two runs of 
C                 forward problem to one. If the size of buffer (NBUF in 
C                 the following) is larger than 2*INFOBI(7)*NEQ, then the second
C                 run can be avoided and the computation can be faster.
C
C    INFOBI(8) -- This is an option for automatic differentiation.
C                 The user can use either TAMC, ADIFOR3.0, Tapenade, or
C                 other automatic differentiation tool to generate
C                 the derivative code. If TAMC is used, the user need
C                 to write a wrapper with the required argument list
C                 to call the TAMC-generated routine. 
C
C            **** Do you use ADIFOR3.0?
C                    Yes -- Set INFOBI(8) = 1, and link the ADIFOR
C                           library when compiling the code.
C                    No  -- Set INFOBI(8) = 0, and provided a wrapper
C                           with the required argument list.  
C
C  RTOLB,ATOLB -- These quantities represent absolute and relative
C                 error tolerances (on local error) which you provide
C                 to indicate how accurately you wish the adjoint solution to
C                 be computed. When you specify INFOBI(1)=1, these two
c                 quantities must be provided.
C
C  NBUF -- This is the size of the buffers (double precision) used to 
C          communicate between forward integration and backward
C          adjoint integration. NBUF should be at least 
C          13*NEQ. Setting NBUF to be N*13*NEQ in your calling program is the
C          most efficient way to use the buffers, where N is the number of 
C          checkpoints to be allowed to stay in the memory. 
C
C          NBUF determines the number of checkpoints which should be set 
C          during the forward integration. At each checkpoint, DASPKADJOINT
C          might need to communicate with the DISK file, which will slow 
C          down the computation. There is a trade-off between efficiency
C          and memory requirements.
C
C  ADJAC -- This is the name of a routine which you may provide (optionally) 
C           for calculating Jacobian partial derivative data involved in 
C           solving linear systems within DASPKADJOINT. This ADJAC is for the 
C           adjoint DAE. The role of ADJAC (and its call sequence) depends on 
C           whether a direct (INFO(12) = 0) or Krylov (INFO(12) = 1) method 
C           is selected.
C 
C         **** INFO(12) = 0 (direct methods):
C           If INFO(5) = 1 is chosen, you must supply an ADJAC routine
C           to compute the matrix A = dG/dY + CJ*dG/dYPRIME for the adjoint 
C           system. It must have the form
C
C           SUBROUTINE ADJAC (T, Y, YPRIME, PD, CJ, RPAR, IPAR, SENPAR, IJAC)
C 
C           where Y and YPRIME are the state variables and their time 
C           derivatives at time T. PD stores the matrix A, which is the 
C           transpose of the Jacobian matrix for the forward integration 
C           except in the case of index-2 DAEs. 
C
C           For index-2 DAEs, the Jacobian is scaled by CJ in the rows
C           corresponding to the index-2 constraints, and scaled by 1/CJ in the
C           column corresponding to the index-2 variables. We assume that the
C           index-2 variables never appear in the index-2 constraints.
C
C           If INFO(5) > 1 is chosen (ADIFOR option), ADJAC is the same as JAC.
C  
C           If the mass matrix F_y' depends either on y or t, and the ADIFOR
C           option is chosen, you MUST add another variable z = y' for all of 
C           the differential variables and change the mass matrix to be 
C           constant. 
C
C         **** INFO(12) = 1 (Krylov methods):
C           In this case, ADJAC has a different argument list from JAC in 
C           DDASPK3.0. The form is 
C                
C              SUBROUTINE ADJAC (
C            *            ADRES, IRES, NY, X, ADY, ADYPRIME, WT, SAVR, DELTA,
C            *            H,CJ, WP, IWP, IERPJ, RPAR, IPAR,SENPAR, ADI,
C            *            ISENFO,SENWRK,ISENWK,RES_ADY)
C
C           where the extra arguments (ADRES, ISENFO,SENWRK,ISENWK,
C           RES_ADY, ADI) are for the residual calculations for the 
C           adjoint equations. If you need to call the RES routine in ADJAC, 
C           the RES routine should be replaced by
C                 
C                DDRESAD(X, ADY, ADYPRIME, CJ, DELTA, IRES, RPAR, IPAR, SENPAR,
C            *         ADRES, NEQ, ISENFO, SENWRK, ISENWK, RES_ADY, ADI)
C               
C           which is a subroutine inside DASPK3.0. You only need to know that
C           DDRESAD does the same thing for the backward integration as RES for
C           the forward integration. You do not need to know the details of 
C           the routine DDRESAD 
C
C  IEOPT(*) -- This array of length 2*NEQ must be provided by the user
C              if INFOBI(2) > 0, or INFO(16) = 1, or INFOBI(5) = 1 is chosen.
C
C              The user MUST specify which equations in the original DAEs 
C              are differential equations, which are index-1 algebraic 
C              constraints and which are index-2 algebraic constraints 
C              by setting
C                 
C               IEOPT(i) = -1, if the ith equation is an index-1 constraint;
C               IEOPT(i) = -2, if the ith equation is an index-2 constraint;
C               IEOPT(i) =  1, if the ith equation is a differential equation.
C
C              and MUST set IEOPT(NEQ+1:2*NEQ) as follows if INFOBI(2) > 0 or 
C              INFOBI(5) = 1 is chosen:
C
C               IEOPT(NEQ+I) = 1, 2 or 3 if Y(I) is a differential variable;
C               IEOPT(NEQ+I) = -1,  if Y(I) is an index-1 algebraic variable;
C               IEOPT(NEQ+I) = -2,  if Y(I) is an index-2 algebraic variable;
C
C              If you have chosen to initialize the forward integration 
C              (INFO(11) > 0), you can just copy the NEQ elements in the 
C              IWORK array starting from LID+1 to the IEOPT array starting 
C              from NEQ+1, where
C                    LID = 40     if INFO(10) = 0 or 2 and 
C                    LID = 40+NEQ if INFO(10) = 1 or 3.
C
C
C  RES_ADP:EXT -- This is the name of a routine which you may supply 
C                 (optionally) that relates to the evaluation of the
C                 adjoint DAEs. If you set INFO(22) = 0 (no parameter
C                 appears in RES), then RES_ADP can be absent (or perhaps a 
C                 dummy routine to satisfy the loader). Othewise if you
C                 set INFO(22) > 0, then you must supply a RES_ADP to compute
C                 the vector-matrix product vF_p, where v is a vector and F_p
C                 is the Jacobian of F with respect to the parameters p. 
C                 RES_ADP must have the form 
C
c                 SUBROUTINE RES_ADP(
C              *             T, Y, YPRIME, CJ, DELTA, ADDELTA, IRES, RPAR, 
C              *             IPAR, SENPAR, ADSENPAR)
C 
C                 where ADSENPAR = ADDELTA*F_p, Y and YPRIME are the forward 
C                 information at time T.
C
C                 If this routine is generated by a reverse mode AD tool, 
C                 SENPAR is the independent variable and DELTA is the 
C                 dependent variable. The dimension of SENPAR should be 
C                 equal to INFO(22).
C
C                 Example using Tapenade as AD tool:
C                 We want to perform automatic differentiation on the function RES
C                 (as described in the DAPSK3.1 docs), but we have to play 
C                 a game with filenames. If we ask Tapenade to create the code
C                 for RES differentiated with respect to SENPAR, it will create
C                 a file named RES_B. Then, later when we want Tapenade to
C                 create the code for RES differentiated with respect to Y, it
C                 will overwrite the file RES_B with the new code. The same 
C                 follows for differentiating RES backwards with respect to
C                 YPRIME. So, in order to keep all codes around, we employ a 
C                 work-around.  We make three copies of the same code - name one 
C                 RES the second RES_OFP, and the third RES_OFYP (and do the same 
C                 for any subroutines called by RES). Then, we differentiate RES 
C                 w.r.t. Y, RES_OFP w.r.t. SENPAR, and RESOFYP w.r.t. YPRIME. 
C                 So, if RES_OFP resides in a file odes.f, then we will use 
C                 tapenade in gradient (reverse) mode, to produce ADDELTA*F_p.
C                 We set "SENPAR" as the independent variable (-vars "SENPAR") and
C                 "DELTA" as the dependent variable (-outvars "DELTA"):
C
C                 java -mx256m -classpath \
C                   "C:\tapenade2.1.3\jars\tapenade.jar" \
C                   -Dali.tapenade_home="C:\TAPENA~1.3" \
C                   -Dali.browser="dummy" \
C                   topLevel.Differentiator \
C                   -backward -head RES_OFP -vars "SENPAR" -outvars "DELTA" \
C                   "C:\DASPK3P1\examples\odes.f"
C                 
C                 The result will be a file res_ofp_b.f containing the function
C
C                 SUBROUTINE RES_OFP_B(t, u, uprime, cj, delta, deltab, ires, rpar, 
C     +                     ipar, senpar, senparb)
C
C                 Even though DELTAB is the input, it is also treated as output.
C                 The implication is that if we pass ADDELTA into RES_OFP_B, it will
C                 be clobbered. Additionally, since RES_OFP presumably alters DELTA, 
C                 we must assume RES_OFP_B alters DELTA. This means RES_ADP must be a 
C                 wrapper for RES_OFP_B that protects both DELTA and ADDELTA.  For this 
C                 particular example, we have a two-dimensional system, so we dimension
C                 DELTA and DELTAB to be 2. Additionally, we need to save the value of 
C                 ADYSENPAR, and add it to the result of RES_OFP_B.

C
C                 SUBROUTINE RES_ADP(
C     *             T, Y, YPRIME, CJ, DELTA, ADDELTA, IRES, RP, 
C     *             IP, SENPAR, ADSENPAR)
C                 IMPLICIT NONE
C                 DOUBLE PRECISION T, Y, YPRIME, DELTA, RP, ADDELTA
C                 DOUBLE PRECISION SENPAR,CJ,ADSENPAR
C                 INTEGER IP,IRES,II
C                 DIMENSION YPRIME(*),DELTA(*),Y(*),RP(*),ADDELTA(*)
C                 DIMENSION SENPAR(*), ADSENPAR(*)
C                 DIMENSION IP(*)
C                 DOUBLE PRECISION DELTAB(2), DELTA_COPY(2), ADSENPAR_SAV(2)
C
C                 DO II = 1, 2
C                   DELTAB(II) = ADDELTA(II)
C                   DELTA_COPY(II) = DELTA(II)
C                   ADSENPAR_SAV(II) = ADSENPAR(II)
C                 END DO
C                 CALL RES_OFP_B(t, y, yprime, cj, delta_copy, deltab,
C     +               ires, rp, ip, senpar, adsenpar)
C                 DO II = 1, 2
C                   ADSENPAR(II) = ADSENPAR(II)+ADSENPAR_SAV(II)
C                 END DO
C                 return
C                 end
C
C  RES_ADY:EXT -- This is the name of a routine which you may supply 
C                 (optionally) that relates to the evaluation of the 
C                 adjoint DAEs. If you set INFOBI(3) = 3 or 7 (user-input
C                 ADRES), then RES_ADY can be absent. Otherwise if you 
C                 choose INFOBI(3) = 0,1,2,4,5,6, then you must supply a 
C                 RES_ADY to compute the vector-matrix product vF_y,
C                 where v is a vector and F_y is the Jacobian of F with 
C                 respect to y. RES_ADY must have the form
C
C                 SUBROUTINE RES_ADY(T, Y, ADY, YPRIME, CJ, DELTA, ADDELTA,
C              *        IRES, RPAR, IPAR, SENPAR)
C
C                 where ADDELTA is a vector (input), ADY = ADDELTA*F_y.
C
C                 To generate this routine, use a reverse mode AD tool to
C                 differentiate RES with Y as the independent variable and
C                 DELTA as the dependent variable.
C                 
C                 Example using Tapenade as AD tool:
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
C
C  NEQAD:IN      -- An integer specifying the size of the adjoint variables.
C                   It must be set according to the following rules:
C
C                   IF INFOBI(4) .EQ. 0 .AND. INFO(22) .GT. 0 THEN
C                      NEQAD = NQ * (NEQ + INFO(22))
C                   ELSE
C                      NEQAD = NQ * NEQ
C                   END IF
C
C  ADY(NEQAD):OUT -- Contains the adjoint variable. If INFO(3) = 0, then ADY may
C                   be of length 0 (but NEQAD must still be specified correctly).
C
C  ADYP(NEQAD):OUT -- Contains the derivative of the adjoint variable. If INFO(3) = 0, 
C                   then ADY may be of length 0 (but NEQAD must still be specified 
C                   correctly).
C
C  G_RES_ADY:EXT -- This is the name of the AD-generated routine which 
C                   you may supply (optionally) that relates to the 
C                   initialization of an index-2 adjoint system.  
C                   If you choose INFOBI(5) = 0 (not index-2), then 
C                   G_RES_ADY can be absent (or perhaps a dummy routine to 
C                   satisfy the loader). Otherwise if you choose INFOBI(5)=1, 
C                   then you must supply a G_RES_ADY to compute the time 
C                   derivative of the vector-matrix product vF_y, which is
C                     u*F_y + v*(dF_y/dt + dF_y/dy * y')
C                   where u and v can be of any vector, dF_y/dy is a tensor.
C
C                   G_RES_ADY must have the form:
C 
C                   SUBROUTINE G_RES_ADY(
C                 *            T, G_t, Y, G_y, ADY, G_ADY, YPRIME, CJ, 
C                 *            DELTA, ADDELTA, G_ADDELTA, IRES, RPAR, IPAR, 
C                 *            SENPAR)
C
C                   where ADY and G_ADY are the output variables and the
C                   rest of the arguments are the input variables. 
C                   G_RES_ADY calculates
C
C                    ADY = ADDELTA * F_y,
C                    G_ADY = G_ADDELTA*F_y + ADDELTA*(dF_y/dt * G_t +
C                                                    + dF_y/dy * G_y)
C                          
C                   where G_t, G_y and G_ADDELTA are seed matrices for 
C                   the forward mode and ADDELTA is a seed matrix for 
C                   the reverse mode if an AD-tool is used.
C
C                   G_RES_ADY can be generated by a combination of reverse 
C                   mode and forward mode AD tools. First, RES_ADY is 
C                   generated by a reverse mode AD tool. Then G_RES_ADY
C                   is generated via ADIFOR for RES_ADY with T, Y and 
C                   ADDELTA as independent variables, and ADY as dependent 
C                   variables. 
C
C  G_RES_ADYP:EXT -- This is the name of the AD-generated routine which 
C                   you may supply (optionally) that relates to the 
C                   initialization of an adjoint system. If INFOBI(3) is not
C                   equal to 2 or 6, then G_RES_ADYP can be absent (or perhaps 
C                   a dummy routine to satisfy the loader). Otherwise if you 
C                   choose INFOBI(3) = 2 or 6, then you must supply a 
C                   G_RES_ADYP to compute the time derivative of the 
C                   vector-matrix product vF_y', which is
C                     u*F_y' + v*(dF_y'/dt + dF_y'/dy * y')
C                   where u and v can be of any vector, dF_y'/dy is a tensor.
C
C                   G_RES_ADYP must have the form:
C  
C                   SUBROUTINE G_RES_ADYP(
C                 *            T, G_t Y, G_y, YPRIME, ADYP, G_ADYP, CJ, 
C                 *            DELTA, ADDELTA, G_ADDELTA, IRES, RPAR, IPAR, 
C                 *            SENPAR)
C
C                    where G_t, G_y and G_ADDELTA are seed matrices for 
C                    the forward mode and ADDELTA is a seed matrix for 
C                    the reverse mode.
C
C                    G_RES_ADYP can be generated by a combination of  
C                    reverse mode and forward mode AD tools. 
C                    First RES_ADYP is generated by a reverse 
C                    mode AD tool. Then G_RES_ADYP is generated 
C                    via ADIFOR for RES_ADYPRIME with T, Y and 
C                    ADDELTA as independent variables, and ADYP 
C                    as dependent variables. 
C
C                    If G_RES_ADYP is constructed by the user, it should 
C                    compute the time derivative of 
C                             ADYP=ADDELTA*F_y'.
C                    which yields
C
C                        G_ADYP = G_ADDELTA*F_y' 
C                               + ADDELTA*(dF_y'/dt*G_t + dF_y'/dy(*)G_y)
C
C                    where dF_y'/dy is a tensor.
C              
C  ADINIT:EXT -- This is the name of the routine which you may supply
C                (optionally) that relates to the initialization of the 
C                adjoint system. If you choose INFOBI(2) < 3, ADINIT can
C                be absent. Otherwise if you choose INFOBI(2) = 3, then you
C                must supply an ADINIT to initialize the adjoint system
C                at time T = TFINAL. The purpose of ADINIT is to
C                calculate the consistent initial conditions for the
C                adjoint DAEs at T = TFINAL. In addition to making the 
C                adjoint DAEs consistent, the initial conditions
C                for the adjoint variables should satisfy the following
C                boundary conditions:
C
C                  ADY*F_y'*Y_p = 0 at t = TFINAL
C
C                In reference [1], we have discussed in detail how to
C                obtain such consistent initial conditions and implemented 
C                it in the DASPKADJOINT program. However, our initialization 
C                is for the general case without any assumptions about the 
C                structure of the original DAEs, and it might not be 
C                efficient for some special problems.
C 
C                If you wish to compute the consistent initial 
C                conditions by a user-input routine, you MUST set
C                INFOBI(2) = 3 and supply ADINIT, which has the form
C
C                 SUBROUTINE ADINIT(T, NEQAD, NQ, Y, YPRIME, ADY, ADYPRIME,
C              *             QSEN, RPAR, IPAR, SENPAR)                    
C
C                where NEQAD=(NEQ+NSENPAR)*NQ is the total number of equations in 
C                the adjoint systems, NQ is the total number of derived 
C                functions, and QSEN stores the values of the derived functions 
C                and their gradients with respect to Y and SENPAR (it is the
C                output of QRES). ADINIT must set ADY and ADYP. (Note that QSEN is
C                stacked like the transpose of ADY).
C
C                ADY is stacked like the return value ADY.  For an ODE system,
C                this means ADINIT returns
C                    dq1(TFINAL)/dy1(TFINAL)
C                      ...
C                    dq1(TFINAL/dyNEQ(TFINAL)
C                    dq1(TFINAL)/dp1
C                    dq1(TFINAL)/dpNSENPAR
C                    dq2(TFINAL)/dy1(TFINAL)
C                      ...
C                      ...
C                    dqNQ(TFINAL)/dy1(TFINAL)
C                    dqNQ(TFINAL/dyNEQ(TFINAL)
C                    dqNQ(TFINAL)/dp1
C                    dqNQ(TFINAL)/dpNSENPAR
C 
C               where the q's are the derived equations, the y's are the 
C               state variables, and the p's are the parameters, NQ is the
C               number of derived equations, NEQ is the number of 
C               state variables, and NSENPAR is the number of parameters
C               to which we are finding the sensitivity. Each dq(TFINAL)/dp should
C               be initialized to zero.
C 
C               For an ODE system, ADYP is set to all zeroes.
C
C  K_RES:EXT -- This is the name of a routine used by DASPKADJOINT for the 
C               forward integration to evaluate the matrix-vector product in 
C               the Krylov iterative method. You must supply this routine 
C               if you select INFO(5) > 2 and INFO(12) = 1.
C
C               See the documentation of DDASPK3.0 for more details.
C             
C  T_RES:EXT -- This is the name of a routine used by DASPKADJOINT for the 
C               forward integration to initialize an index-2 DAE. You must 
C               supply it if you choose INFO(11) = 4 or 5 and INFO(20) > 2. 
C               See the documentation of DDASPK3.0 for more details.
C       
C  G_RES_ADP:EXT -- This is the name of the AD-generated routine which 
C                   you may supply (optionally) that relates to the 
C                   sensitivity evaluation of an index-2 system. 
C                   If you choose INFOBI(5) = 0 (not index-2), then 
C                   G_RES_ADY can be absent (or perhaps a dummy routine to 
C                   satisfy the loader). Otherwise if you choose INFOBI(5) = 1, 
C                   then you must supply a G_RES_ADP to compute the time 
C                   derivative of the vector-matrix product vF_p, which is
C                     u*F_p + v*(dF_p/dt + dF_p/dy * y')
C                   where u and v can be of any vector, dF_p/dy is a tensor.
C
C                   G_RES_ADP must have the form:
C 
C                   SUBROUTINE G_RES_ADP(
C                 *            T, G_t Y, G_y, YPRIME, DELTA, ADDELTA, 
C                 *            G_ADDELTA, CJ, IRES, RPAR, IPAR, 
C                 *            SENPAR, ADP, G_ADP)
C
C                   where ADP and G_ADP are the output variables and the
C                   rest of the arguments are the input variables. 
C                   G_RES_ADP calculates
C
C                    ADP = ADDELTA * F_p,
C                    G_ADP = G_ADDELTA*F_p + ADDELTA*(dF_p/dt * G_t +
C                                                    + dF_p/dy * G_y)
C                          
C                   where G_t, G_y and G_ADDELTA are seed matrices for 
C                   the forward mode and ADDELTA is a seed matrix for 
C                   the reverse mode if an AD-tool is used.
C
C                   G_RES_ADP can be generated by a combination of reverse 
C                   mode and forward mode AD tools. First, RES_ADP is 
C                   generated by a reverse mode AD tool. Then G_RES_ADP
C                   is generated via ADIFOR for RES_ADP with T, Y and 
C                   ADDELTA as independent variables, and ADP as dependent 
C                   variables. 
C
C       ***Note**   The purpose of this routine is for the evaluation of 
C                   vector-matrix product v{\dot h_p}, where h is the 
C                   index-2 constraints. Therefore, if index-2 constraints
C                   do not have any sensitivity parameters, the user can
C                   supply a simple routine which only assigns G_ADP to be 
C                   zero.
C                   
C                   If ADIFOR3.0 (or later version) is used, the script file
C                   "spec.ad3" for generating RES_ADY (or RES_ADYP, RES_ADP)
C                   can be written as:
C
C                    Function
C                      JtM
C                    MaxColsM
C                      1
C                    Independent
C                      y (or yprime for ADYP, senpar for ADP)
C                    Dependent
C                      Delta
C                    Top
C                      res
C                    ExceptionDefaultMode
C                      Performance
C                    Cache
C                      JtM_Cache
C                    OutputDir
C                      JtM_Output
C
C                  The script file "spec.ad3" for generating G_RES_ADY 
C                  (or G_RES_ADYP, G_RES_ADP)
C                  by forward mode can be written as 
C                    
C                   Function
C                     JM
C                   MaxColsM
C                     1
C                   Independent
C                     t
C                     y
C                     ady (or adyp, adp)
C                   Dependent
C                     addelta
C                   Top
C                     res_ady (or res_adyp, res_adp)
C                   ScalarDerivatives
C                     true
C                   ExceptionDefaultMode
C                     Performance
C                   Cache
C                     JM_Cache
C                   OutputDir
C                     JM_Output  
C                   
C------OUTPUT - AFTER ANY RETURN FROM DDASPKADJOINT----------------------
C
C  The principal aim of the code is to return the sensitivity of the
C  derived function with respect to the parameters p at T = TFINAL. Unlike DASPK, 
C  you cannot compute sensible sensitivities at intermediate times between T0 
C  and TFINAL. 
C  However, if INFO(3)=1,2 then the full system variables (Y and YPRIME) and the 
C  adjoint variables (ADY and ADYP) and T will contain intermediate values. 
C  They are the only variables that contain return values and no other variable 
C  should be saved or examined (QSEN will not be set).
C  As with DASPK, the IDID parameter should be checked aftern any return. 
C
C    Y(*) -- if INFO(3)=1,2 then Y and YPRIME contain the values of the original
C            system's state variables at time T.
C            Otherwise, it contains the computed solution approximation at TFINAL.
C
C    YPRIME(*) -- if INFO(3)=1,2 then Y and YPRIME contain the values of the original
C            system's state variables at time T.
C            Otherwise, it contains the computed solution approximation at TFINAL.
C
C    QSEN(*) -- contains the values of the derived functions and the sensitivity
C             information.
C
C             QSEN(1:NQ) contains the values of the derived functions.
C
C             If INFOBI(6) = 0 is chosen and there are some parameters in 
C             the initial conditions, QSEN(NQ+1:NQ*(NEQ+1)) contains the 
C             values of ADY*F_y'. Otherwise QSEN(NQ+1:NQ*(NEQ+1)) contains 
C             the values of the adjoint variables ADY at t=t_0.
C             
C             If INFOBI(6) = 1, QSEN contains the sensitivities of the 
C             derived functions g with respect to p, starting from 
C             NQ*(NEQ+1)+1,. The sensitivity of the ith function with 
C             respect to the jth parameter in p is stored at 
C
C                  QSEN(NQ*(NEQ+1) + (j-1)*NQ + i).
C
C             For those parameters that appear in the initial conditions,
C             if INFOBI(6) = 0 is chosen, the user needs to calculate
C             the sensitivities by
C                          QSEN(NQ+1:NQ*(NEQ+1)) * (Y_0)_p
C             after calling DASPKADJOINT.
C
C    IDID -- reports what the code did. Please refer to the documentation
C          of DASPK3.0 for more details.
C
C
C    RWORK, IWORK -- contain information which is usually of no interest
C                   to the user.  However you may be interested in the 
C                   performance data for the backward integration. 
C
C                   Starting from IWORK(39)+1, the IWORK(*) array contains
C                   an integer array to support continutation in backward mode,
C                   and then the whole integer work array for the backward integration.
C                   If checkpoints are used, IWORK is organized as
C                        IWORKF
C                        ICONT
C                        IWORKB
C                        ISENWK (Sensitivity analysis integer work vector)
C                   where IWORKF is at least 40 entries long.
C                   If checkpoints are not used, IWORK is organized as
C                        IWORKF/B
C                        ISENWK
C                        ICONT
C                   IWORK(39) contains the index of the last entry in IWORKF
C                   That means ICONT starts at IWORK(39)+1.
C                   The first entry in ICONT contains the index of the last entry of ICONT
C                   That means IWORKB starts at IWORK(IWORK(39)+1)+1
C                   To retrieve information from the backward mode IWORK array (IWORKB), you
C                     must use the following code to access various statistics:
C                     LICONT = IWORK(39)
C                     LIWORKB= IWORK(LICONT+1)
C                     IWORK(LIWORKB+11) contains the number of steps,
C                     IWORK(LIWORKB+12) contains the number of residual
C                                         evaluations.
C                     IWORK(LIWORKB+13) contains the number of Jacobian
C                                         evaluations,
C                     To compute the length of IWORK actually required,
C                       TOTAL = IWORK(17) ! length of first section of IWORK (IWORKF)
C                       TOTAL = TOTAL + 90 ! length of ICONT
C                       IF (LIWORKB .NE. 1) THEN
C                         TOTAL = TOTAL + IWORK(LIWORKB+17) ! len of IWORKB
C                       END IF
C
C                   Some of the information in IWORK contains statistics about 
C                   the real work array, RWORK.
C                   If checkpointing is used, RWORK is organized as
C                        RWORKF
C                        RWORKB
C                        unused space
C                        BUFFER (begins at LRWEND+1 and ends at LRW)
C                   If checkpointing is not used, RWORK is organized as
C                        RWORKF/B
C                        unused space
C                        BUFFER (begins at LRWEND+1 and ents at LRW)
C                   where LRWEND is stored at IWORK(LICONT+18)
C                     To compute the length of RWORK actually required,
C                        TOTAL = IWORK(18) ! length of first chunk of RWORK (RWORKF)
C                        LRWORKB = IWORK(40)
C                        IF (LRWORKB.NE.1) THEN
C                          TOTAL = TOTAL + IWORK(LIWORKB+18) ! length of RWORKB
C                        END IF
C                        TOTAL = TOTAL + LRW - IWORK(LICONT+18) + 1
C                   etc.
C                   (you can refer to the documentation of DASPK3.0
C                    for more performance data).
C
C
C                  Regarding RWORK:
C                   Starting from IWORK(40)+1, RWORK(*) array contains
C                   the whole real work space for the backward integration.
C
C    ADY, ADYP -- if INFO(3)=1,2 then ADY and AYP contain the values of the adjoint
C                  variables at time T.
C                  ADY is stacked like the transpose of QSEN.
C                  For example, in the case where ADY is initialized as Q_Y,
C                  DASPK_Adjoint calls QRES, which fills in QSEN. QSEN is then
C                  stored in the ADY and ADYP sections of the work vector with the 
C                  following code:
C              c...    default initialization: ADY = Q_Y, ADYP = 0.0 
C              c...    for ODE y' = f(t, y)
C                      do iq = 1, nq
C                        do j = 1, neq
C                          rwork(lady+(iq-1)*neqpns + j - 1) = QSEN(j*nq + iq)
C                      end do
C
C              c...    quadrature part
C                        do j = neq+1, neqpns
C                          rwork(lady + (iq-1)*neqpns + j-1) = 0.0d0
C                        end do
C                      end do
C                      do i = 1, neqad
C                        rwork(ladyp+i-1) = 0.0d0
C                      end do
C                  where NQ is the number of derived functions Q, NEQ is the number
C                  of equations in Y, and NEQPNS set according to the rules:
C                       IF (INFOB(4) .EQ. 0 .AND. INFO(22) .GT. 0) THEN
C                          NEQPNS = NEQ + INFO(22)
C                       ELSE
C                          NEQPNS = NEQ
C                       END IF
C
C-----------------------------------------------------------------------
C
C***REFERENCES
C  1. Yang Cao, Shengtai Li, Linda Petzold, and Radu Serban, 
C     Adjoint Sensitivity Analysis for Differential-Algebraic Equations:
C     The Adjoint DAE System and its Numerical Solution, 
C     SIAM J. Sci. Comput. 24:3 (2003), 1076-1089.
C  2. Shengtai Li and Linda R. Petzold, Software and Algorithms for 
C     Sensitivity Analysis of Large-Scale Differential Algebraic Systems,
C     J. Comput. Appl. Math. 125 (2000) 131-145.
C
C=========================================================================
