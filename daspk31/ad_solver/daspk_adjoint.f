      SUBROUTINE DDASPKADJOINT(
     *     RES, NEQ, T, Y, YPRIME, TOUT, TFINAL, INFOI, RTOL, ATOL,
     *     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL, SENPAR,
     *     ADRES, NQ, QRES, QSEN, INFOBI, RTOLB, ATOLB, NBUF, 
     *     ADJAC, IEOPT, RES_ADP, RES_ADY,
     *     NEQAD, ADY, ADYP,
     *     G_RES_ADY, G_RES_ADYP,
     *     ADINIT, K_RES, T_RES, G_RES_ADP)
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
      IMPLICIT NONE
      
C     Parameters of the subroutine
      INTEGER   INFOI, IWORK, IEOPT, INFOBI, IPAR
      INTEGER   LRW, LIW, IDID, NEQ, NQ, NBUF
      INTEGER   NEQAD
      DOUBLE PRECISION T, Y, YPRIME, TOUT, TFINAL, RWORK
      DOUBLE PRECISION RTOL, ATOL, RPAR, SENPAR, QSEN
      DOUBLE PRECISION RTOLB, ATOLB, ADY, ADYP
      DIMENSION INFOI(30),IWORK(LIW), IEOPT(*), INFOBI(*)
      DIMENSION Y(*),YPRIME(*)
      DIMENSION RWORK(LRW)
      DIMENSION RTOL(*),ATOL(*), rtolb(*), atolb(*)
      DIMENSION RPAR(*),IPAR(*),SENPAR(*), QSEN(*)
      DIMENSION ADY(*), ADYP(*)
      EXTERNAL  RES, JAC, PSOL, RES_ADY, ADRES, ADJAC, QRES
      EXTERNAL  RES_ADP, K_RES, T_RES, G_RES_ADY, G_RES_ADYP, ADINIT
      EXTERNAL  G_RES_ADP
      
C     Transient local variables
      INTEGER   II, I, J, K, ITEMP, IQ
      INTEGER   I1, I2, I3, I4, I5
      INTEGER   INDEX1, IRES, NPHI
      INTEGER   LF_P2, LF_P3, LF_P4, IREC1, IREC2
      INTEGER   INFO24
      Character MSG*80
      DOUBLE PRECISION TAB, TAE, ADYA, ADYPA
      DIMENSION ADYA(NEQAD), ADYPA(NEQAD)
      INTEGER   LTLOC, LYLOC, LYPLOC ! These are used in 2 loops - they need to be
                                     ! passed to FORWARD_INTEGRATE with values in them.

C     Variables that are fairly constant and can be recomputed
C     on every call      
      INTEGER   ISIZEI, ISIZED
      INTEGER   ISTOL, ICOPT, MEVA, IQUADV, NSENPAR, NPARA
      INTEGER   INDEX2
      INTEGER   NEQPNS, LIWKS, LSEN
      DOUBLE PRECISION TEND, TGOALB, TBRET
      
C     Variables that need to be computed once, then stored 
C     and retrieved on continuation calls
      DOUBLE PRECISION TSTART, INITSTEP, CJ
      INTEGER   LISE, LSE, onePass, NCHK, NCHKSTP
      INTEGER   LICONT, LRCONT, LIWORKB
      INTEGER   LYF, LYPF, LPHIS, LRWKS, LTCHK, LF_P
      INTEGER   LADY, LADYP, LY, LT, LYP, LRWEND, LRWORKB, LRWF
      INTEGER   TDIR, LPHI, LIWB, LRWB, LADI, LTMP, IFWRITE
      
C     Variables that keep track of where we are in the overall
C     computation
      INTEGER   chkPos, ICHK
      INTEGER   INFOB, INFO
      DIMENSION INFOB(30), INFO(30)
      
C     Variables used only the first time this is called.
      INTEGER   NBUFLOC, NBUFH, LGRWK, LGMX, LGRWKB, lg_x
      INTEGER   LYSAVE, LTPRE
      DOUBLE PRECISION g_t
            
C     Timing variables used in the backward integration
      DOUBLE PRECISION TOUTF, TOUTB, TB, TF, TNEXT
      
C     Constants
      INTEGER LML,LMU,LMTYPE,LIWM,LMXORD,LJCALC,LPHASE,LK,LKOLD,
     *        LNS,LNSTL,LNST,LNRE,LNJE,LETF,LNCFN,
     *        LNCFL, LNIW, LNRW, LNNI, LNLI, LNPS,
     *        LNSE, LMITER, LMAXL, LKMP, LNRMAX, LLNWP,
     *        LLNIWP, LLOCWP, LLCIWP, LKPRIN, LMXNIT,
     *        LMXNJ, LMXNH, LLSOFF, LNPD, LNY, 
     *        LNLIS, LICNS
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
      INTEGER LTSTOP, LHMAX, LH, LTN, 
     *   LHOLD, LS, LROUND, LEPLI, LSQRN, LRSQRN,
     *   LEPCON, LSTOL, LEPIN, LPRT, LNZMX,
     *   LALPHA, LBETA, LGAMMA, LPSI, LSIGMA, LDELTA
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4, 
     *   LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12,
     *   LEPCON=13, LSTOL=14, LEPIN=15, LPRT=16, LNZMX=17,
     *   LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51)
C
C    Pointers into RCONT (a subsection of RWORK at the beginning of the buffer)
C
      INTEGER RCONTLEN
      PARAMETER (RCONTLEN=10)
      INTEGER LTSTART, LINITSTEP, LCJ, LTOUTB, LTNEXT, LTB, LTBRET
      PARAMETER (LTSTART=1, LINITSTEP=2, LCJ=3, LTOUTB=4, LTNEXT=5, 
     *   LTB=6, LTBRET=7)
C 
C    Pointers into ICONT (a subsection of IWORK between the IWORKF and IWORKB parts)
C
      INTEGER ICONTLEN
      INTEGER LLIWORKB, LLISE, LLSE, LonePass, LNCHK, LNCHKSTP,
     *        LLYF, LLYPF, LLPHIS, LLRWKS, LLTCHK, LLF_P,
     *        LLADY, LLADYP, LLY, LLT, LLYP, LLRWEND, LLRWF,
     *        LTDIR, LLPHI, LLIWB, LLRWB, LLADI, LLTMP, LIFWRITE,
     *        LINFOB, LINFO, LCHKPOS, LICHK, LLRCONT
      PARAMETER (ICONTLEN=90)
      PARAMETER (LLIWORKB=1,LLISE=2,LLSE=3,LonePass=4,LNCHK=5,
     *        LNCHKSTP=6,
     *        LLYF=7, LLYPF=8, LLPHIS=9, LLRWKS=10, LLTCHK=11, 
     *        LLF_P=12,
     *        LLADY=13, LLADYP=14, LLY=15, LLT=16, LLYP=17, 
     *        LLRWEND=18, LLRWF=19, LTDIR=20, LLPHI=21,
     *        LLIWB=22, LLRWB=23, LLADI=24, LLTMP=25, LIFWRITE=26,
     *        LCHKPOS=27,LICHK=28,LLRCONT=29,LINFOB=31,LINFO=61)
     
      LOGICAL ldebug
C      data    ldebug/.true./
      data    ldebug/.false./

      IF (INFOI(1) .EQ. 0) THEN
C
C-----------------------------------------------------------------------
C     This block is executed for the initial call only.
C     It contains checking of inputs and initializations.
C-----------------------------------------------------------------------

c       Initialization for the adifor3.0 tape
        if (infobi(8) .eq. 1) call ad3Initialize()

        TSTART = T
      ELSE
C
C-----------------------------------------------------------------------
C     This block is executed for continuation calls only.
C     It contains initializations from saved values.
C-----------------------------------------------------------------------
C
        LICONT  = IWORK(39)
        LRCONT  = IWORK(LICONT+LLRCONT)
        TSTART  = RWORK(LRCONT+LTSTART)
      END IF
      
      ires = 0
      if (nq .lt. 1) then
         print *, ' Improper value for NQ in DASPKADJOINT. NQ =', NQ
         IDID = -33
         RETURN
      END IF
      NSENPAR = infoi(22)
      NPARA   = infoi(19)
      if (nsenpar .gt. 0 .and. (infobi(3).gt.3)) then
         print *, ' when info(22)>0, infob(3) must be less than 4.'
         idid = -33
         return
      end if
c         
      do i = 1, 30
         info(i)  = infoi(i)
         infob(i) = 0
      end do
      info(3)  = 0 !ensure this is a legal value for DDASPK
      info(19) = 0
      istol  = infobi(1)
      itemp = 1
      if (istol .gt. 1 .or. istol .lt. 0) goto 999
      icopt  = infobi(2)
      itemp = 2
      if (icopt .gt. 3 .or. icopt .lt. 0) goto 999
      meva   = infobi(3)
      itemp = 3
      if (meva .gt. 7 .or. meva .lt. 0) goto 999     
      iquadv = infobi(4)
      itemp = 4
      if (iquadv .gt. 1 .or. iquadv .lt. 0) goto 999     
      index2 = infobi(5)
      itemp = 5
      if (index2 .gt. 1 .or. index2 .lt. 0) goto 999
      itemp = 6
      if (infobi(6) .gt. 1 .or. infobi(6) .lt. 0) goto 999
      itemp = 7
      if (infobi(7) .lt.0 .or. infobi(7) .gt. 2000) goto 999

      IF (INFOI(1) .EQ. 0) THEN
c
c...    Compute the consistent initial conditions for the forward 
c       integration at tstart
c    
      	if (info(11) .gt. 0) then
        	info(14) = 1           ! compute the initial conditions only 
         	call DDASPK(
     *        RES, NEQ, T, Y, YPRIME, TFINAL, INFO, RTOL,ATOL,
     *        IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL, 
     *        SENPAR, RES_ADY, k_res, T_RES, RES_ADP)
         	if (idid .lt. 0) return
         	if (ldebug) print *, 'Initialization done for the forward'
         	info(11) = 0           ! consistent initial conditions obtained
      	end if
c
c...    Compute the pointers for iworkb, and rworkb
        call daspkPointer(neq, info, iwork, rwork, liw, lrw, 
     *       lise, lse, idid)
        if (idid .lt. 0) return
      END IF
      
      if (iquadv .eq. 0 .and. nsenpar .gt. 0) then
         neqpns = neq + nsenpar
      else
         neqpns = neq
      end if
      if (neqad .ne. nq*neqpns) then
         print *, ' NEQAD is ', NEQAD, ', but should be ', NQ*NEQPNS
         idid = -33
         return
      end if

c
c...  Set points for IWORK(*)
      if (info(5) .eq. 2 .and. info(12) .eq. 0) then
c
c...  Adifor options, the 3*neq is the last part of the Iwork(*)
         liwks = iwork(lniw) - 3*neq + 1
      else
c
c...  finite difference options
         liwks = iwork(lniw) + 1
      end if
    
      IF (INFOI(1) .EQ. 0) THEN
        initstep = 0.d0
c
c...    save the length of the rwork for forward problem
        lgrwk = iwork(lnrw)
        infob(2)  = info(2)
        infob(4)  = 1             ! set up the tstop for each checkpoint
        infob(5)  = info(5)       ! Jacobian evaluation method
        infob(6)  = info(6)
        infob(9)  = info(9)       ! mxord is the same as for state variables
C
C...    always using the finite difference for the adjoint KMV
C
        infob(12) = info(12)
        if (infob(12) .eq. 1) infob(5) = 0
        infob(13) = info(13)
        infob(15) = info(15)
        infob(16) = info(16)      ! error control
        if (meva .eq.2 .or. meva.eq.6) then
           infob(16) = 1
        end if
        infob(18) = info(18)      ! exclude the quadrature variables from 
                                  !   the error test.
        infob(24) = -nq           ! backward integration
c
c...    adjoint system for index-2 DAEs
        if (index2 .eq. 1) infob(25) = 1
        if (iquadv .eq. 0) then
           infob(22) = nsenpar    ! sensitivity parameters in RES
           infob(28) = nsenpar + info(28)
        else
           infob(22) = 0
           infob(28) = info(28)
        end if
        infob(30) = meva          ! evaluation method for adjoint equations
c
c...    Compute the work space for the backward integration
        nbufloc = 0
        nbufh = 0                 ! buffer for checkpoint data
        nchk = 0                  ! number of check points
        if (meva .lt. 4) then
           call daspkPointer(neqad, infob, 
     *          iwork, rwork, liw, lrw, 
     *          lise, lse, idid)
           if (idid .lt. 0) return
           lgrwkb = iwork(lnrw)
           if (infobi(7) .gt. 0) then
c...  number of steps in forward integration is known
              if (meva .eq. 2) then
                 lgmx = max(lgrwkb, lgrwk) + nsenpar*nq + 4*neqad + 
     *                  1 + RCONTLEN
              else
                 lgmx = max(lgrwkb, lgrwk) + nsenpar*nq + 2*neqad + 
     *                  1 + RCONTLEN
              end if
              nbufloc = lrw - lgmx
           end if
           if (meva .eq. 2) then
              lgmx = lgrwkb + lgrwk + nsenpar*nq + 4*neqad + 
     *               1 + 2*neq + RCONTLEN
           else
              lgmx = lgrwkb + lgrwk + nsenpar*nq + 2*neqad + 
     *               1 + 2*neq + RCONTLEN
           end if
           nbufh = lrw - lgmx
        end if
c
      onePass = 0               ! indicator where only one forward integration
                                ! is needed (onePass = 1)
        if (meva .gt. 3) then
c...  for linear problem
           onePass = 1
           nchkstp = 0
           liworkb = 0
        else
c...  for nonlinear problem
           if (infobi(7).gt.0.and.nbufloc.gt.
     *         (infobi(7)+1)*(2*neq+1)) then
c...  if the number of steps of the forward step is known and the buffer is
c     larger enough to hold the data for all steps
              onePass = 1
              nchkstp = infobi(7)
              liworkb = 0
           else
c...  number of time steps between two checkpoints, 
c     checkpoints are set at every 12 steps
              nchkstp = 3*(nbufh-5-2*neq)/(13*neq + 39)
c     
c...    maximum number of checkpoints allowed in the dynamic memory
              nchk = nchkstp/3
              if (nchk .le. 0) then
                print *, ' RWORK is too short to hold the buffer.'
                idid = -33
                return
              end if
              licont  = liwks + 4*nchk
              liworkb = licont + icontlen
           end if
        end if
C
C...    set pointers in RWORK(*)
c...    set pointers for saving Y and YPRIME, adY and adYPRIME
c
c..   lyf, lypf: pointers used for backward integration 
        lyf = 1
        lypf = 1
        lphis = 1
        lrwks = 1
        ltchk = 1
c
c...  save the forward information at the last part of the rwork(*)
        if (onePass .eq. 0) then
           lyf     = lrw - neq
           lypf    = lyf - neq
           lf_p    = lypf - nsenpar*nq
        else
           lf_p    = lrw - nsenpar*nq
        end if
        if (meva .eq. 2 .or. meva .eq. 6) then
c...  nonlinear DAE and the F_y' is not constant
           lady    = lf_p - 2*neqad
           ladyp   = lady - 2*neqad
        else
           lady    = lf_p - neqad
           ladyp   = lady - neqad
        end if
      
c
c...  check point information, start from "ltchk" (onePass=0) or
C     "ladyp" (onePass=1)
        if (onePass .eq. 0) then
           lphis = ladyp - 7*nchk*neq
           lrwks = lphis - 35*nchk
           ltchk = lrwks - (nchk+1)
           lrwend = ltchk - 1
        else
           lrwend = ladyp - 1
        end if      
c
c...  "lrwend" is the last rwork(*) space that can be used by DASPK integration
c     either forward or backward

c
c...  lyp, ly, lt save the data for the intermediate forward integrations.
        lyp = lrwend
        ly = lyp
        lt = ly
        if (meva .lt. 4) then
c...  nonlinear DAE
           lyp = lrwend - (nchkstp + 1)*neq
           ly  = lyp -  (nchkstp + 1)*neq
           lt  = ly - (nchkstp + 1)
           lrwend = lt - 2
        end if
        
c...   save space for rcont (this puts it at the beginning of the buffer, I think)
        lrcont = lrwend - rcontlen
        lrwend = lrcont - 1
c
c...  lrworkb = beginning of the rwork(*) space for backward integration
        if (onePass .eq. 0) then
           lrworkb = lgrwk + 1
           lrwf    = lrworkb
        else
           lrworkb = 1
           lrwf    = lrwend
        end if
c
c...  initial stepsize (double)
        if (info(8) .eq. 1) initstep = rwork(3);
c         
        if (meva .lt. 4) then     ! forward information is needed
           if (onePass .eq. 1) then
c...  only one forward integration is required
              ltloc = lt
              lyloc = ly
              lyploc = lyp
c
c...  t = tstart information
              rwork(ltloc) = tstart
              do i = 0, neq-1
                 rwork(lyloc+i)  = y(i+1)
                 rwork(lyploc+i) = yprime(i+1)
              end do
              ltloc  = ltloc + 1
              lyloc  = lyloc + neq
              lyploc = lyploc + neq
           else            
c
c...  required intermediate forward integration
c
c..   t = tstart information
              ichk = 0            ! order number of the checkpoints
c     
c...    Save the current values of Y and YPRIME, tstart if necessary
              rwork(ltchk+ichk) = tstart
              do i = 0, 3
                 iwork(liwks+i) = 1
              end do
              do i = 0, 34
                 rwork(lrwks+i) = 0.0d0
              end do
              do i = 1, neq
                 rwork(lphis + i-1) = Y(i)
                 rwork(lphis + i-1 + neq) = YPRIME(i)
              end do
              ichk = ichk + 1
              ifwrite = 0         ! no outfile is required
              chkPos = 0
c     
c...    Pointers for the forward integration
           end if
c     
c...    Pointers for the forward integration
           lphi = ldelta + 3*neq
           IF (INFO(12) .NE. 0) lphi = lphi + neq
           IF (INFO(16) .NE. 0 .OR. (INFO(18).EQ.1.AND.INFO(28).GT.0)) 
     *        lphi = lphi + NEQ
        end if
c
c...    integration direction
        if (t .gt. TFINAL) then
           tdir = -1
        else 
           tdir = 1
        end if

        IF (onePass .eq. 1) THEN
C          ICONT begins after the standard IWORK, and then ISENWK.
C          ISENWK contains the IEOPT values for the sensitivity info and
C          is of length <= NEQPNS+NEQ+2.  
C          If there is a bug in the future, it may be caused by this computation.
C          This is designed to ensure each subarray of IWORK does not overlap with
C          any others, so if there is a problem, this is the first thing 
C          to check.     
           licont  = LISE+NEQPNS+NEQ+2
        END IF
        IWORK(39)              = LICONT
        IWORK(40)              = LRWORKB
        IWORK(LICONT+LLIWORKB) = LIWORKB
        IWORK(LICONT+LonePass) = onePass
        IWORK(LICONT+LLISE)    = LISE
        IWORK(LICONT+LLSE)     = LSE
        IWORK(LICONT+LNCHK)    = NCHK
        IWORK(LICONT+LNCHKSTP) = NCHKSTP
        IWORK(LICONT+LLYF)     = LYF
        IWORK(LICONT+LLYPF)    = LYPF
        IWORK(LICONT+LLPHIS)   = LPHIS
        IWORK(LICONT+LLRWKS)   = LRWKS
        IWORK(LICONT+LLTCHK)   = LTCHK
        IWORK(LICONT+LLF_P)    = LF_P
        IWORK(LICONT+LLADY)    = LADY
        IWORK(LICONT+LLADYP)   = LADYP
        IWORK(LICONT+LLY)      = LY
        IWORK(LICONT+LLT)      = LT
        IWORK(LICONT+LLYP)     = LYP
        IWORK(LICONT+LLRWF)    = LRWF
        IWORK(LICONT+LLRWEND)  = LRWEND
        IWORK(LICONT+LTDIR)    = TDIR
        IWORK(LICONT+LLPHI)    = LPHI
        IWORK(LICONT+LLRCONT)  = LRCONT
                
        RWORK(LRCONT+LINITSTEP) = INITSTEP
        RWORK(LRCONT+LTSTART)   = TSTART
                
        CALL FORWARD_INTEGRATE(
     *        RES, NEQ, T, Y, YPRIME, TFINAL, INFO, RTOL,ATOL,
     *        IDID, RWORK, LRWF, IWORK, LIW, RPAR, IPAR, JAC, PSOL, 
     *        SENPAR, RES_ADY, k_res, T_RES, RES_ADP,
     *        onePass, INFOBI, LTLOC, LYLOC, LYPLOC, ICHK, chkPos,
     *        LPHIS, LJCALC, LYF, LYPF, LRWKS, LPHI, TDIR,
     *        MEVA, NCHKSTP, NCHK, IFWRITE, LTCHK, LADYP, LIWKS,
     *        LALPHA, LNST, LHOLD, LKOLD)
     
        IWORK(LICONT+LIFWRITE) = IFWRITE
        
      ELSE
!       We are in a continuation call
        LICONT  = IWORK(39)
        LRWORKB = IWORK(40)
        LIWORKB = IWORK(LICONT+LLIWORKB)
        onePass = IWORK(LICONT+LonePass)
        LISE    = IWORK(LICONT+LLISE)
        LSE     = IWORK(LICONT+LLSE)
        NCHK    = IWORK(LICONT+LNCHK)
        NCHKSTP = IWORK(LICONT+LNCHKSTP)
        LYF     = IWORK(LICONT+LLYF)
        LYPF    = IWORK(LICONT+LLYPF)
        LPHIS   = IWORK(LICONT+LLPHIS)
        LRWKS   = IWORK(LICONT+LLRWKS)
        LTCHK   = IWORK(LICONT+LLTCHK)
        LF_P    = IWORK(LICONT+LLF_P)
        LADY    = IWORK(LICONT+LLADY)
        LADYP   = IWORK(LICONT+LLADYP)
        LY      = IWORK(LICONT+LLY)
        LT      = IWORK(LICONT+LLT)
        LYP     = IWORK(LICONT+LLYP)
        LRWF    = IWORK(LICONT+LLRWF)
        LRWEND  = IWORK(LICONT+LLRWEND)
        TDIR    = IWORK(LICONT+LTDIR)
        LPHI    = IWORK(LICONT+LLPHI)
        IFWRITE = IWORK(LICONT+LIFWRITE)
        
        LRCONT   = IWORK(LICONT+LLRCONT)
        INITSTEP = RWORK(LRCONT+LINITSTEP)
      END IF ! if this is the first time DASPKAdjoint is called

c------------------------------------------------------------------------------
c                 initialization for the backward integration
c------------------------------------------------------------------------------
c
c...  lsen is the pointer in QSEN for Q_P
      lsen = nq*(neq+1)+1
      tend = tfinal - tstart
      
      IF (INFOI(1) .EQ. 0) THEN
C       beginning of initialization for backward integration on first call
        do i = lsen, lsen + npara*nq - 1
           qsen(i) = 0.0d0
        end do
c
c...    Compute the Q and Q_Y, Q_P, stored in QSEN(nq, (neq+np+1)
        info24 = info(24)
        info(24) = -nq
        call DSENSDB(
     *       QRES, NEQ, Tfinal, Y, YPRIME, QSEN, INFO, 
     *       RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)
        info(24) = info24
c
c...    Initialize the INFOB for backward integration            
        if (istol .eq. 0) then
           if (info(2) .eq. 0) then
              rtolb(1) = 2.0d0*rtol(1)
              atolb(1) = 2.0d0*atol(1)
           else
              do i = 1, neq
                 rtolb(i) = 2.0d0*rtol(i)
                 atolb(i) = 2.0d0*atol(i)
              end do
           end if
        end if
c
c...  Note that the work space for the backward integration starts from 
c     liworkb+1 and lrworkb + 1
c
c...    Initialize the iworkb array, the same as iwork
        if (onePass .eq. 0) then
           do i = 1, licns-1
              iwork(liworkb+i) = iwork(i)
           end do
           if (info(6) .eq. 1) then
              iwork(liworkb+lml) = iwork(lmu)
              iwork(liworkb+lmu) = iwork(lml)
           end if
        else if (info(6) .eq. 1) then
c...  banded Jacobian information
           itemp = iwork(lml)
           iwork(lml) = iwork(lmu)
           iwork(lmu) = itemp
        end if      
c
c...    initialization options
        infob(11) = 0
        if (icopt .eq. 1) then
           infob(11) = 1
        else if (icopt .eq. 2) then
           infob(11) = 5          ! initialization separately
        end if
        if (infob(11) .gt. 0 .or. infob(16).eq.1) then
           do i = 1, neq
              iwork(liworkb+licns+i-1) = ieopt(i)
           end do
           do i = neq+1, neqpns
              iwork(liworkb+licns+i-1) = 1
           end do
        end if
        if (infob(11) .gt. 3) then
           do i = 1, neq
              if (ieopt(i+neq) .eq. -2) then
                 iwork(liworkb+licns+i-1+neqpns) = 1
              else
                 iwork(liworkb+licns+i-1+neqpns) = 0
              end if
           end do
           do i = neq+1, neqpns
              iwork(liworkb+licns+i-1+neqpns) = 0
           end do
        end if
        liwb  = liw - liworkb
        lrwb  = lrw - lrworkb
c
c...    Compute the pointers of ISENWK and SENWRK for initialization 
C       of the adjoint system
        call daspkPointer(neqad, infob, 
     *       iwork(liworkb+1), rwork(lrworkb+1), liwb, lrwb, 
     *       lise, lse, idid)
        if (idid .lt. 0) then
           print *, ' Error return from backward integration:'
           if (LDEBUG) THEN
             print *, ' IWORK STARTED FROM: ', LIWORKB+1, 
     *          ', RWORK STARTED FROM:', LRWORKB+1
           END IF
           return 
        end if
c
c...    adjoint system for index-2 DAEs
        if (index2 .eq. 1 .or. icopt .eq. 1 .or. icopt .eq. 2) then
c
c...  check if ieopt(*) is input correctly
           do i = 1, 2*neq
              if(ieopt(i).eq.0 .or. ieopt(i).lt.-2 .or. 
     *                                    ieopt(i).gt.3)then
                 MSG='DASPKADJOINT -- ILLEGAL VALUE FOR IEOPT(*) AT I1.'
              CALL XERRWD(MSG,60,4,0,1,I,0,0,0.0D0,0.0D0)
                 idid = -33
                 return
              END IF
           end do
c
c...   forward ID(*)
           do i = 1, neq
              iwork(liworkb+lise+i+1) = ieopt(i+neq) 
           end do
c
c...   Backward ID(*)
           do i = 1, neq
              iwork(liworkb+lise+neqpns+1+i) = ieopt(i)
           end do
           do i = neq+1, neqpns
              iwork(liworkb+lise+neqpns+1+i) = 1
           end do 
        end if        
        if (icopt .eq. 1 .or. icopt .eq. 2) then
c
c...  require G_x information during initialization
           lg_x = iwork(liworkb+lnrw) - neqad + 1
           ladi = lg_x - 6*neqpns
           if (icopt .eq. 2 .or. meva .eq. 2) ladi = ladi - neqpns
        else
           ladi = 0
        end if
        ltmp = lrworkb - neqpns - 1
c     
c...    transform the backward integration into forward
        tb = 0.0d0
c
c...    save TFINAL at 19th element of rworkb.
        rwork(lrworkb+19) = tfinal
        if (icopt .eq. 3) then    ! user-input
C          ADINIT returns values in RWORK(LADY) and RWORK(LADYP)
           call adinit(tfinal, neqad, nq,Y,YPRIME,rwork(ladY),
     *          rwork(ladYP),QSEN,RPAR,IPAR,senpar)
        else if (icopt .eq. 0) then 
c...    default initialization: ADY = Q_Y, ADYP = 0.0 
c...    for ODE y' = f(t, y)
           do iq = 1, nq
              do j = 1, neq
                 rwork(lady+(iq-1)*neqpns + j - 1) = QSEN(j*nq + iq)
              end do
c
c...    parameter part 
              do j = neq+1, neqpns
                 rwork(lady + (iq-1)*neqpns + j-1) = 0.0d0
              end do
           end do
           do i = 1, neqad
              rwork(ladyp+i-1) = 0.0d0
           end do
        else if (icopt .eq. 1) then ! index-1 or index-0 case
           infob(14) = 1
c
c...    TSTOP
           rwork(lrworkb+1) = tend
c
c...  forward information
           do i = 1, neq
              rwork(lrworkb + ladi + i-1) = y(i)
              rwork(lrworkb + ladi + i-1 + neqpns) = yprime(i)
           end do
c
C...  First initialization with adRes + g_x
           infob(29) = 1          ! needed g_x
c     
c...  Initialization for backward integration.
           do i = 1, neqad
              rwork(ladY+i-1)  = 0.0d0
              rwork(ladyp+i-1) = 0.0d0
           end do
           do iq = 1, nq
              do j = 1, neq
                 rwork(ladyp+(iq-1)*neqpns + j-1) = QSEN(j*nq + iq)
c
c...    Passing g_x through rworkb(lg_x)
                 rwork(lrworkb+lg_x + (iq-1)*neq + j-1) = 
     *                                                QSEN(j*nq + iq)
c                 if (ieopt(j) .lt. 0)  
c     *                rwork(lady+(iq-1)*neqpns + j-1) = QSEN(j*nq + iq)
              end do
           end do
           if (meva.eq.2 .or. meva.eq.6) then
              call DDASPK(
     *             G_RES_ADYP, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *             tend, INFOb, RTOLb,ATOLb,
     *             IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1),
     *             LIWB, RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *             G_RES_ADY,RES_ADP)
           else
              call DDASPK(
     *           adres, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           end if 
           if (idid .lt. 0) then
              print *, ' ERROR RETURN FROM BACKWARD INITIALIZATION:'
              IF (LDEBUG) THEN
                print *, ' IWORK STARTED FROM: ', LIWORKB+1, 
     *                 ', RWORK STARTED FROM:', LRWORKB+1
              END IF
              return
           end if
 80        continue
           if (ldebug) print *, ' initialization stage 1 pass! \n'
c
c...    If some sensitivity parameters are in the RES routine
           if (nsenpar .gt. 0) then
c
c...    Evaluate g_p - ady*(F_p) at t = tfinal
              do k = 0, nq-1
                 do i = 0, neq-1 
                    rwork(ldelta+i) = rwork(lady+i+k*neqpns)
                 end do
                 do i = 0, nsenpar-1
                    rwork(lf_p+i) = 0.0d0
                 end do
                 call res_adp(t, y, YPRIME, cj, rwork(ltmp), 
     *              rwork(ldelta), ires, rpar, ipar, senpar, 
     *              rwork(lf_p))
                 do i = 0, nsenpar-1
                    qsen(lsen + i*nq+k) = qsen(lsen+i*nq+k)-
     *                                     rwork(lf_p+i)
                 end do
              end do
           end if
C         
C...    second initialization with adRes
           infob(29) = 2          ! no g_x involved
c
c...    Work only when ady = 0 for the differential variables in the 
c       first pass, or when the number of algebraic variables is equal
c       to the number of algebraic equations.
           do i = 1, neqad
             rwork(lady+i-1) = rwork(ladyp+i-1)
              rwork(ladyp+i-1) = 0.0d0
           end do
c
           if (meva.eq.2 .or. meva.eq.6) then
              call DDASPK(
     *           G_RES_ADYP, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           else
              call DDASPK(
     *           adres, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           end if 
           if (idid .lt. 0) return
           if (ldebug) print *, ' initialization stage 2 pass!'
c
C...  Reset INFO(14) to zero
           infob(14) = 0
           infob(29) = 0
        else if (icopt .eq. 2) then ! index-2 case
c
c...  Hessenberg index-2 or mixed index-1 and index-2 case:
c-----------------------------------------------------------------------
c...  1. solving for the index-2 variables ady = (B^*C^*)^{-1}g_y
c
c               \dot\lam_1 = \lam_1 * A  + \mu_1 * C
c                      0   = \dot\lam_1 * B  - \lam_1 * \dot B - g_y 
c
c     for mu_1 = g_y(CB)^{-1}, 
c         \dot\lam_1 = g_y(CB)^{-1}C, 
c     with initial values
c         lam_1 = 0 and \dot\lam = 0.
c     The lam will be fixed and  \dot\lam_1 and \mu_1  will be computed.
c
c
           infob(14) = 1
c...  TSTOP
           rwork(lrworkb+1) = tend
c
c...  forward information
           do j = 1, neq
              rwork(lrworkb + ladi + j-1) = y(j)
              rwork(lrworkb + ladi + j-1 + neqpns) = YPRIME(j)
           end do
           do i = 1, nq
              do j = 1, neq
c
c...  Passing g_x through rworkb(lg_x)
                 rwork(lrworkb+lg_x + (i-1)*neq + j-1) = QSEN(j*nq + i)
              end do
           end do
           infob(29) = 3          ! g_y needed
           do i = 1, neqad
              rwork(lady+i-1) = 0.0d0
              rwork(ladyp+i-1) = 0.0d0
           end do
           call DDASPK(
     *        adres, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *        tend, INFOb, RTOLb,ATOLb,
     *        IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *        RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *        G_RES_ADY,RES_ADP)
           if (idid .lt. 0) return
           if (ldebug) print *, ' 1st pass for index-2 initialization:'
           cj = 1.0d0
c
c...    save mu_1 for later use
           lysave = ldelta + neqad
           do i = 1, neqad
              rwork(lysave + i - 1) = rwork(lady+i-1)
           end do
c
c       check if there is index-1 variable. If true
c...    solving the adjoint equation for G = \int g
c                    adxp = A^* adx + C^* ady + g_x
c                      0  = B^* adxp - \dot B^* adx
c       for adxp and ady with 
c              adx = -previous adxp
c       fixed during this stage.
c
           index1 = 0
           j = 0
           do while (index1 .eq. 0 .and. j .lt. neq)
              j = j + 1
              if (ieopt(j) .eq. -1) index1 = 1
           end do
           do i = 1, neqad
              rwork(lady+i-1) = - rwork(ladyp+i-1)
           end do
           if (index1 .eq. 1 .and. nsenpar .gt. 0) then
              infob(29) = 4       ! g_x needed
              do i = 1, neqad
                 rwork(ladyp+i-1) = 0.0d0
              end do            
              if (meva.eq.2 .or. meva.eq.6) then
                 call DDASPK(
     *              G_RES_ADYP, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *              tend, INFOb, RTOLb,ATOLb,
     *              IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1),LIWB, 
     *              RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *              G_RES_ADY,RES_ADP)
              else
               call DDASPK(
     *              adres, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *              tend, INFOb, RTOLb,ATOLb,
     *              IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1),LIWB, 
     *              RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *              G_RES_ADY,RES_ADP)
              end if 
              if (idid .lt. 0) return
              if (ldebug) print *, 
     *          ' 2nd pass for index-2 initialization:'
           end if
c
c...    If some sensitivity parameters are in the RES routine
           g_t = 1.0d0
           if (nsenpar .gt. 0) then
c
c...  Evaluate g_p - ady*(F_p) at t = tfinal, only for index upto 1
              do k = 0, nq-1
                 do i = 0, neq-1
                    if (ieopt(i+1) .gt. -2) then
                       rwork(ldelta+i) = rwork(lady+i+k*neqpns)
                    else
                       rwork(ldelta+i) = 0.0d0
                    end if
                 end do
                 do i = 0, nsenpar-1
                    rwork(lf_p+i) = 0.0d0
                 end do
                 call res_adp(t, y, YPRIME, cj, rwork(ltmp), 
     *              rwork(ldelta),ires, rpar, ipar, senpar, 
     *              rwork(lf_p))
                 do i = 0, nsenpar-1
                    qsen(lsen + i*nq+k) = qsen(lsen+i*nq+k)-
     *                   rwork(lf_p+i)
                 end do
c
c...  Added another term to the sensitivity: -\mu_1*\dot(h}_p
c   
c---------------------------------------------------------------  
                 i1 = lrworkb+1+ldelta
                 i2 = lysave + (k-1)*neqpns ! addelta               
                 do i = 1, neq
                    if (ieopt(i+1) .eq. -2) then
                       rwork(ldelta+i-1) = rwork(i2+i-1)
                    else
                       rwork(ldelta+i-1) = 0.0d0
                    end if
                    rwork(i1+i-1) = 0.0d0
                 end do
                 do i = 1, nsenpar
                    rwork(ladyp+i-1) = 0.0d0
                    rwork(lf_p+i-1) = 0.0d0
                 end do                 
                 call G_RES_ADP(tb, g_t, y, YPRIME, yprime, 
     *              cj, rwork(ltmp), rwork(ldelta), rwork(i1), ires, 
     *              rpar, ipar, senpar, rwork(ladyp), rwork(lf_p))
                 do i = 0, nsenpar-1
                    qsen(lsen + i*nq+k) = qsen(lsen+i*nq+k)+
     *                                     rwork(lf_p+i)
                 end do
c            
c---------------------------------------------------------------
              end do
           end if
c
c...  2. calculuate v_1 = g_x - lam * F_x - mu_1 * \dot C
c                   lam = - \dot\lam_1 = -g_y(CB)^{-1}C    
c
           do i = 1, nq
              i1 = lady   + (i-1)*neqpns ! ady
              i2 = lysave + (i-1)*neqpns ! addelta
              i3 = ldelta + (i-1)*neqpns ! g_ady
              i4 = ladyp  + (i-1)*neqpns ! adyp
              i5 = lrworkb+1+ldelta + (i-1)*neqpns ! g_addelta
c
              do j = 1, neq
                 rwork(i4+j-1) = 0.0d0               
                 if (ieopt(j) .lt. 0) then
c...  mu = 0, 
                    rwork(i1+j-1) = 0.0d0
                 end if
              end do
c.............v_2 = lam * F_x, stored in ladyp
              call res_ady(tb, y, rwork(i4), yprime, cj, rwork(ltmp), 
     *           rwork(i1), ires, rpar, ipar, senpar)
c
              do j = 1, neq
                 rwork(i1 + j - 1) = 0.0d0
c..............lambda = 0, 
                 if (ieopt(j) .gt. -2) then
                    rwork(i2+j-1) = 0.0d0 ! only mu_1 is retained
                 end if
                 rwork(i3 + j - 1) = 0.0d0
                 rwork(i5 + j - 1) = 0.0d0
              end do
c
c.............v_1 = mu_1 * \dot C, stored in g_ady = rwork(i3)
              call G_RES_ADY(tb, g_t, y, YPRIME, rwork(i1), rwork(i3), 
     *           yprime, cj, rwork(ltmp), rwork(i2), rwork(i5), 
     *           ires, rpar, ipar, senpar)
c     
              do j = 1,neq
c
c............ v = g_x - lam * F_x - mu_1*\dot C, stored in lady
c
c             since C = h_x = -F_x, v = g_x - lam * F_x + \mu_1*\dot F_x
c
                 if (ieopt(j) .gt. 0) then
                    rwork(i1+j-1) = QSEN(j*nq + i) - rwork(i4+j-1) + 
     *                   rwork(i3+j-1) 
                 end if
              end do
           end do
c
c...  3. calculate adx = (I-C*(B*C*)^{-1}B*) \mu
c     
c...   solving for adx = \mu M^{-1} P
c     adxp M = A^* adx + C^* ady + v
c         0  = B^* adxp - \dot B^* adx
c     for adxp and ady with adx = 0 fixed during this stage.
c     
           infob(29) = 4          ! g_x needed
           do i = 1, nq
              i1 = lady   + (i-1)*neqpns ! ady
              i4 = ladyp  + (i-1)*neqpns ! adyp
              do j = 1, neq
c     
c...  Passing v to g_x
                 if (ieopt(j) .gt. 0) then
                    rwork(lrworkb+lg_x + (i-1)*neq + j-1) = 
     *                   rwork(i1+j-1)
                 end if
                 rwork(i1+j-1) = 0.0d0
                 rwork(i4+j-1) = 0.0d0
              end do                  
           end do
           if (meva.eq.2 .or. meva.eq.6) then
              call DDASPK(
     *           G_RES_ADYP, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1),LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           else
              call DDASPK(
     *           adres, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1),LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           end if 
           if (idid .lt. 0) return
           if (ldebug) print *, ' 3rd pass for index-2 initialization:'
c
c...  If some sensitivity parameters are in the RES routine
           if (nsenpar .gt. 0) then
c
c...  Evaluate \mu_2 * h_p at t = tfinal, only for index-2 equation
              do k = 0, nq-1
                 do i = 0, neq-1
                    if (ieopt(i+1) .eq. -2) then
                       rwork(ldelta+i) = rwork(lady+i+k*neqpns)
                    else
                       rwork(ldelta+i) = 0.0d0
                    end if
                 end do
                 do i = 0, nsenpar-1
                    rwork(lf_p+i) = 0.0d0
                 end do
                 call res_adp(t, y, YPRIME, cj, rwork(ltmp),
     &              rwork(ldelta), 
     &              ires, rpar, ipar, senpar, rwork(lf_p))
                 do i = 0, nsenpar-1
                    qsen(lsen + i*nq+k) = qsen(lsen+i*nq+k)-
     *                      rwork(lf_p+i)
                 end do
              end do
           end if
c
c...  4. solving the adjoint equation for g(x,y,p) 
c                adxp = A^* adx  + C^* ady
c                   0 = B^* adxp - \dot B* adx
c     for adxp and ady, with adx fixed during the solving.
           do i = 1, neqad
              rwork(lady+i-1) = rwork(ladyp+i-1)
           end do
           infob(29) = 2           ! no g_x is needed
           if (meva.eq.2 .or. meva.eq.6) then
              call DDASPK(
     *           G_RES_ADYP, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           else
              call DDASPK(
     *           adres, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           tend, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *           RPAR, IPAR, adJAC, PSOL,SENPAR,RES_ADY,K_res,
     *           G_RES_ADY,RES_ADP)
           end if 
           if (idid .lt. 0) return
           if (ldebug) print *, ' 4th pass for index-2 initialization:'
 100       continue
c
C...  Reset INFO(14) to zero
           infob(14) = 0
           infob(29) = 0
        end if
        do j = neq+1, neqpns
c
c...  quadrature part
           do iq = 1, nq
              rwork(lady + (iq-1)*neqpns + j-1) = 0.0d0
           end do
        end do
c
C...    Reset INFO(11) to zero: consistent initial condition obtained
        infob(11) = 0
        lrwb = lrwend - lrworkb
c------------------------------------------------------------------------------
c                  prepare for the backward integration
c------------------------------------------------------------------------------
c...  Compute the pointers again for ISENWK and SENWRKB with new INFOB(11)
        call daspkPointer(neqad, infob, 
     *     iwork(liworkb+1), rwork(lrworkb+1), liwb, lrwb, 
     *     lise, lse, idid)
        if (idid .lt. 0) then
           print *, 'Returning from initialization for',
     *              ' backward integration'
           IF (LDEBUG) THEN
             print *, ' IWORK STARTED FROM: ', LIWORKB+1, 
     *          ', RWORK STARTED FROM:', LRWORKB+1
           END IF
           return
        end if
        if (meva .lt. 4) then                ! forward information is needed
c     
c...  set up the pointers for the internal forward integration data
          Ltpre=  lt - 1
c...  number of steps between two checkpoints , ISENWK(1) 
C..   Note, if you change this code, change the code calling
C..   dchi (to find return values of Y and YPRIME), because it
C..   is assuming isenwk(1) == nchkstp+1
          iwork(liworkb+lise) = nchkstp+1
c
c...  Difference between LSE and LNRW in backward integration, Stored in
C     ISENWK(2)
           iwork(liworkb+lise+1) = ltpre - lrworkb - lse
  	       if (iwork(liworkb+lise+1) .lt. 0) then
	          print *, ' length of rwork is not enough, short of ',
     *	            iwork(liworkb+lise+1)
	          idid = -33
	          return
	       end if
c...  ladi is a pointer in backward integration
            ladi = iwork(liworkb+lnrw) - 6*neqpns + 1
           if (meva .eq. 2) ladi = ladi - neqpns
c
c...  Store the interpolation time
           rwork(Ltpre) = tb   ! used in the adres evaluation
           if (meva .eq. 2) then
              do i = 1, neq
                 rwork(lrworkb + ladi + i - 1) = y(i)
                 rwork(lrworkb + ladi + neqpns + i - 1) = yprime(i)
              end do
           end if            
        end if
c
c...    adjoint system for index-2 DAEs
        if (index2 .eq. 1) then
           do i = 1, neq
              iwork(liworkb+lise+i+1) = ieopt(i+neq) 
           end do
           do i = 1, neq
              iwork(liworkb+lise+neqpns+1+i) = ieopt(i)
           end do
           do i = neq+1, neqpns
              iwork(liworkb+lise+neqpns+1+i) = 1
           end do
        end if
        
c...    Compute ady*F_p if necessary
        if (nsenpar .gt. 0 .and. info(28) .eq. infob(28)) then
c
c...    compute ady*F_p by adjoint code, saved in rwork(lf_p)
c       only for iquadv = 1
           do k = 0, nq-1
              do i = 0, neq-1
                 rwork(ldelta+i) = rwork(lady+i +k*neqpns)
              enddo
              do i = 0, nsenpar-1
                 rwork(lf_p+i+k*nsenpar) = 0.0d0
              end do
              call res_adp(t, y, YPRIME, cj, rwork(ltmp), 
     &           rwork(ldelta), 
     &           ires, rpar, ipar, senpar, rwork(lf_p+k*nsenpar))
           end do
        end if

        if (onePass .eq. 0) then
           chkPos = chkPos - 1
           rwork(ltchk+ichk) = TFINAL
        end if
        TOUTB = TB
        tnext  = tstart

        IWORK(LICONT+LLIWB) = LIWB
        IWORK(LICONT+LLRWB) = LRWB
        IWORK(LICONT+LLADI) = LADI
        IWORK(LICONT+LLTMP) = LTMP
        IWORK(LICONT+LCHKPOS) = chkPos
        IWORK(LICONT+LICHK)   = ICHK
        DO I = 1, 30
          IWORK(LICONT+LINFOB+I-1) = INFOB(I)
          IWORK(LICONT+LINFO+I-1)  = INFO(I)
        END DO
        
        RWORK(LRCONT+LCJ) = CJ
        RWORK(LRCONT+LTOUTB) = TOUTB
        RWORK(LRCONT+LTNEXT) = TNEXT
        RWORK(LRCONT+LTB)    = TB
        TBRET                = TB
        RWORK(LRCONT+LTBRET) = TBRET
C       end of initialization for backward integration on first call
      ELSE
C       beginning of initialization for backward integration on continuation call
        LIWB   = IWORK(LICONT+LLIWB)
        LRWB   = IWORK(LICONT+LLRWB)
        LADI   = IWORK(LICONT+LLADI)
        LTMP   = IWORK(LICONT+LLTMP)
        chkPOS = IWORK(LICONT+LCHKPOS)
        ICHK   = IWORK(LICONT+LICHK)
        DO I = 1, 30
          INFOB(I) = IWORK(LICONT+LINFOB+I-1)
          INFO(I)  = IWORK(LICONT+LINFO+I-1)
        END DO
        
        CJ    = RWORK(LRCONT+LCJ)
        TOUTB = RWORK(LRCONT+LTOUTB)
        TB    = RWORK(LRCONT+LTB)
        TBRET = RWORK(LRCONT+LTBRET)
        TNEXT = RWORK(LRCONT+LTNEXT)
C       beginning of initialization for backward integration on continuation call
      END IF ! end initialization stuff 

C     We will need to know how large each element in the checkpoint file is.
C     We query C for this info.      
      CALL GETSIZEOF(1,isized)
      CALL GETSIZEOF(0,isizei)

C     so far, so good

c------------------------------------------------------------------------------
c                       begin backward integration
c------------------------------------------------------------------------------
      if ((tfinal-tout)*tdir .lt. 0) then
        print *, 'tout = ', tout, ' occurs after final time ', tfinal
        idid = -33
        return
      end if
      tgoalb  = tfinal - tout
      if ((tgoalb-tend)*tdir .gt. 0) then
        tgoalb = tend
      end if

      IF (INFOI(3) .GT. 0) THEN
        IF (INFOI(1) .EQ. 0) THEN
C         This means that we are still in the first call, and that
C         we need to return the "initial conditions" at time tfinal.
          IF (INFOI(3) .EQ. 2 .AND. (TOUT-TFINAL)*tdir .ge. 0) THEN 
            DO J = 1, NEQAD
               ADY(J)  = RWORK(lADY+J-1)
               ADYP(J) = RWORK(lADYP+J-1)
            END DO
            T = TFINAL
            INFOI(1) = 1
            RETURN
          END IF
        ELSE ! infoi(3) > 0 and infoi(1) ~= 0
          IF (INFOI(3) .eq. 2 .and. (tgoalb-tbret)*tdir .le. 0) THEN
            print *, 'TOUT out of range: previous return value was',
     *               TFINAL-TBRET, ', TOUT is ', TFINAL-TOUT
            idid = -33
            return
          end if
          IF (INFOI(3) .eq. 2 .and. (tgoalb-tb)*tdir .le. 0) THEN
C           This could happen if DASPK is taking large steps, but the
C           user wants small ones. We are simplying going to interpolate 
C           the result.
C           Compute the return values for ADY and ADYP
            TAB = TBRET
            do j = 1, NEQad
                ADYA(j) = ADY(j)
                ADYPA(j) = ADYP(j)
            end do
            call dchi2(tgoalb, NEQAD, tab, tb, ADYA, 
     *             RWORK(LADY), ADYPA, RWORK(LADYP), ADY, ADYP)
C           Compute the return values for Y and YPRIME
C           This code is similar to that in ddresad, which 
C           DDASPK uses to compute the values of Y and YPRIME
C           for the adjoint residual routine.
            call dchi(TFINAL-TGOALB, NEQ, NCHKSTP+1, RWORK(LT), 
     *                RWORK(LY), RWORK(LYP), Y, YPRIME)
C           Do some bookkeeping
            TBRET = TGOALB
            RWORK(LRCONT+LTBRET) = TBRET
            T = TFINAL - TBRET
            RETURN
          END IF
          IF (INFOI(3) .eq. 1 .and. (TBRET-TB)*tdir .lt. 0) THEN
C           This could happen if the user switched from infoi(3) = 2 to 1
C           (i.e. the user switched from requesting output at set times
C            to requesting the next DASPK-chosen output time).
             DO J = 1, NEQAD
                ADY(J)  = RWORK(lADY+J-1)
                ADYP(J) = RWORK(lADYP+J-1)
             END DO
             call dchi(TFINAL-TB, NEQ, NCHKSTP+1, RWORK(LT), 
     *                 RWORK(LY), RWORK(LYP), Y, YPRIME)
C            Do some bookkeeping
             TBRET = TB
             T = TFINAL - TBRET
             RWORK(LRCONT+LTBRET) = TBRET
            RETURN
          END IF
        END IF
      END IF

!       print *, '\ncalling with t', t
      
!       print *, 'nchk is ', nchk, ' chkPos is ', chkPos

      do while (tb .lt. tend)
         if (onePass.eq.0 .and. ((tb-toutb)*tdir .ge. 0)) then  
                                ! forward information is required
c-----------------------------------------------------------------------------
c                  internal forward integration
c-----------------------------------------------------------------------------
            if (ifwrite.eq.1 .and. ichk.eq.0) then
c     
c...  Restore the information from the file
               IF (chkPos .LT. 0) THEN
                    IDID = -67
                    WRITE(*,*) 'ERROR - returning bc IFPOS went neg'
                    RETURN
               END IF
               irec1 = chkPos*((ladyp-ltchk)*isized+4*nchk*isizei)
               irec2 = irec1 + (ladyp-ltchk)*isized
               call adread(rwork(ltchk), 1, ladyp-ltchk, irec1)
               call adread(iwork(liwks), 0, 4*nchk, irec2)
               chkPos = chkPos - 1
               ichk = NCHK
            end if     
            ichk = ichk - 1
            if (ichk.gt.0 .or. (chkPos.ge.0 .and. ifwrite.eq.1))then
c     
c...  restore the information of PHI(*), Order, Time step at tf_beg.
c     do the forward integration by DASPK
c
c...  forward time
               tf = rwork(ltchk+ichk)
c
c...  integer information
               do i = 1, 4
                  iwork(5+i) = iwork(liwks+i-1+ichk*4)
               enddo
c
c...  double precision information
               do i = 3, 7
                  rwork(i) = rwork(lrwks + i-3 + ichk*35)
               end do
               do i = 1, 30
                  rwork(lalpha+i-1) = rwork(lrwks+i+4+ichk*35)
               end do
               nphi = iwork(7)
               if (nphi .lt. 5) nphi = nphi + 1
               do j = 0, neq-1
                  do i = 0, nphi
                     rwork(lphi + j + i*neq) =
     *                    rwork(lphis + j + i*neq + ichk*7*neq)
                  end do
                  yprime(j+1)=rwork(lphis+j+(nphi+1)*neq+7*ichk*neq)
               end do
               iwork(ljcalc) = -1 ! re-evaluate the Jacobian
               info(1) = 1      ! continual integration
            else                ! new integration
               info(1) = 0
               tf = rwork(ltchk)
               do j = 1, neq
                  y(j)  = rwork(lphis + j - 1)
                  yprime(j) = rwork(lphis + j - 1 + neq)
               end do
               if (info(8) .eq. 1) rwork(3) = initstep
            end if
c     
c...  loop of forward DASPK time step
            toutf = tfinal - tb
            if (info(1) .eq. 0) toutf = tfinal
            IWORK(LNST)=0
            IWORK(LNRE)=0
            IWORK(LNJE)=0
            IWORK(LETF)=0
            IWORK(LNCFN)=0
            ltloc = LT
            lyloc = ly
            lyploc = lyp
            if (ldebug)print *, 
     *          '  internal forward integration at ', tf
            do while (iwork(lnst) .lt. nchkstp .and. (tf-toutf)*tdir<0)
c     
c...  save the Y, YPRIME, T information in the buffer (rworkb(*))     
               rwork(LTloc) = tf
               if (info(1) .eq. 1) then
                  do i = 0, neq-1
                     rwork(lyloc+i)  = rwork(lphi + i)
                  end do
               else
!                 This will be called only once, when the last checkpoint is
!                 read in. (so we are starting from the t0 in the forward problem)
                  do i = 0, neq-1
                     rwork(lyloc+i)  = y(i+1)
                  end do
               end if
               do i = 0, neq-1
                  rwork(lyploc+i) = yprime(i+1)
               end do
               ltloc  = ltloc + 1
               lyloc  = lyloc + neq
               lyploc = lyploc + neq
               call DDASPK(
     *              RES, NEQ, TF, Y, YPRIME, TOUTF, INFO, RTOL,ATOL,
     *              IDID, RWORK, LRWF, IWORK, LIW, RPAR,IPAR,JAC,PSOL, 
     *              SENPAR, RES_ADY, k_res,T_RES,RES_ADP)
               if (idid .lt. 0) return
               if (ldebug) write(*,19)iwork(lnst),tf,rwork(lhold),
     $              iwork(lkold),'f'
 19         format(1x,I5,1x,e13.7,1x,e12.6,1x,i2,A1)
            end do
            rwork(LTloc) = tf
            do i = 1, neq
               rwork(lyloc+i-1)  = y(i)
               rwork(lyploc+i-1) = yprime(i)
            end do
            tnext = rwork(ltchk+ichk)
            toutb = tfinal - tnext
            rwork(lrworkb+1) = toutb
         end if                 ! (onePass .eq. 0) and (tb .ge. toutb)
         if (onePass .ne. 0 .and. (tb .ge. toutb)) then
c           We need toutb set only if we are re-entering daspkAdjoint
c           in the (middle) of the while (tb < toutb) loop
            toutb = tfinal - tnext
            rwork(lrworkb+1) = toutb
         end if
c     
c...  loop of backward integration
         if (ldebug) print *,'backward integration:', 
     *        tfinal-tb, tfinal - toutb, chkPos, tnext
!          print *,'backward integration:', 
!      *        tfinal-tb, tfinal - toutb, chkPos, tnext
         if (nsenpar .gt. 0 .or. ldebug .or. 
     *                      infoi(3) .eq. 1 ) then
            infob(3) = 1
         else
            infob(3) = 0
         end if
         if (infoi(3) .eq. 2 .and. ((tgoalb-toutb)*tdir < 0)) then
C          If we are trying to get output at specific points and we
C          will get to the output point before we get to toutb,
C          then we must shorten our time of integration.
           infob(3) = 1
         end if
         if (ldebug) print *,'b integ:', 
     *        tb, tgoalb, toutb, tnext, rwork(lT)
         do while ((tb-toutb)*tdir < 0)
            if (infoi(3) .eq. 2 .and. ((tgoalb-toutb)*tdir < 0)) then
C             If we are trying to get the output at specific points
C             and there is a chance that the next step could bring
C             us past toutb, then we must save the info at this
C             step so we can interpolate the answer.
              TAB = TB
              do j = 1, NEQad
                ADYA(j) = RWORK(ladY+j-1)
                ADYPA(j) = RWORK(ladYP+j-1)
              end do
            end if
            call DDASPK(
     *           adRES, NEQad, Tb, rwork(ladY), rwork(ladYP),
     *           toutb, INFOb, RTOLb,ATOLb,
     *           IDID, rwork(LRWORKB+1), LRWB, iwork(LIWORKB+1), LIWB, 
     *           RPAR, IPAR, ADJAC, PSOL, SENPAR, RES_ADY, K_res, 
     *           G_RES_ADY, RES_ADP)
            if (idid .lt. 0) return
            if (ldebug) write(*,19) iwork(liworkb+lnst),tfinal-tb,
     *                  rwork(lrworkb+lhold), iwork(liworkb+lkold),'b'
            if (nsenpar .gt. 0 .and. info(28) .eq. infob(28)) then
c
c...  do the integral of ady*F_p: (a+b)*\delta t/2, for iquadv=1
               lf_p2 = ldelta + neq 
               do k = 0, nq-1
                  do i = 0, neq-1
                     rwork(ldelta+i) = rwork(lady+i+k*neqpns)
                  enddo
                  do i = 0, nsenpar-1
                     rwork(lf_p2 + i) = 0.0d0
                  end do
c
c...  compute ady*F_p by adjoint code, saved in rwork(lf_p2)
                  call res_adp(tfinal-tb, rwork(lrworkb+ladi),
     *                 rwork(lrworkb+ladi+neqpns), 
     *                 cj, rwork(ltmp), rwork(ldelta), 
     *                 ires, rpar, ipar, senpar, 
     *                 rwork(lf_p2))
c
c...  Compute the integral \int ady*F_p
                  do i = 0, nsenpar-1
                     ii = i + k*nsenpar
                     qsen(lsen + i*nq + k) = qsen(lsen + i*nq + k) -
     *                    (rwork(lf_p2+i) + rwork(lf_p+ii))/2.0d0
     *                    *rwork(lrworkb+lhold)
c
c...  replace the old one by the current one
                     rwork(lf_p+ii) = rwork(lf_p2+i)
                  end do
               end do
            end if
C           We may need to return now
            if ((INFOI(3) .eq. 1) .or.
     *          (INFOI(3) .eq. 2 .and. (tb-tgoalb)*tdir .ge. 0)) then
C              Return Y, YPRIME, ADY, and ADYP at time TB.
               IF (INFOI(3) .eq. 1) then
!                Find the value of Y at time TB
                 call dchi(TFINAL-TB, NEQ, NCHKSTP+1, RWORK(LT),
     *                     RWORK(LY), RWORK(LYP), Y, YPRIME)
C                Record the values of ADY and ADYP at time TB
                 DO J = 1, NEQAD
                    ADY(J)  = RWORK(lADY+J-1)
                    ADYP(J) = RWORK(lADYP+J-1)
                 END DO
               END IF
               TBRET = TB
               if (infoi(3) .eq. 2) THEN
                 if ((tb-tgoalb)*tdir .gt. 0) then
C                  We need to interpolate the value
                   call dchi2(tgoalb, NEQAD, tab, tb, ADYA, 
     *             RWORK(LADY), ADYPA, RWORK(LADYP), ADY, ADYP)
                   TBRET = TGOALB
                 else
                   DO J = 1, NEQAD
                      ADY(J)  = RWORK(lADY+J-1)
                      ADYP(J) = RWORK(lADYP+J-1)
                   END DO
                   TBRET = TGOALB
                 end if
                 call dchi(TFINAL-TBRET, NEQ, NCHKSTP+1, RWORK(LT),
     *                     RWORK(LY), RWORK(LYP), Y, YPRIME)
               end if
               RWORK(LRCONT+LTOUTB)  = TOUTB
               RWORK(LRCONT+LTNEXT)  = TNEXT
               RWORK(LRCONT+LTB)     = TB
               RWORK(LRCONT+LTBRET)  = TBRET
               IWORK(LICONT+LICHK)   = ICHK
               IWORK(LICONT+LCHKPOS) = chkPos
               DO I = 1, 30
                 IWORK(LICONT+LINFOB+I-1) = INFOB(I)
                 IWORK(LICONT+LINFO+I-1) = INFO(I)
               END DO
               T = TFINAL - TBRET
               INFOI(1) = 1
               RETURN
            end if
         end do                 ! do while(tb .lt. toutb)
c 100     tb = toutb
      end do                    ! do while(tb .lt. tend)
 997  continue
 
      if (infoi(3) .gt. 0) then
        print *, 'Output time ', tout, ' out of bounds'
        idid = -33
        return
      end if
C This is where the backward loop ends. 

c------------------------------------------------------------------------------
c                            restore Y, YPRIME and QSEN
c------------------------------------------------------------------------------
c
c...  Save the pointers of backward information
      IWORK(39) = licont
      IWORK(40) = lrworkb
      if (iquadv .eq. 0) then
c
c...  Save the integral \int ady*F_p
         do k = 0, nq-1
            do i = 0, nsenpar-1
               qsen(lsen + i*nq + k) = qsen(lsen + i*nq + k) -
     *              rwork(lady+neq+i+k*neqpns)
            end do
         end do
      end if
c     
c...  Copy the result to QSEN
      do iq = 1, nq
         do i = 1, neq
            QSEN(i*nq + iq) = rwork(ladY+i-1+(iq-1)*neqpns)
         end do
      end do
c
c...  restore the final state for Y and YPRIME
      if (onePass .eq. 0) then
         do i = 1, neq
            Y(i)  = rwork(lyf + i - 1)
            YPRIME(i) = rwork(lypf+ i - 1)
         end do
      end if
      
      CALL computeSensitivity(NQ, NEQ, LDELTA, MEVA, NEQPNS, 
     *                        IRES, TSTART, RWORK, RPAR, IPAR,
     *                        SENPAR, LADI, LRWORKB, CJ, LADY,
     *                        Y, QSEN, LSEN, INFOBI, NPARA, 
     *                        ADRES, LTMP)

      if (ifwrite.eq.1) call adclose()
      if (infobi(8).eq.1) call ad3finalize()
      return
 999  continue
c
c...  error messages
      MSG = 'DASPK--  ELEMENT (=I1) OF INFOBI VECTOR IS NOT VALID'
      CALL XERRWD(MSG,50,1,0,1,ITEMP,0,0,0.0D0,0.0D0)
      return
      end
C     END SUBROUTINE DASPKADJOINT
C
C -------------------------------------------------------------------------------- C
C                           SUBROUTINE FORWARD_INTEGRATE                           C
C -------------------------------------------------------------------------------- C
      SUBROUTINE FORWARD_INTEGRATE(
     *        RES, NEQ, T, Y, YPRIME, TFINAL, INFO, RTOL,ATOL,
     *        IDID, RWORK, LRWF, IWORK, LIW, RPAR, IPAR, JAC, PSOL, 
     *        SENPAR, RES_ADY, k_res, T_RES, RES_ADP,
     *        onePass, INFOBI, LTLOC, LYLOC, LYPLOC, ICHK, chkPos,
     *        LPHIS, LJCALC, LYF, LYPF, LRWKS, LPHI, TDIR,
     *        MEVA, NCHKSTP, NCHK, IFWRITE, LTCHK, LADYP, LIWKS,
     *        LALPHA, LNST, LHOLD, LKOLD)
      
      IMPLICIT NONE
      INTEGER INFO, NEQ, IDID, LRWF, IWORK, LIW, IPAR
      INTEGER onePass, INFOBI, LTLOC, LYLOC, LYPLOC, ICHK, chkPos
      INTEGER LPHIS, LJCALC, LYF, LYPF, LRWKS, LPHI
      INTEGER MEVA, NCHKSTP, NCHK, IFWRITE, LTCHK, LADYP, LIWKS
      INTEGER LALPHA, TDIR, LNST, LHOLD, LKOLD
      DOUBLE PRECISION T, TFINAL, Y, YPRIME, RTOL, ATOL, RWORK
      DOUBLE PRECISION RPAR, SENPAR
      DIMENSION INFO(*), Y(*), YPRIME(*), RWORK(*), IWORK(*)
      DIMENSION RPAR(*), IPAR(*), SENPAR(*), INFOBI(*)
      EXTERNAL RES, JAC, PSOL, RES_ADY, K_RES, T_RES, RES_ADP
      
C     Local variables
      INTEGER I, J, IREC1, IREC2, NPHI
      logical ldebug
      DATA LDEBUG /.false./

c------------------------------------------------------------------------------
c                      Begin forward integration
c------------------------------------------------------------------------------
      info(3) = 1               ! out at each time step
c
c...  Loop of the integration of DASPK
      if (ldebug) then
         write(*,*) ' begin forward integration : '
         write(*,*) '-steps-----time------timestep---order'
      end if
      do while ((t - TFINAL)*tdir < 0)
         call DDASPK(
     *        RES, NEQ, T, Y, YPRIME, TFINAL, INFO, RTOL,ATOL,
     *        IDID, RWORK, LRWF, IWORK, LIW, RPAR, IPAR, JAC, PSOL, 
     *        SENPAR, RES_ADY, k_res, T_RES, RES_ADP)
         if (idid .lt. 0) return
         if (ldebug) then
            write(*,19) iwork(lnst), t, rwork(lhold), iwork(lkold),'f'
 19         format(1x,I5,1x,e13.7,1x,e12.6,1x,i2,A1)
         end if
c
c...  The checkpoint is determined by the number of time steps
         if (onePass .eq. 1) then
c...  only one Pass of forward integration is needed
            if (iwork(lnst) .gt. infobi(7)) then
               idid = -33
               print *, ' The number of the time steps ',
     *              'input in infobi(7) is not correct.',
     *              IWORK(LNST), ' steps are required.'
               return
            end if
            rwork(LTloc) = t
            do i = 0, neq-1
               rwork(lyloc+i)  = y(i+1)
               rwork(lyploc+i) = yprime(i+1)
            end do
            ltloc  = ltloc + 1
            lyloc  = lyloc + neq
            lyploc = lyploc + neq
         else if (meva.lt.4 .and. mod(iwork(lnst), NCHKSTP).eq.0) then
c
c...  Save the data at the checkpoint, to the buffer first, 
c     store to the disk file if the buffer overflowed.
c
            if (ichk .eq. NCHK) then
c...  buffer overflow
               if (ifwrite.eq.0) then
                  call adopen() ! open temporary file for write         
                  ifwrite = 1   ! require to write to the file
               end if
c
c...  Save the data to the file by  C routine
c
c              adwrite(var, sizeof(var), dimension(var), position)
c
               call adwrite(rwork(ltchk), 1, ladyp-ltchk,irec1)
               call adwrite(iwork(liwks), 0, 4*nchk, irec2)
               chkPos = chkPos + 1
               ichk = 0
            end if
            if (ldebug) print *, ' --- save at t=', t, ' checkpoint', 
     &           ichk
            rwork(ltchk+ichk) = t
            do i = 1, 4
               iwork(liwks+i-1+ichk*4) = iwork(5+i)
            enddo
            nphi = iwork(7)
            if (nphi .lt. 5) nphi = nphi + 1
            do i = 3, 7
               rwork(lrwks + i-3 + ichk*35) = rwork(i)
            end do
            do i = 1, 30
               rwork(lrwks + i+4 + ichk*35) = rwork(lalpha+i-1)
            end do
            do j = 0, neq-1
               do i = 0, nphi
                  rwork(lphis + j + i*neq + ichk*7*neq) = 
     *                 rwork(lphi + j + i*neq)
               end do
               rwork(lphis+j+(nphi+1)*neq+7*ichk*neq) = YPRIME(j+1)
            end do
            
            iwork(ljcalc) = -1
            ichk = ichk + 1
         end if
      end do                    ! end of the forward integration

      if (onePass .eq. 0) then
c...  Save the final state for Y and YPRIME for back integration
         do i = 1, neq
            rwork(lyf +i-1) = Y(i)
            rwork(lypf+i-1) = YPRIME(i)
         end do
      end if
      end ! end subroutine FORWARD_INTEGRATE
C
C -------------------------------------------------------------------------------- C
C                           SUBROUTINE computeSensitivity                          C
C -------------------------------------------------------------------------------- C
      SUBROUTINE computeSensitivity(NQ, NEQ, LDELTA, MEVA, NEQPNS, 
     *                              IRES, TSTART, RWORK, RPAR, IPAR,
     *                              SENPAR, LADI, LRWORKB, CJ, LADY,
     *                              Y, QSEN, LSEN, INFOBI, NPARA, 
     *                              ADRES, LTMP)
c----------------------------------------------------------------------------
c              sensitivity calculation
c----------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER  NQ, NEQ, LDELTA, MEVA, NEQPNS, IRES
      INTEGER  IPAR, LADI, LRWORKB, LADY, LSEN, INFOBI
      INTEGER  NPARA, LTMP
      DOUBLE PRECISION TSTART, RWORK, RPAR, SENPAR, CJ
      DOUBLE PRECISION Y, QSEN
      DIMENSION RWORK(*), IPAR(*), RPAR(*), SENPAR(*), INFOBI(*)
      DIMENSION QSEN(*), Y(*), CJ(*)
      EXTERNAL  ADRES
      
c     Local vars
      INTEGER  K, LF_P2, LF_P3, LF_P4, J, I

c...  compute the sensitivity and store it in QSEN
      do k = 0, nq - 1
         lf_p2 = ldelta + neq
c     
c...  compute ADY*F_y' by adjoint code, saved in rwork(lf_p2)
         if (meva.eq.3.or.meva.eq.7) then
c     
c...  use input adres
            lf_p3 = ldelta + neq + neqpns
            lf_p4 = ldelta + neq + 2*neqpns
            ires = 0
            do j = 0, neq-1
               rwork(lf_p3+j) = rwork(lady+j+k*neqpns)
               rwork(lf_p4+j) = 0.0d0
            end do
            call adres(tstart,rwork(lf_p4),rwork(lf_p3),
     *           CJ, rwork(lf_p2),IRES,RPAR,IPAR,SENPAR, 
     *           rwork(lrworkb+ladi))            
         else if (meva.eq.1.or.meva.eq.5 .or.
     *           meva.eq.2.or.meva.eq.6) then
c     
c...  input RES_ADYP
            do j = 0, neq-1
               rwork(ldelta+j) = rwork(lady+j+k*neqpns)
               rwork(lf_p2+j) = 0.0d0
            enddo
            call adres(tstart, rwork(lrworkb+ladi),
     *           rwork(lrworkb+ladi+neqpns),rwork(lf_p2), 
     *           cj, rwork(ltmp), rwork(ldelta), 
     *           ires, rpar, ipar, senpar)
         else if (meva.eq.0.or.meva.eq.4) then
            do j = 0, neq-1
               rwork(lf_p2 + j) = rwork(lady+j+k*neqpns)
            end do              
         end if
         if (infobi(6) .eq. 1) then
            do i = 0, npara-1
c     
c...  calculate ADY*F_y'*(Y_0)_p
               do j = 1, neq
                  qsen(lsen+i*nq+k) = qsen(lsen+i*nq+k) + 
     *                 rwork(lf_p2 + j - 1)*y((i+1)*neq+j)
               end do
            end do
         else
            do i = 1, neq
               QSEN(i*nq + k+1) = rwork(lf_p2+i-1)
            end do
         end if
      end do
      end  ! end of subroutine computeSensitivity
C
C
C -------------------------------------------------------------------------------- C
C                           SUBROUTINE daspkPointer                                C
C -------------------------------------------------------------------------------- C
      subroutine daspkPointer(
     *     neq, info, iwork, rwork, liw, lrw, 
     *     lise, lse, idid)
      IMPLICIT NONE
      integer info(*), iwork(*), idid, lise, lse
      real*8  rwork(*)
      INTEGER NEQ, LIW, LRW
      INTEGER MYID, NUMPROCS, NP, MYNEQ, NPB, MXORD
      INTEGER ICNFLAG, NONNEG, LID, LSOFF, LENIC, LENID, ITMP1
      INTEGER LENPD, MBAND, MSAVE, LENRW, LENWP, LENIWP, LENIW
      INTEGER MAXL, ITMP2, ITMP3, ITMP, NY, NAVG, ICNFLG, LSAVR
      INTEGER LE, LWT, LVT, LPHI, LWM, LPHI_LEN, I, LTEMP
      
      INTEGER LML, LMU, LMTYPE, 
     *   LIWM, LMXORD, LJCALC, LPHASE, LK, LKOLD,
     *   LNS, LNSTL, LNST, LNRE, LNJE, LETF, LNCFN,
     *   LNCFL, LNIW, LNRW, LNNI, LNLI, LNPS,
     *   LNSE, LMITER, LMAXL, LKMP, LNRMAX, LLNWP,
     *   LLNIWP, LLOCWP, LLCIWP, LKPRIN,
     *   LMXNIT, LMXNJ, LMXNH, LLSOFF, LNPD, LNY, 
     *   LNLIS, LICNS
      INTEGER LTSTOP, LHMAX, LH, LTN,
     *   LHOLD, LS, LROUND, LEPLI, LSQRN, LRSQRN,
     *   LEPCON, LSTOL, LEPIN, LPRT, LNZMX,
     *   LALPHA, LBETA, LGAMMA, LPSI, LSIGMA, LDELTA

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
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4,
     *   LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12,
     *   LEPCON=13, LSTOL=14, LEPIN=15, LPRT=16, LNZMX=17,
     *   LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51)
      CHARACTER MSG*80

      MYID = INFO(26)
      IF (INFO(27) .EQ. 0) INFO(27) = 1 ! for the multiprocessing
      NUMPROCS = INFO(27)
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
            NP = NAVG + 1
         ELSE
            NP = NAVG
         END IF
         MYNEQ = IWORK(LNY)*(NP+1)
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
         NP = 0
      END IF
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
        LSOFF = IWORK(LLSOFF)
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
      IF (INFO(11) .GT. 3) LENID = 2*NY
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
               ITMP1 = NY*(3*NY)
            ENDIF
         ELSE
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
         IWORK(LLNIWP) = LENIWP
         LENIW = 40 + LENIC + LENID + LENIWP
      ELSE IF (INFO(12) .EQ. 1)  THEN
         IWORK(LMTYPE) = INFO(5)
         IF (INFO(5) .GT. 1) THEN
            IWORK(LMTYPE) = 2
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
      IF (NP .GT. 0 .AND. INFO(20) .EQ. 5) THEN
C      iwork space for the Sparse Compressed Row format 
         LENIW = LISE + NY + 1 + INFO(22) + RWORK(LNZMX)*(2*NY*NY)
      END IF
C       
      IF(INFO(16) .NE. 0 .OR. 
     *     (INFO(18).EQ.1.AND.INFO(28).GT.0)) LENRW = LENRW + MYNEQ
      ITMP2 = 0
      IF (INFO(11) .GT. 3) THEN
C     Work space for the index-2 initialization
         IF (INFO(20) .LT. 2) THEN
            ITMP2 = 2*NY
            IF (NP .GT. 0) ITMP2 = 2*MYNEQ + 4*NY + 2*INFO(24)
         ELSE
            ITMP2 = NY
            IF (NP .GT. 0) ITMP2 = 2*NY+INFO(22)
         END IF
      END IF
      ITMP3 = 0
      IF(NP.NE.0) THEN
C     Work space for the sensitivity evaluations
         IF (INFO(20) .LT. 2) THEN
            ITMP3 = 4*NY + 2*INFO(24)
         ELSE IF (INFO(20) .EQ. 3) THEN
            ITMP3 = NP*(2*NY+MAX(INFO(24),NY)+INFO(22))
         ELSE IF(INFO(20) .EQ. 4) THEN
            ITMP3 =  NP  
         ELSE IF(INFO(20) .EQ. 5) THEN
C     Space for    (DF/DP), (DF/DY,DF/DY')
            ITMP3 = NP*NY + RWORK(LNZMX)*(2*NY*NY)+NP 
         END IF
      END IF
      ITMP = MAX0(ITMP1, ITMP2, ITMP3)
      LENRW = LENRW + ITMP
C     Integer work space for seed matrix of ADIFOR with sparseLinC option.
C
      IF ((INFO(12).EQ.0.AND.INFO(5).EQ.2) .OR. ((NP.GT.0 .OR. 
     *     INFO(11).EQ.5) .AND. INFO(20).EQ.5) ) 
     *     LENIW = LENIW + 3*NY + INFO(22)
C     Backward adjoint method
      IF (INFO(24) .LT. 0) THEN
C     The 5*MYNEQ space is for the seeded matrix
C
         LENRW = LENRW + 6*NY
         IF (INFO(11) .GT. 0) THEN
C
C...  Require G_Y information in ADRES evaluation during initialization
            LENRW = LENRW + MYNEQ
         END IF
         IF (INFO(30) .EQ. 2 .OR. INFO(11) .EQ. 5) LENRW = LENRW + NY
         LENIW = LENIW + 2
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
      LISE = IWORK(LLCIWP) + IWORK(LLNIWP)
      IF (INFO(24).GE.0 .OR. (INFO(24).LT.0.AND.
     *     INFO(30).NE.2.AND.INFO(30).NE.6) ) THEN         
         LSAVR = LDELTA
         IF (INFO(12) .NE. 0) LSAVR = LDELTA + MYNEQ
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
         IF (INFO(12) .EQ. 0) THEN
            LSAVR = LDELTA
            LE  = LSAVR + 2*MYNEQ
         ELSE
            LSAVR = LDELTA + 2*MYNEQ
            LE  = LSAVR + MYNEQ
         END IF
         LWT = LE + 2*MYNEQ
         LVT = LWT
         IF (INFO(16) .NE. 0 .OR. 
     *        (INFO(18).EQ.1.AND.INFO(28).GT.0)) LVT = LWT + MYNEQ
         LPHI = LVT + 2*MYNEQ
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
      idid = 1
      return
c-------------------------------------------------------------
 703  MSG = 'DASPK--  MAXORD (=I1) NOT IN RANGE'
      CALL XERRWD(MSG,34,3,0,1,MXORD,0,0,0.0D0,0.0D0)
      GO TO 750
 704  MSG='DASPK--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)'
      CALL XERRWD(MSG,60,4,0,2,LENRW,LRW,0,0.0D0,0.0D0)
      GO TO 750
 705  MSG='DASPK--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)'
      CALL XERRWD(MSG,60,5,0,2,LENIW,LIW,0,0.0D0,0.0D0)
      GO TO 750
 717  MSG = 'DASPK--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ'
      CALL XERRWD(MSG,52,17,0,1,IWORK(LML),0,0,0.0D0,0.0D0)
      GO TO 750
 718  MSG = 'DASPK--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ'
      CALL XERRWD(MSG,52,18,0,1,IWORK(LMU),0,0,0.0D0,0.0D0)
      GO TO 750
 728  MSG = 'DASPK-- NEQ MUST EQUAL #OF STATE VAR.*(#OF SENS.VAR.+1)'
      CALL XERRWD(MSG,58,24,0,0,0,0,1,0.D0,0.0D0)
      GO TO 750
 750  continue
      idid = -33
      return
      end ! end subroutine daspkPointer
C
C
C -------------------------------------------------------------------------------- C
C                           SUBROUTINE DSENDDB                                     C
C -------------------------------------------------------------------------------- C
      SUBROUTINE DSENSDB(
     *     QRES, NEQ, T, Y, YPRIME, QSEN, INFO, 
     *     RWORK, IWORK, IRES, RPAR, IPAR, SENPAR)
C
C***BEGIN PROLOGUE DSENSD
C***DATE WRITTEN   980813   (YYMMDD)
C***REVISION DATE  981124   (Compatible with DASPK3.0)
C
C***PURPOSE  This subroutine is auxiliary to DASPK and can
C            be called directly by the user between calls to
C            DASPK.  Its purpose is to compute the sensitivity
C            of a derived quantity specified by the user
C            in the external subroutine QRES with respect
C            to the problem parameters, given the solution and
C            sensitivities output by DASPK at time T. 
C            This version is in double precision.
C            
C-------------------------------------------------------------------
C *Arguments:
C
C  QRES:EXT       This is the name of a subroutine which you
C                 provide to define the  derived quantity Q(t,y,y',p).
C                 The subroutine you provide must have the form
C
C                 SUBROUTINE QRES(T, Y, YPRIME, Q, IRES, RPAR, IPAR, SENPAR)
C
C                 to define the derived quantity
C                 Q = Q(T, Y, P)
C
C                 If INFO(20)=3 or 4 (ADIFOR is selected), QRES is
C                 the name of ADIFOR-generated routine. It must has form
C
C                 SUBROUTINE QRES(T,Y,G_Y,YPRIME,DELTA,G_DELTA,IRES,
C                *     RPAR,IPAR,SENPAR,G_SENPAR)                 
C
C
C  T:IN           This is the current value of the independent variable
C                 on output from DASPK.
C
C  Y(*):IN        This array contains the solution components and
C                 sensitivities at T, on output from DASPK.
C
C  YPRIME(*):IN   This array contains contains the solution and
C                 sensitivity derivatives at T, on output from DASPK.
C            
C  QSEN(*):OUT    This array of dimension NQ(Np+NY+1), where Np is the
C                 number of parameters (defined in INFO(22) in DASPK),
C                 will contain the values of Q and gradient Q_y, Q_p
C                 
C
C  INFO(*):IN     This is the vector of code options that you provided
C                 to DASPK.  Note that DSENSD will follow the options
C                 that you specified for DASPK.  In particular,
C                 if INFO(20) = 0, sensitivities of the derived quantities
C                 will be computed by second order central differences,
C                 if INFO(20) = 1, they will be computed by first order
C                 differences.  
C                 If INFO(20) = 2, you should not be using this routine
C                 and should compute the sensitivities
C                 of the derived quantities yourself, given T, Y, YPRIME
C                 from DASPK.  To use this routine successfully,
C                 you must set INFO(24)=NQ as input to DASPK.
C                 If INFO(20) > 3, sensitivities of the derived quantities
C                 will be computed by ADIFOR with seed matrix option.
C
C  RWORK(*):IN    This is the real work array on output from DASPK.
C                 Note that if you are using this routine (DSENSD)
C                 to compute sensitivities of derived quantities, 
C                 you will need to augment the length of the RWORK
C                 array which is specified in the documentation to
C                 DASPK by an additional 2*NQ locations if
C                 INFO(20) = 0 or 1.
C                                   
C  IWORK(*):IN    This is the integer work array on output from DASPK.
C
C  RPAR,IPAR:IN   These are real and integer parameter arrays which
C                 you can use for communication between your calling
C                 program and the RES and QRES
C
C  SENPAR(*):IN   This is the real work array which is used in DASPK to
C                 store problem parameters of the sensitivities.
C-------------------------------------------------------------------------
C
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       EXTERNAL QRES
       INTEGER NEQ, INFO, IWORK, IRES, IPAR
       DOUBLE PRECISION T, Y, YPRIME, QSEN, RWORK, RPAR, SENPAR
       DIMENSION Y(*),YPRIME(*),QSEN(*)
       DIMENSION INFO(*),RWORK(*),IWORK(*),RPAR(*),IPAR(*),SENPAR(*)
       INTEGER LE, LWT, LQ, NQ, NY, IP, IDROW, I, IDROWB, LSAVR
       DOUBLE PRECISION UROUND, SQUR, DEL, YSAVE, DELINV, PSAVE
C
C      Set pointers into IWORK as in DASPK.
C
      INTEGER LML, LMU, LMTYPE, 
     *   LIWM, LMXORD, LJCALC, LPHASE, LK, LKOLD,
     *   LNS, LNSTL, LNST, LNRE, LNJE, LETF, LNCFN,
     *   LNCFL, LNIW, LNRW, LNNI, LNLI, LNPS,
     *   LNSE, LMITER, LMAXL, LKMP, LNRMAX, LLNWP,
     *   LLNIWP, LLOCWP, LLCIWP, LKPRIN,
     *   LMXNIT, LMXNJ, LMXNH, LLSOFF, LNPD, 
     *   LNY, LICNS
      PARAMETER (LML=1, LMU=2, LMTYPE=4, 
     *   LIWM=1, LMXORD=3, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *   LNS=9, LNSTL=10, LNST=11, LNRE=12, LNJE=13, LETF=14, LNCFN=15,
     *   LNCFL=16, LNIW=17, LNRW=18, LNNI=19, LNLI=20, LNPS=21,
     *   LNSE=22, LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26, LLNWP=27,
     *   LLNIWP=28, LLOCWP=29, LLCIWP=30, LKPRIN=31,
     *   LMXNIT=32, LMXNJ=33, LMXNH=34, LLSOFF=35, LNPD = 36, 
     *   LNY=37, LICNS=41)
C
C     Set pointers into RWORK as in DASPK.
C
      INTEGER LTSTOP, LHMAX, LH, LTN, LCJ, LCJOLD,
     *   LHOLD, LS, LROUND, LEPLI, LSQRN, LRSQRN,
     *   LEPCON, LSTOL, LEPIN, LPRT,
     *   LALPHA, LBETA, LGAMMA, LPSI, LSIGMA, LDELTA
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4, LCJ=5, LCJOLD=6,
     *   LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12,
     *   LEPCON=13, LSTOL=14, LEPIN=15, LPRT=16,
     *   LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51)
C
C     Set up for multiprocessing
C
C.....Compute the number of equations for this processor
C
      LSAVR = LDELTA
      IF (INFO(12) .NE. 0) LSAVR = LDELTA + NEQ
      LE  = LSAVR + NEQ
      LWT = LE + NEQ
      LQ  = LWT + NEQ
C
C     Find the dimension of the derived quantities
      NQ = ABS(INFO(24))
      NY = IWORK(LNY)
C
C     sensitivity evaluated by ADIFOR
C
      IRES = 0
      IF (INFO(20).GT.2) THEN    ! seed matrix option
         CALL JbAdfSM(
     *        T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR,QRES,NY,INFO(19),
     *        RWORK(LQ))
         RETURN
      ELSE IF (INFO(20) .EQ. 2) THEN
         IRES = 1
         CALL QRES(T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR)
         RETURN
      ENDIF
C
C     forward finite difference scheme.
C
C     Initial call to user defined QRES routine.
      CALL QRES(T,Y,YPRIME,QSEN,IRES,RPAR,IPAR,SENPAR)
C
C     Iterate on the number of parameters, NP.
C
      UROUND = RWORK(LROUND)
      SQUR = SQRT(UROUND)*1000
C
C...  Compute Q_x
      DO IP = 1, NY
         IDROW = IP*NQ
C
C       Determine the perturbation.
C
         DEL=SQUR*MAX(ABS(Y(IP)), ABS(1.D0/RWORK(LWT+IP-1)))
         DEL=(Y(IP)+DEL)-Y(IP)
C
         YSAVE = Y(IP)
         Y(IP) = Y(IP) + DEL
C
C     Call RES with perturbed values.
C
         CALL QRES(T,Y,YPRIME,RWORK(LQ),IRES,RPAR,IPAR,SENPAR)
C
C     Approx. first order sensitivity, and restore Y, YP, parameter values.
C
         DELINV = 1.0D0/DEL
         DO I = 1,NQ
            QSEN(IDROW+I)=DELINV*(RWORK(LQ+I-1) - QSEN(I))
         ENDDO
C
C       Restore y, and parameter values.  Proceed to
C       next Ny sensitivity segment.
         Y(IP) = YSAVE
      ENDDO
C
C...  Compute Q_p
      IDROWB = NQ*(NY+1)
      DO IP = 1, INFO(22)
         IDROW = IDROWB + (IP-1)*NQ
	 IF (SENPAR(IP).EQ.0.D0) THEN
	    DEL = SQUR
	 ELSE
            DEL = SQUR*ABS(SENPAR(IP))
	 END IF
         DEL = SENPAR(IP) + DEL - SENPAR(IP)
         PSAVE = SENPAR(IP)
         SENPAR(IP) = SENPAR(IP) + DEL
         CALL QRES(T,Y,YPRIME,RWORK(LQ),IRES,RPAR,IPAR,SENPAR)
         DELINV = 1.0D0/DEL
         DO I = 1, NQ
            QSEN(IDROW+I) = DELINV*(RWORK(LQ+I-1) - QSEN(I))
         END DO
         SENPAR(IP) = PSAVE
      END DO
C         
      RETURN
      END
C-------------------------END OF DSEND ROUTINE----------------------
      SUBROUTINE JbAdfSM(
     *     T,Y,YP,DELTA,IRES,RPAR,IPAR,senpar,G_RES,NY,isenfo,wrk)
c=======================================================================
c
C  This version is for using the work array
C
c   How to use this routine:
C    1. Put all the routines related to 
C       SUBROUTINE qRES(T, Y, YPRIME, Q, IRES, RPAR, IPAR, SENPAR)
C       in a file called "qres.f"
C    2. Create file "qres.cmp" with one line:
C        qres.f
C    3. Create file "qres.adf":
C        AD_PROG = qres.cmp
C        AD_TOP = qres
C        AD_IVARS = y,senpar
C        AD_OVARS = Q
C        AD_PREFIX = q
C        AD_OUTPUT_DIR = jac
C
C    4. run Adifor generate q_qres(...) with
C       % Adifor AD_SCRIPT=qres.adf
C    
C    5. Call this routine with proper parameters
C
C**** Copyright (C) 1999, Shengtai Li
C         
C=======================================================================   
      IMPLICIT NONE
      real*8  T                 ! independent variable, time
      real*8  Y(*)              ! Solution
      real*8  YP(*)             ! yp(i) = dy(i)/dt
      real*8  cj                ! const coming from the DASPK
      integer ires              ! error indicator 
      real*8  rpar(*)           ! real parameter of the problems 
      real*8  senpar(*)        ! sensitivity parameter array
      integer ipar(*)           ! integer parameter of the problems
      integer isenfo(*)
      external g_res            ! residual generated by the Adifor
      integer ny                ! number of equations for state variables only
      real*8  delta(*)          ! work array for residuals, size = NY
      real*8  wrk(*)            ! work array, size = width*(3*NY)
c === local variable
      integer i, j, ig_y, ig_delta, ipos, length, itotal, np,nq,
     *        ig_senpar
c
c === define the width of the seed matrix
c
      np = NY + isenfo(4)
      nq = abs(isenfo(6))
c
c === set index for the work array
c
      length = nq*np
      ig_y = 0
      ig_senpar = ig_y + ny*np
      ig_delta = ig_senpar + isenfo(4)*np
      itotal = ig_delta + length
c
c === initialize the seed matrix
c
      do i = 1, itotal
         wrk(i) = 0.0d0
      end do
      DO I=1, ny
         ipos = I + (I-1)*np
         wrk(ig_y + ipos )=1.0D0
      end do
      DO I=1, isenfo(4)
         ipos = I + (I-1)*np + ny*np + ny
         wrk(ig_senpar + ipos )=1.0D0
      end do
c
c === call Adifor generated routine
c
      call g_res(np, t, y, wrk(ig_y+1), yp, 
     *     delta, wrk(ig_delta+1), ires,
     *     rpar, ipar, senpar, wrk(ig_senpar+1))
c
c...  transpose (ny, nq) --> (nq, ny)
      do i = 1, np
         do j=1,Nq
            delta(j + I*NQ)=wrk(ig_delta+ (j-1)*np + i)
         end do
      end do
c
      RETURN
      END ! end of subroutine JbAdfSM

C=======================================================================   
c...  Cubic Hermit Interpolation
C     This subroutine was copied from the one in DDASPK, and
C     then edited so it gets two 1-D arrays in place of 2D arrays.
C     It finds y(t) and yp(t) using the values of ya and ypa
C     It should be the case that tab < t < tae (or tab > t > tae).
C     And yab is y(tb), yae is y(te), ypab is yp(tb), yae is yp(te)
C=======================================================================   
      SUBROUTINE dchi2(t, NY, tab, tae, yab, yae, ypab, ypae, y, yp)
      implicit none
      integer NY, j
      double precision t, tab, tae, yab(NY), yae(NY)
      double precision ypab(NY), ypae(NY), y(NY), yp(NY)
      double precision dy, y1tl, y1tr, c1, c2, c3, c4, qdt, dif, dt
      
      dt = tae - tab
      qdt = 1.0d0/dt
      dif = t - tab
      do j = 1, NY
         dy = (yae(j) - yab(j))*qdt
         y1tl = ypab(j)
         y1tr = ypae(j)
         c1 = yab(j)
         c2 = y1tl
         c3 = (3.0d0*dy - y1tr - 2.0d0*y1tl)*qdt
         c4 = (y1tr + y1tl - 2.0d0 * dy) * qdt * qdt
         y(j) = c1 + dif * (c2 + dif*(c3 + dif*c4))
         yp(j)= c2 + dif * (2.0d0 * c3 + dif * 3.0d0 * c4)
      end do
      return
      end ! End of subroutine dchi2
