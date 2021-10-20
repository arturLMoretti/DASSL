#----------------------------------------------------------------------------
# This makefile compiles and loads the DDASPK example program 
# If necessary, change the constants COMP and FFLAGS below for the
# compiler to be used.
#
# (C) Shengtai Li
#----------------------------------------------------------------------------
#
#    TEST PROBLEM FOR THE DASPK3.0 AND DASPKADJOINT
#
#-----------------------------------------------------------------------------
.SUFFIXES: .F .f .o .c

# PATH to the fortran 77 compiler
COMP = g77
ADFC = adf77
# PATH to TAMC preprocessor command
TAMC = tamc
TAMCFLAGS = -module res -pure
# PATH to Adifor preprocessor command
ADF   = Adifor

# ADIFOR library for SparseLinC options
ADLIB    =  $(AD_LIB)/lib/libSparsLinC-$(AD_OS).a 
ADLIB3   =  $(AD_LIB)/lib/$(AD_OS)/libSparsLinC-$(AD_OS)86.a \
            $(AD_LIB)/lib/$(AD_OS)/libad.a
# Compiler flags
FFLAGS = -g -O2
CFLAGS = -g -O2

#
# change the following lines for other platform
#
ifeq ($(AD_OS),IRIX)
ADLIB  = $(AD_LIB)/lib/libSparsLinC-$(AD_OS)-o32.a 
FFLAGS = -O2
CFLAGS = -O2
endif

SOLVR  = ../../solver
# SOLVR  = $(HOME)/Matlab/solverc
ADSLVR = ../../ad_solver
PRECON = ../../preconds

#--------------------------------------------------------------------
# DASPK using ADIFOR packages (no MPI packages)
#--------------------------------------------------------------------
OBJS = $(SOLVR)/ddaspk.o $(SOLVR)/daux.o $(SOLVR)/dlinpk.o \
       $(SOLVR)/mpi_dummy.o
OBJSB = $(ADSLVR)/daspk_adjoint.o $(ADSLVR)/adwrite.o $(ADSLVR)/ad3setup_dummy.o
OBJSB3 = $(ADSLVR)/daspk_adjoint.o $(ADSLVR)/adwrite.o $(ADSLVR)/ad3setup.o \
         adroutine_template.o 

# derived informations for heat equation and ILU preconditioner
SENSD = $(SOLVR)/dsensd.o reshiludr.o

# Adifor-generated routines for different examples
#
AD_HEATD   = res4sen/h_resh.o res4sen/g_resh.o jacsp/i_resh.o \
             res4sen/h_reshdr.o
AD_HEATD2   = JM_Output/gh_resh.o JM_Output/g_resh.o \
              JM_Output/gh_accum.o JM_Output/g_accum.o\
              SJ_Output/sg_resh.o SJ_Output/sg_accum.o\
             JM_Output/gq_reshdr.o JM_Output/gq_accum.o
AD_HEATILU = res4sen/h_reshilu.o res4sen/g_reshilu.o jacsp/i_reshilu.o \
	     res4sen/h_reshiludr.o
AD_HOTILU = res4sen/h_reshot.o res4sen/g_reshot.o jacsp/i_reshot.o 
AD_WEBD = res4sen/h_resweb.o jac/k_resweb.o jacsp/i_resweb.o
AD_WEBILU = res4sen/h_reswebilu.o jacsp/i_reswebilu.o
AD_PEND = res4sen/h_respend.o jacsp/i_respend.o jac/j_respend.o \
          jac/g_respend.o res4sen/g_g_respend.o
AD_WHEEL = res4sen/h_wheel.o  jac/j_wheel.o jac/g_wheel.o res4sen/g_g_wheel.o
AD_CRANK = res4sen/h_crank.o  jac/j_crank.o jac/g_crank.o res4sen/g_g_crank.o \
	   res4sen/h_ddot.o  jac/j_ddot.o jac/g_ddot.o res4sen/g_g_ddot.o
AD_FEKETE = res4sen/h_fekete.o  jac/j_fekete.o jac/g_fekete.o \
	    res4sen/g_g_fekete.o res4sen/l_fekete.o
#
# TAMC generated routines
AD_PENDB  = jac/j_respend.o res4sen/g_g_respend.o \
            resad_pend.o respend_ady.o adj/g_respend_ady.o adj/g_resad_pend.o
AD_FEKETEB =  jac/j_fekete.o  fekete_ady.o  resad_fekete.o \
              adj/g_fekete_ady.o  adj/g_resad_fekete.o \
              res4sen/l_fekete.o res4sen/h_qres_fekete.o

HEAT = dheat.o $(AD_HEATD) $(OBJS) $(PRECON)/dbanpre.o $(SOLVR)/dsensd.o \
       reshdr.o 
HEAT2 = dheat2.o $(AD_HEATD2) $(OBJS) $(PRECON)/dbanpre.o $(SOLVR)/dsensd.o \
       reshdr.o 
HEATB = dheat_ad.o $(OBJS) $(OBJSB) jacsp/i_resh.o resh_ady.o resh_adp.o 
HEATILU = dheatilu.o $(AD_HEATILU) $(OBJS) $(SENSD) \
          $(PRECON)/dilupre_adf.o $(PRECON)/dsparsk.o 
HEATILUB = dheatilu_ad.o $(OBJS) jacsp/i_reshilu.o  \
          $(PRECON)/dilupre_adf.o $(PRECON)/ad_dilupre_adf.o \
          $(PRECON)/dsparsk.o $(OBJSB) reshilu_ady.o reshilu_adp.o

HEATB2 = dheat_ad2.o $(OBJS) $(OBJSB3) SJ_Output/sg_resh.o \
         SJ_Output/sg_accum.f \
         JtM_Output/ady_resh.o JtM_Output/adp_resh.o JtM_Output/adp_accum.o \
         JtM_Output/ady_accum.o
HEATILUB2 = dheatilu_ad2.o $(OBJS) SJ_Output/sg_reshilu.o  \
          SJ_Output/sg_accum.f \
          $(PRECON)/dilupre_adf.o $(PRECON)/ad_dilupre_adf.o \
          $(PRECON)/dsparsk.o $(OBJSB3) \
          JtM_Output/ady_reshilu.o JtM_Output/adp_reshilu.o \
          JtM_Output/adp_accum.o JtM_Output/ady_accum.o

HOTB = dhot_ad.o $(OBJS) jacsp/i_reshot.o  \
          $(PRECON)/dilupre_adf.o $(PRECON)/ad_dilupre_adf.o \
          $(PRECON)/dsparsk.o $(OBJSB) reshot_ady.o reshot_adp.o

BGILUB = burger.o $(OBJS) jacsp/i_res_burger.o  \
          $(PRECON)/dilupre_adf.o $(PRECON)/ad_dilupre_adf.o \
          $(PRECON)/dsparsk.o $(OBJSB) res_burger_ady.o res_burger_adp.o

WEB = dweb.o $(AD_WEBD) $(OBJS) $(PRECON)/drbdpre.o $(PRECON)/drbgpre.o 

WEBB = dweb_ad.o $(OBJS)  resweb_ady.o resweb_adp.o resweb_adyp.o \
       $(OBJSB) jacsp/i_resweb.o

WEBILU = dwebilu.o $(AD_WEBILU) $(OBJS) $(PRECON)/dilupre_adf.o \
         $(PRECON)/dsparsk.o 

PEND = pend.o $(AD_PEND) $(OBJS)
PENDB= pendb.o $(AD_PENDB) $(OBJS) $(OBJSB)
WHEEL = wheel_driver.o wheel.o $(AD_WHEEL)  $(OBJS)
CRANK = crank_driver.o crank.o $(AD_CRANK) $(OBJS)
FEKETE = fekete_driver.o fekete.o $(AD_FEKETE) $(OBJS) $(SOLVR)/dsensd.o
FEKETEB = feketeb_driver.o fekete.o $(AD_FEKETEB) $(OBJS) $(OBJSB) 
TESTAD = testad.o  $(OBJS) $(OBJSB)
TESTAD2 = testad2.o  $(OBJS) $(OBJSB)

# -----------------------------default --------------------------------------
#
default : 
	@echo "Usage: make OBJ "
	@echo " "
	@echo " where OBJ can be: "
	@echo "   solution and sensitivity by forward method: "
	@echo "         heat   --- Heat equation (index-0) by Direct method "
	@echo "         heatilu -- Heat equation by Krylov method with ILU  "
	@echo "                    preconditioner " 
	@echo "         web    --- food web (index-1) problem by Direct method "
	@echo "         webilu --- food web problem by Krylov method with ILU "
	@echo "                    preconditioner "
	@echo "         wheel  --- wheel set problem (mixed index-1 and index-2) "
	@echo "                    by Direct method "
	@echo "         crank ---- crank problem (mixed index-1 and index-2) by "
	@echo "                    Direct method "
	@echo "         fekete --- Fekete equation (index-2) by Direct method "
	@echo "         pend   --- single pendulum (index-2) by Direct method "
	@echo " "
	@echo "   sensitivity by adjoint method by TAMC: "
	@echo "         heatb  --- Heat equation by direct method "
	@echo "         heatilub - Heat equation by Krylov method with ILU  "
	@echo "                    preconditioner " 
	@echo "         webb   --- food web problem by direct method "
	@echo "         feketeb -- fekete equation by direct method "
	@echo "         pendb  --- pendulum by direct method "
	@echo "         testad --- test problem for simple index-2 examples "
	@echo "         testad2 -- test problem for well-poseness."
	@echo "   sensitivity by adjoint method with ADIFOR3.0: "
	@echo "         heatb2 --- Heat equation by direct method "
	@echo "         heatilub2- Heat equation by Krylov method with ILU  "
	@echo "                    preconditioner " 
#
# ----------------- Heat equation ( index-0 example )-------------------------
#
#    command : make heat; make heatb; make heatilu; make heatilub
#    test for 
#            sensitivity by forward method, 
#            sensitivity by adjoint method method, 
#
#   forward method
heat :  $(HEAT) 
	$(COMP) $(FFLAGS) -o heat $(HEAT) $(ADLIB)
heat2 :  $(HEAT2) 
	$(COMP) $(FFLAGS) -o heat2 $(HEAT2) $(ADLIB3)

#   adjoint method
heatb :  $(HEATB) 
	$(COMP) $(FFLAGS) -o heatb $(HEATB) $(ADLIB)
heatb2 : $(HEATB2) 
	-rm -f adlog*
	$(COMP) $(FFLAGS) -o heatb2 $(HEATB2) $(ADLIB3)
JM_Output/gh_resh.f JM_Output/g_resh.f : resh.f
	Adifor3.0 -spec spec_hf.ad3 resh.f
	Adifor3.0 -spec spec_hf2.ad3 resh.f
JM_Output/gh_reshdr.f : reshdr.f
	Adifor3.0 -spec spec_hfq.ad3 reshdr.f

JtM_Output/ady_resh.f JtM_Output/ady_accum.f : resh.f
	Adifor3.0 -spec spec_h.ad3 resh.f
JtM_Output/adp_resh.f JtM_Output/adp_accum.f : resh.f
	Adifor3.0 -spec spec_hp.ad3 resh.f
SJ_Output/sg_resh.f sg_accum.f : resh.f
	Adifor3.0 -spec spec_hfs.ad3 resh.f
#   support routines
res4sen/h_resh.f : resh.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=heat.cmp AD_IVARS=u,uprime,senpar \
               AD_PMAX=2
res4sen/h_reshdr.f : reshdr.f
	$(ADF) AD_SCRIPT=res4sensm_qres.adf AD_PROG=qres_heat.cmp \
                         AD_IVARS=u,uprime,senpar  AD_PMAX=2
res4sen/g_resh.f : resh.f
	$(ADF) AD_SCRIPT=res4senmv.adf AD_PROG=heat.cmp AD_IVARS=u,uprime,senpar
jacsp/i_resh.f : resh.f
	$(ADF) AD_SCRIPT=resjacsp.adf AD_PROG=heat.cmp AD_IVARS=u,uprime
resh_ady.f : resh.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		resh.f; mv resh_ad.f resh_ady.f
resh_adp.f : resh.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		resh.f; mv  resh_ad.f resh_adp.f

#   Krylov iterative method
heatilu: $(HEATILU)
	$(COMP) $(FFLAGS) -o heatilu $(HEATILU) $(ADLIB)  -lm
#   Krylov option for adjoint method
heatilub: $(HEATILUB)
	$(COMP) $(FFLAGS) -o heatilub $(HEATILUB) $(ADLIB3)  -lm
heatilub2: $(HEATILUB2)
	-rm -f adlog*
	$(COMP) $(FFLAGS) -o heatilub2 $(HEATILUB2) $(ADLIB3)  -lm
#   support routines
res4sen/h_reshilu.f res4sen/g_reshilu.f jacsp/i_reshilu.f : reshilu.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=heatilu.cmp \
                         AD_IVARS=u,uprime,senpar  AD_PMAX=2
	$(ADF) AD_SCRIPT=res4senmv.adf AD_PROG=heatilu.cmp \
                         AD_IVARS=u,uprime,senpar
	$(ADF) AD_SCRIPT=resjacsp.adf  AD_PROG=heatilu.cmp AD_IVARS=u,uprime
res4sen/h_reshiludr.f : reshiludr.f
	$(ADF) AD_SCRIPT=res4sensm_qres.adf AD_PROG=qres_heatilu.cmp \
                         AD_IVARS=u,uprime,senpar  AD_PMAX=2
JtM_Output/ady_reshilu.f : reshilu.f
	Adifor3.0 -spec spec_h.ad3 reshilu.f
JtM_Output/adp_reshilu.f : reshilu.f
	Adifor3.0 -spec spec_hp.ad3 reshilu.f
SJ_Output/sg_reshilu.f : reshilu.f
	Adifor3.0 -spec spec_hfs.ad3 reshilu.f

reshilu_ady.f : reshilu.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		reshilu.f; mv reshilu_ad.f reshilu_ady.f
reshilu_adp.f : reshilu.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		reshilu.f; mv  reshilu_ad.f reshilu_adp.f

hotb: $(HOTB)
	$(COMP) $(FFLAGS) -o hotb $(HOTB) $(ADLIB)  -lm
jacsp/i_reshot.f : reshot.f
	$(ADF) AD_SCRIPT=resjacsp.adf  AD_PROG=hot.cmp AD_IVARS=u,uprime
reshot_ady.f : reshot.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		reshot.f; mv reshot_ad.f reshot_ady.f
reshot_adp.f : reshot.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		reshot.f; mv  reshot_ad.f reshot_adp.f


#   Krylov option for adjoint method
burgerb: $(BGILUB)
	$(COMP) $(FFLAGS) -o burgerb $(BGILUB) $(ADLIB)  -lm
#   support routines
jacsp/i_res_burger.f : res_burger.f
	$(ADF) AD_SCRIPT=resjacsp.adf  AD_PROG=burger.cmp AD_IVARS=u,uprime
res_burger_ady.f : res_burger.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		res_burger.f; mv res_burger_ad.f res_burger_ady.f
res_burger_adp.f : res_burger.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		res_burger.f; mv  res_burger_ad.f res_burger_adp.f

#
# --------------------- Food-web problem (index-1 example) -------------------
#
#    command : make web; make webb; make webilu
#    test for 
#            initialization, 
#            sensitivity by forward method, 
#            sensitivity by adjoint method method, 
#
#     forward method
web :  $(WEB)
	$(COMP) $(FFLAGS) -o web $(WEB) $(ADLIB)  -lm 
#     adjoint method
webb : $(WEBB)
	$(COMP) $(FFLAGS) -o webb $(WEBB) $(ADLIB)
#     support routines for the adjoint method
resweb_ady.f : resweb.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		resweb.f; mv resweb_ad.f resweb_ady.f
resweb_adp.f : resweb.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		resweb.f; mv  resweb_ad.f resweb_adp.f
resweb_adyp.f : resweb.f
	$(TAMC) $(TAMCFLAGS) -input uprime -output delta -admark adyp -reverse \
		resweb.f; mv  resweb_ad.f resweb_adyp.f
#
#     support routines for the forward method
res4sen/h_resweb.f : resweb.f
	$(ADF) AD_SCRIPT=res4sensm_web.adf AD_PROG=web.cmp AD_IVARS=u,uprime,senpar \
               AD_PMAX=20
jac/k_resweb.f : resweb.f
	$(ADF) AD_SCRIPT=kmv.adf AD_PROG=web.cmp AD_IVARS=u,uprime
jacsp/i_resweb.f : resweb.f
	$(ADF) AD_SCRIPT=resjacsp.adf AD_PROG=web.cmp AD_IVARS=u,uprime
#
#     krylov method with ILU preconditioner
#
webilu : $(WEBILU)
	$(COMP) $(FFLAGS) -o webilu $(WEBILU) $(ADLIB)  -lm
#      support routines for the ILU preconditioner
res4sen/h_reswebilu.f res4sen/i_reswebilu.f : reswebilu.f
	$(ADF) AD_SCRIPT=res4sensm_web.adf AD_PROG=webilu.cmp \
                         AD_IVARS=u,uprime,senpar
	$(ADF) AD_SCRIPT=resjacsp.adf AD_PROG=webilu.cmp \
                         AD_IVARS=u,uprime
#
# ----------------------- wheel set problem (mixed index-1 and index-2) -------
#
#    command : make wheel
#    test for 
#            initialization, 
#            sensitivity by forward method, 
#
#     Only the forward sensitivity is computed
wheel :  $(WHEEL)
	$(COMP) $(FFLAGS) -o wheel $(WHEEL) $(ADLIB)  -lm
res4sen/h_wheel.f : wheel.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=wheel.cmp
jac/j_wheel.f : wheel.f
	$(ADF) AD_SCRIPT=resjac.adf AD_PROG=wheel.cmp AD_PMAX=17
#     support routine for the initialization if no sensitivity
jac/g_wheel.f : wheel.f
	$(ADF) AD_SCRIPT=res2nd_st.adf AD_PROG=wheel.cmp
res4sen/g_wheel.f : wheel.f
	$(ADF) AD_SCRIPT=res4senmv.adf AD_PROG=wheel.cmp
#     support routine for the initialization if no sensitivity
res4sen/g_g_wheel.f : res4sen/g_wheel.f
	$(ADF) AD_SCRIPT=res2nd.adf AD_PROG=res2nd_wheel.cmp
#
# ----------------------- CRANK problem (mixed index-1 and index-2) ------------
#    command : make crank
#    test for 
#            initialization, 
#            sensitivity by forward method, 
#
crank :  $(CRANK)
	$(COMP) $(FFLAGS) -o crank $(CRANK) $(ADLIB) -lm
res4sen/h_crank.f res4sen/h_ddot.f : crank.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=crank.cmp
jac/j_crank.f jac/j_ddot.f : crank.f
	$(ADF) AD_SCRIPT=resjac.adf AD_PROG=crank.cmp AD_PMAX=24
jac/g_crank.f jac/g_ddot.f : crank.f
	$(ADF) AD_SCRIPT=res2nd_st.adf AD_PROG=crank.cmp
res4sen/g_crank.f : crank.f
	$(ADF) AD_SCRIPT=res4senmv.adf AD_PROG=crank.cmp
res4sen/g_g_crank.f res4sen/g_g_ddot.f : res4sen/g_crank.f
	$(ADF) AD_SCRIPT=res2nd.adf AD_PROG=res2nd_crank.cmp
#
# ----------------------- Fekete problem (index-2 example) ----------------------
#
#    command : make fekete;  make feketeb
#    test for 
#            initialization, 
#            sensitivity by forward method, 
#            sensitivity by adjoint method
#
#      forward mode
fekete :  $(FEKETE)
	$(COMP) $(FFLAGS) -o fekete $(FEKETE) $(ADLIB) -lm
res4sen/h_fekete.f : fekete.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=fekete.cmp AD_PMAX=8
jac/j_fekete.f : fekete.f
	$(ADF) AD_SCRIPT=resjac.adf AD_PROG=fekete.cmp AD_PMAX=160
jac/g_fekete.f : fekete.f
	$(ADF) AD_SCRIPT=res2nd_st.adf AD_PROG=fekete.cmp
res4sen/g_fekete.f : fekete.f
	$(ADF) AD_SCRIPT=res4senmv.adf AD_PROG=fekete.cmp
res4sen/g_g_fekete.f : res4sen/g_fekete.f
	$(ADF) AD_SCRIPT=res2nd.adf AD_PROG=res2nd_fekete.cmp
res4sen/l_fekete.f : fekete.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=fekete.cmp \
               AD_TOP=init AD_IVARS=senpar AD_OVARS=y,yprime \
               AD_PREFIX=l AD_PMAX=8
#    adjoint mode
feketeb :  $(FEKETEB)
	$(COMP) $(FFLAGS) -o feketeb $(FEKETEB) $(ADLIB) -lm
fekete_ady.f : fekete.f
	$(TAMC) $(TAMCFLAGS) -input y -output delta -admark ady -reverse \
	fekete.f; mv fekete_ad.f fekete_ady.f
adj/g_fekete_ady.f adj/g_resad_fekete.f : fekete_ady.f
	$(ADF) AD_SCRIPT=index2ady.adf AD_PROG=index2ad_fekete.cmp
#    supportive routines for the adjoint mode
feketep : feketep.o jacsp/i_fekete.o fekete_ady.o fekete.o
	$(COMP) $(FFLAGS) -pg -o feketep feketep.o jacsp/i_fekete.o \
        fekete_ady.o  fekete.o resad_fekete.o $(ADLIB) -lm
jacsp/i_fekete.f : fekete.f
	$(ADF) AD_SCRIPT=resjacsp.adf AD_PROG=fekete.cmp
res4sen/h_qres_fekete.f : qres_fekete.f
	$(ADF) AD_SCRIPT=res4sensm_qres.adf AD_PROG=qres_fekete.cmp \
                          AD_IVARS=y,senpar AD_PMAX=168
#
# ----------------------- Pendulum problem (index-2 formulation) -------------
# 
#    command : make pend;  make pendb
#    test for 
#            initialization, 
#            sensitivity by forward method, 
#            sensitivity by adjoint method
#
pend :  $(PEND)
	$(COMP) $(FFLAGS) -o pend $(PEND) $(ADLIB)  -lm
res4sen/g_respend.f : 
	$(ADF) AD_SCRIPT=res4senmv.adf AD_PROG=pend.cmp
res4sen/g_g_respend.f : res4sen/g_respend.f
	$(ADF) AD_SCRIPT=res2nd.adf AD_PROG=res2nd_pend.cmp
pendb : $(PENDB)
	$(COMP) $(FFLAGS) -o pendb $(PENDB) $(ADLIB)  -lm
respend_ady.f : respend.f
	$(TAMC) $(TAMCFLAGS) -input y -output delta -admark ady -reverse \
	respend.f; mv respend_ad.f respend_ady.f
adj/g_respend_ady.f adj/g_resad_pend.f : respend_ady.f
	$(ADF) AD_SCRIPT=index2ady.adf AD_PROG=index2ad_pend.cmp
#
# ----------------------------- simple index-2 example -----------------------
#
#    command : make testad
#    test for 
#            sensitivity by forward method, 
#            sensitivity by adjoint method
#
testad : $(TESTAD)
	$(COMP) $(FFLAGS) -o testad $(TESTAD) $(ADLIB) -lm
res4sen/h_respend.f  jac/j_respend.f jacsp/i_respend.f : respend.f
	$(ADF) AD_SCRIPT=res4sensm.adf AD_PROG=pend.cmp
	$(ADF) AD_SCRIPT=resjac.adf  AD_PROG=pend.cmp AD_PMAX=5
	$(ADF) AD_SCRIPT=resjacsp.adf AD_PROG=pend.cmp
jac/g_respend.f : respend.f
	$(ADF) AD_SCRIPT=res2nd_st.adf AD_PROG=pend.cmp

#
# ----------------------------- simple index-2 example ------------------------
#
#    command : make testad2
#    test for 
#            well-poseness of the initialization problem
#
testad2 : $(TESTAD2)
	$(COMP) $(FFLAGS) -o testad2 $(TESTAD2) $(ADLIB) -lm

clean : 
	-rm -rf core *~

new :
	make clean
	-rm -f *.o heat heatilu web webilu heatb heatilub webb heatb2 heatilub2
	-rm -f pend pendb testad crank fekete feketeb wheel testad2
	-rm -rf jac jacsp res4sen adj AD_cache
	-rm -rf *_Output *_Cache adlog*
	-rm -f *_ady.f *_adyp.f *_adp.f

JtM_Output/ady_reshilu.o : JtM_Output/ady_reshilu.f
	$(ADFC) $(FFLAGS) -c -o $@ $<
JtM_Output/adp_reshilu.o : JtM_Output/adp_reshilu.f
	$(ADFC) $(FFLAGS) -c -o $@ $<
# Rule for compiling a Fortran source file:
.f.o: ;	$(COMP) $(FFLAGS) -c $*.f -o $*.o
.c.o:
	$(CC) $(FFLAGS) -c $*.c -o $*.o

