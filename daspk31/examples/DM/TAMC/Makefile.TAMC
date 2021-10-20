#
# This makefile compiles and loads the DDASPK example program. It was originally
# written for an SGI machine, but has been updated to run on Linux and Cygwin.
# If necessary, change the constants COMP and FFLAGS below for the
# compiler to be used.

.SUFFIXES: .F .f .o .c

# Warning: f95 is not a viable altnernative to g77. At least for the brp
# example, it will compile and run, but it will produce garbage.
COMP = g77
CC   = gcc
# The -Df2cFortran flag is for Cygwin
CFLAGS = -g -O2 -Df2cFortran 
FFLAGS = -g -O2

# PATH to TAMC preprocessor command
TAMC = tamc
TAMCFLAGS = -module res -pure

SOLVR = ../../../solver
ADSLVR = ../../../ad_solver
ADSLVR_OBJS = $(ADSLVR)/ad3setup_dummy.o
PRECON = ../../../preconds

#----------------------------------------------------------------
# DASPK without involving the ADIFOR and MPI packages, but including
# calls to the Tapenade automatic differentiation tool.
#
OBJS = $(SOLVR)/ddaspk.o $(SOLVR)/daux.o $(SOLVR)/dlinpk.o \
       $(SOLVR)/mpi_dummy.o $(SOLVR)/adf_dummy.o

OBJSB = $(ADSLVR)/daspk_adjoint.o $(ADSLVR)/adwrite.o $(ADSLVR_OBJS)

SENSD = $(SOLVR)/dsensd.o

HEATB = dheat_ad.o resh.o $(OBJS) $(PRECON)/dbanpre.o $(OBJSB) resh_ady.o resh_adp.o

HEATILUB = dheatilu_ad.o reshilu.o $(OBJS) $(PRECON)/dilupre.o $(PRECON)/ad_dilupre.o\
          $(PRECON)/dsparsk.o $(OBJSB) reshilu_ady.o reshilu_adp.o

WEBB = dweb_ad.o $(OBJS) $(PRECON)/drbdpre.o $(PRECON)/drbgpre.o $(OBJSB) \
       resweb.o resweb_ady.o resweb_adp.o resweb_adyp.o

# -----------------------------default --------------------------------------
#
default : 
	@echo "Usage: make OBJ "
	@echo " "
	@echo " where OBJ can be: "
	@echo "   sensitivity by adjoint method: "
	@echo "         heatb  --- Heat equation by direct method "
	@echo "         heatilub - Heat equation by Krylov method with ILU  "
	@echo "                    preconditioner " 
	@echo "         webb   --- food web problem by direct method "

#
# test with sensitivity to the derived information

# heatb and support routines for the adjoint method		
heatb:	$(HEATB)
	$(COMP) $(FFLAGS) -o heatb $(HEATB)
resh_ady.f : resh.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		resh.f; mv resh_ad.f resh_ady.f
resh_adp.f : resh.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		resh.f; mv  resh_ad.f resh_adp.f

# heatilub and support routines for the adjoint method		
heatilub: $(HEATILUB) 
	$(COMP) $(FFLAGS) -o heatilub $(HEATILUB) 
reshilu_ady.f : reshilu.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		reshilu.f; mv reshilu_ad.f reshilu_ady.f
reshilu_adp.f : reshilu.f
	$(TAMC) $(TAMCFLAGS) -input senpar -output delta -admark adp -reverse \
		reshilu.f; mv reshilu_ad.f reshilu_adp.f

# webb and support routines for the adjoint method		
webb :  $(WEBB) 
	$(COMP) $(FFLAGS) -o webb $(WEBB)
resweb_ady.f : resweb.f
	$(TAMC) $(TAMCFLAGS) -input u -output delta -admark ady -reverse \
		resweb.f; mv resweb_ad.f resweb_ady.f
resweb_adp.f : resweb.f
	$(TAMC) $(TAMCFLAGS) -input senpara -output delta -admark adp -reverse \
		resweb.f; mv  resweb_ad.f resweb_adp.f
resweb_adyp.f : resweb.f
	$(TAMC) $(TAMCFLAGS) -input uprime -output delta -admark adyp -reverse \
		resweb.f; mv  resweb_ad.f resweb_adyp.f

clean : 
	-rm -rf core core.* *~ *.stackdump
	
cleaner :
	make clean
	-rm -rf $(SOLVR)/*.o $(ADSLVR)/*.o $(PRECON)/*.o
	-rm -rf *.o heatb heatilub webb heatb.exe heatilub.exe webb.exe

new :
	make cleaner
	-rm -rf fort.*
#	-rm -f *_ady.f *_adp.f

# Rule for compiling a Fortran source file:
.f.o: ;	$(COMP) $(FFLAGS) -c $*.f -o $*.o
.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o
