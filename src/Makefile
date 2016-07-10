####################################################################
#
#	Makefile to compile, print or check the LMG model
#       programs
#
#       by Currix TM.
#
#
BINPATH = ../bin
####################################################################
####################################################################
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain U2-U1
#
AVEC_Lipkin_U1_SRC = nrtype.f90 defparam_lipkin.f90 avalavec_lipkin_U1.f90 
#
#
# ALGEBRAIC MODEL ENERGIES AND EIGENVECTORS (MODEL HAMILTONIAN) Chain U2-SO2
#
AVEC_Lipkin_SO2_SRC = nrtype.f90 defparam_lipkin.f90 avalavec_lipkin_SO2.f90 
#
####################################################################
####################################################################
FC = gfortran
#
FOPT   = -Wall
MODLIB   = -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
FOPT   += $(MODLIB)
#
LAPACK95 = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95
## Cluster
LAPACK77 = -llapack
BLAS     = -lblas
#
LIBS     = $(LAPACK95)  $(LAPACK77) $(BLAS)
#
######################################################################
######################################################################
#
all: avec_lipkin_gfort_u1 avec_lipkin_gfort_so2
#
.PHONY : all
#
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain U2-U1)
avec_lipkin_gfort_u1: $(AVEC_Lipkin_U1_SRC) Makefile
	$(FC) -o $(BINPATH)/$@ $(FOPT) $(AVEC_Lipkin_U1_SRC) $(LIBS)
#################################################################
# HAMILTONIAN EIGENVALUES AND EIGENVECTORS (MODEL HAMILTONIAN, Chain U2-SO2)
avec_lipkin_gfort_so2: $(AVEC_Lipkin_SO2_SRC) Makefile
	$(FC) -o $(BINPATH)/$@ $(FOPT) $(AVEC_Lipkin_SO2_SRC) $(LIBS)
#################################################################
clean:
	@echo "Cleaning object files"
	@rm -f *.o *.mod $(BINPATH)/avec_lipkin_gfort_u1 $(BINPATH)/avec_lipkin_gfort_so2




