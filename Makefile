# Makefile to compile tdci
# /wsu/el7/pgi/2018-187/linux86-64/18.7/lib/

SHELL := /bin/bash

FC = pgf95
FFLAGS = -tp x64 -O3 
LIBDIR = /wsu/el7/pgi/2018-187/linux86-64/18.7/lib/
LAPACK_LIBS = -L${LIBDIR} -llapack -lblas
OPT_FLAGS = -Minfo -Mneginfo -time -fast -Mconcur=allcores -mp=allcores -Munroll -Mvect
ALL_OBJECTS = tdci.o units.o setup_variables.o control_variables.o variables_global.o write_info.o read_integrals.o initialize.o getfield.o util.o getham0_cisd.o getham0.o getham.o propagate.o Zpropagate.o analysis.o davidson_ip.o qcmatrixio.o


all : update tdci 

tdci : $(ALL_OBJECTS)
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(ALL_OBJECTS) $(LAPACK_LIBS) -o tdci

tdci.o : tdci.f90 variables_global.o initialize.o write_info.o getfield.o util.o getham0_cisd.o getham0.o getham.o analysis.o propagate.o Zpropagate.o davidson_ip.o	
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c tdci.f90 -o tdci.o

#davidson_cisd.o : mod_davidson_cisd.f90 variables_global.o util.o getham0_cisd.o
#	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_davidson_cisd.f90 -o davidson_cisd.o

davidson_ip.o : mod_davidson_ip.f90 variables_global.o util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_davidson_ip.f90 -o davidson_ip.o

Zpropagate.o : mod_Zpropagate.f90 analysis.o util.o variables_global.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_Zpropagate.f90 -o Zpropagate.o

propagate.o : mod_propagate.f90 analysis.o util.o variables_global.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_propagate.f90 -o propagate.o

analysis.o : mod_analysis.f90 util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_analysis.f90 -o analysis.o

getham.o : mod_getham.f90 getham0.o getham0_cisd.o variables_global.o util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_getham.f90 -o getham.o 

getham0.o : mod_getham0.f90 variables_global.o read_integrals.o util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_getham0.f90 -o getham0.o

getham0_cisd.o : mod_getham0_cisd.f90 variables_global.o read_integrals.o util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_getham0_cisd.f90 -o getham0_cisd.o

getfield.o : mod_getfield.f90 variables_global.o util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_getfield.f90 -o getfield.o

initialize.o : mod_initialize.f90 read_integrals.o variables_global.o util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_initialize.f90 -o initialize.o

read_integrals.o : mod_readintegrals.f90 variables_global.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_readintegrals.f90 -o read_integrals.o

write_info.o : mod_write.f90 units.o setup_variables.o control_variables.o variables_global.o util.o
	$(FC) -c mod_write.f90 -o write_info.o

util.o : mod_util.f90 read_integrals.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) -c mod_util.f90 -o util.o

variables_global.o : mod_variables_global.f90 units.o setup_variables.o control_variables.o 
	$(FC) -c mod_variables_global.f90 -o variables_global.o

control_variables.o : mod_variables_control.f90
	$(FC) -c mod_variables_control.f90 -o control_variables.o

setup_variables.o : mod_variables_setup.f90 
	$(FC) -c mod_variables_setup.f90 -o setup_variables.o

units.o : mod_variables_units.f90
	$(FC) -c mod_variables_units.f90 -o units.o

update : 
	./CHANGE_DATE

qcmatrixio.o : qcmatrixio.F
	pgfortran -c qcmatrixio.F -o qcmatrixio.o

clean :
	@( echo " *** CLEANING *** ")
	rm $(ALL_OBJECTS)



