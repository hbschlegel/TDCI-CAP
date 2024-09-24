# Makefile to compile tdci
# /wsu/el7/pgi/2018-187/linux86-64/18.7/lib/

SHELL := /bin/bash


# Directories
# Source, object, binary, and module directories
SRC = src
OBJ = obj
BIN = bin
MOD = $(OBJ)/modules

# Source and object files
F90_SOURCES = $(wildcard $(SRC)/*.f90)
F_SOURCES = $(wildcard $(SRC)/*.F)
SOURCES = $(F90_SOURCES) $(F_SOURCES)

F90_OBJECTS = $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(F90_SOURCES))
F_OBJECTS = $(patsubst $(SRC)/%.F,$(OBJ)/%.o,$(F_SOURCES))
OBJECTS = $(F90_OBJECTS) $(F_OBJECTS)

$(info SOURCES: $(SOURCES))
$(info OBJECTS: $(OBJECTS))

#FC = pgf95
FC = nvfortran
FFLAGS = -tp px -O3 -g
LAPACK_LIBS = -llapack -lblas
OPT_FLAGS = -Minfo -Mneginfo -time -fast -Mconcur=allcores -mp=allcores -Munroll -Mvect
O_FLAGS = -module $(MOD) -I$(MOD)

#ALL_OBJECTS = $(OBJ)/tdci.o $(OBJ)/variables_units.o $(OBJ)/variables_setup.o $(OBJ)/variables_control.o $(OBJ)/variables_global.o $(OBJ)/write_info.o $(OBJ)/readintegrals.o $(OBJ)/initialize.o $(OBJ)/getfield.o $(OBJ)/util.o $(OBJ)/sort.o $(OBJ)/io_binary.o $(OBJ)/getham0_cisd.o $(OBJ)/getham0.o $(OBJ)/getham.o $(OBJ)/propagate.o $(OBJ)/Zpropagate.o $(OBJ)/analysis.o $(OBJ)/davidson_ip.o $(OBJ)/qcmatrixio.o



all : update $(BIN)/tdci 

$(BIN)/tdci : $(OBJECTS)
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(OBJECTS) $(LAPACK_LIBS) $(O_FLAGS) -o $(BIN)/tdci

# General pattern rule for object files
$(OBJ)/%.o : $(SRC)/%.f90
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

# Require all other objects be built
$(OBJ)/tdci.o : $(SRC)/tdci.f90 $(filter-out $(OBJ)/tdci.o,$(OBJECTS))
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

#$(OBJ)/davidson_cisd.o : $(SRC)/davidson_cisd.f90 $(OBJ)/variables_global.o $(OBJ)/util.o $(OBJ)/getham0_cisd.o
#	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $(SRC)/davidson_cisd.f90 -o $(OBJ)/davidson_cisd.o

$(OBJ)/davidson_ip.o : $(SRC)/davidson_ip.f90 $(OBJ)/variables_global.o $(OBJ)/util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/Zpropagate.o : $(SRC)/Zpropagate.f90 $(OBJ)/analysis.o $(OBJ)/util.o $(OBJ)/variables_global.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/propagate.o : $(SRC)/propagate.f90 $(OBJ)/analysis.o $(OBJ)/util.o $(OBJ)/variables_global.o $(OBJ)/sort.o $(OBJ)/io_binary.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/analysis.o : $(SRC)/analysis.f90 $(OBJ)/util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/getham.o : $(SRC)/getham.f90 $(OBJ)/getham0.o $(OBJ)/getham0_cisd.o $(OBJ)/variables_global.o $(OBJ)/util.o $(OBJ)/sort.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/getham0.o : $(SRC)/getham0.f90 $(OBJ)/variables_global.o $(OBJ)/readintegrals.o $(OBJ)/util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/getham0_cisd.o : $(SRC)/getham0_cisd.f90 $(OBJ)/variables_global.o $(OBJ)/readintegrals.o $(OBJ)/util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/getfield.o : $(SRC)/getfield.f90 $(OBJ)/variables_global.o $(OBJ)/util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/initialize.o : $(SRC)/initialize.f90 $(OBJ)/readintegrals.o $(OBJ)/variables_global.o $(OBJ)/util.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/readintegrals.o : $(SRC)/readintegrals.f90 $(OBJ)/io_binary.o $(OBJ)/variables_global.o $(OBJ)/qcmatrixio.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/write_info.o : $(SRC)/write_info.f90 $(OBJ)/variables_units.o $(OBJ)/variables_setup.o $(OBJ)/variables_control.o $(OBJ)/variables_global.o $(OBJ)/util.o
	$(FC) $(O_FLAGS) -c $< -o $@ 

$(OBJ)/util.o : $(SRC)/util.f90 $(OBJ)/readintegrals.o
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(LAPACK_LIBS) $(O_FLAGS) -c $< -o $@

$(OBJ)/variables_global.o : $(SRC)/variables_global.f90 $(OBJ)/variables_units.o $(OBJ)/variables_setup.o $(OBJ)/variables_control.o 
	$(FC) $(O_FLAGS) -c $< -o $@ 

$(OBJ)/variables_control.o : $(SRC)/variables_control.f90
	$(FC) $(O_FLAGS) -c $< -o $@ 

$(OBJ)/variables_setup.o : $(SRC)/variables_setup.f90 
	$(FC) $(O_FLAGS) -c $< -o $@ 

$(OBJ)/variables_units.o : $(SRC)/variables_units.f90
	$(FC) $(O_FLAGS) -c $< -o $@ 

$(OBJ)/sort.o : $(SRC)/sort.f90
	$(FC) $(O_FLAGS) -c $< -o $@ 

$(OBJ)/io_binary.o : $(SRC)/io_binary.f90
	$(FC) $(O_FLAGS) -c $< -o $@ 

update : 
	mkdir -p $(BIN)
	mkdir -p $(OBJ)
	mkdir -p $(MOD)
	./CHANGE_DATE

$(OBJ)/qcmatrixio.o : $(SRC)/qcmatrixio.F
	pgfortran $(O_FLAGS) -c $< -o $@

clean :
	@( echo " *** CLEANING *** ")
	rm -r obj/
	rm -r bin/
	mkdir obj/
	mkdir bin/



