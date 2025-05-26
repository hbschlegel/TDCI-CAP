# Makefile to compile tdci

SHELL := /bin/bash

# Directories
SRC := src
OBJ := obj
MOD := $(OBJ)/modules
DEP := $(OBJ)/dep
BIN := bin

HDF5_DIR := dep/hdf5-1.14.4-3/hdf5-install
HDF5_INC  := -I$(HDF5_DIR)/include
HDF5_LIB  := -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5

# Source and object files
F90_SOURCES := $(wildcard $(SRC)/*.f90)
F_SOURCES   := $(wildcard $(SRC)/*.F)
SOURCES     := $(F90_SOURCES) $(F_SOURCES)

F90_OBJECTS := $(patsubst $(SRC)/%.f90,$(OBJ)/%.o,$(F90_SOURCES))
F_OBJECTS   := $(patsubst $(SRC)/%.F,$(OBJ)/%.o,$(F_SOURCES))
OBJECTS     := $(F90_OBJECTS) $(F_OBJECTS)

# If the working directory 'dirty' (does not exactly match the commit), append '+'
dirty := $(shell if [ -z "$$(git status --porcelain --untracked-files=no 2>/dev/null)" ]; then echo ""; else echo "+"; fi)
# Include git version number and dirty indicator
GIT_HASH := $(shell git rev-parse --short HEAD)$(dirty)

# Compiler and flags
FC          := nvfortran
FFLAGS      := -tp=host -O3 -g -Mpreprocess -DGIT_HASH=\"$(GIT_HASH)\"
OPT_FLAGS   := -Minfo -Mneginfo -time -fast -Mconcur=allcores -mp=allcores -Munroll -Mvect
O_FLAGS     := -module $(MOD) -I$(MOD) $(HDF5_INC)
LIBS        := -llapack -lblas $(HDF5_LIB)
DEPFLAGS    := -MMD -MF $(DEP)/$*.d

DEPLIST := $(DEP)/deps.mk



# Default target
.PHONY : all
all : update dep $(DEPLIST) $(BIN)/tdci $(BIN)/tdci_core

.PHONY : dep
dep :
	@echo "==> running ./configure to install HDF5 and makedepf90"
	@./configure
	@echo "==> dependencies ready"


$(BIN)/tdci : $(SRC)/tdci.sh
	cp $(SRC)/tdci.sh $@
	chmod +x $@

$(BIN)/tdci_core : $(OBJECTS)
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(OBJECTS) $(LIBS) $(O_FLAGS) -o $@

$(DEP)/deps.mk : $(SOURCES) | $(DEP)
	@echo '  [makedepf90] scanning sources …'
	@dep/makedepf90/bin/makedepf90 -b $(OBJ)/ -I$(SRC) -I$(MOD) $(HDF5_INC) $(SOURCES) > $@

# Pattern rules
$(OBJ)/%.o : $(SRC)/%.f90 | $(DEP)
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(O_FLAGS) -c $< -o $@

$(OBJ)/%.o : $(SRC)/%.F | $(DEP)
	$(FC) $(FFLAGS) $(OPT_FLAGS) $(O_FLAGS) -c $< -o $@

$(MOD)/%.mod : $(OBJ)/%.o ;

# Directory setup
update : | $(BIN) $(OBJ) $(MOD) $(DEP)

$(BIN) $(OBJ) $(MOD) $(DEP) :
	@mkdir -p $@

# Include generated dependency files
-include $(DEPLIST)

.PHONY : clean
clean :
	@echo '*** CLEANING ***'
	@rm -rf $(OBJ) $(BIN)



