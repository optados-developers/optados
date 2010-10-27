DEBUG ?= 0

ifeq ($(DEBUG), 1)
 FFLAGS_OPT=-O0 -fbounds-check -pg -g -ff2c #-fno-underscoring
 NAME_OPT=.debug
else
 FFLAGS_OPT=-O3 -funroll-loops
 NAME_OPT=
endif

F90=gfortran

#FFLAGS_OPT=-O3 -funroll-loops 
#FFLAGS_OPT=-O0 -fbounds-check -pg -g
#FFLAGS_OPT=-O0 -fbounds-check -pg -g -ff2c #-fno-underscoring
#FFLAGS_OPT=-O0  -fprofile-arcs -ftest-coverage

FFLAGS= $(FFLAGS_OPT) -fconvert=big-endian -frecord-marker=4

# Where would you like the executables?
BIN_DIR=./
#${HOME}/bin
EXTENSION=$(NAME_OPT).x86_64


# Put object names here
OBJS=algorithms.o comms.o parameters.o io.o cell.o constants.o dos.o

all : optados

optados : $(OBJS)
	$(F90) $(FFLAGS) optados.f90 $(OBJS) -o $(BIN_DIR)/optdos$(EXTENSION) 

algorithms.o : algorithms.f90 io.o constants.o
	$(F90) -c $(FFLAGS) algorithms.f90

cell.o : cell.f90 constants.o
	$(F90) -c $(FFLAGS) cell.f90

constants.o : constants.f90
	$(F90) -c $(FFLAGS) constants.f90

dos.o : dos.f90 constants.o
	$(F90) -c $(FFLAGS) dos.f90

comms.o : comms.F90 constants.o
	$(F90) -c $(FFLAGS) comms.F90

io.o : io.F90 constants.o
	$(F90) -c $(FFLAGS) io.F90

parameters.o : parameters.f90 algorithms.o constants.o io.o
	$(F90) -c $(FFLAGS) parameters.f90


# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD


