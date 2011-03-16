DEBUG ?= 1
G95 ?= 1
IFORT ?= 0


ifeq ($(DEBUG), 1)
 FFLAGS_OPT=-O0 -C -pg -g  #-fno-underscoring
 NAME_OPT=.debug
else
 FFLAGS_OPT=-O3 -DMPI -convert big_endian 
 NAME_OPT=
endif

F90=mpif90

#FFLAGS_OPT=-O3 -funroll-loops 
#FFLAGS_OPT=-O0 -fbounds-check -pg -g
#FFLAGS_OPT=-O0 -fbounds-check -pg -g -ff2c #-fno-underscoring
#FFLAGS_OPT=-O0  -fprofile-arcs -ftest-coverage

FFLAGS= $(FFLAGS_OPT) 

ifeq ($(G95), 1)
 F90=g95
 FFLAGS= 
 NAME_OPT=.g95
endif

# Where would you like the executables?
BIN_DIR=./
#${HOME}/bin
EXTENSION=$(NAME_OPT).x86_64


# Put object names here
OBJS=algorithms.o cell.o comms.o constants.o core.o dos.o dos_utils.o electronic.o io.o jdos.o jdos_utils.o optics.o parameters.o pdos.o

all : optados

optados : optados.f90 $(OBJS)
	$(F90) $(FFLAGS)  -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -Wl,--end-group optados.f90 $(OBJS) -o $(BIN_DIR)/optdos$(EXTENSION) 

algorithms.o : algorithms.f90 io.o constants.o
	$(F90) -c $(FFLAGS) algorithms.f90

cell.o : cell.f90 comms.o constants.o io.o algorithms.o
	$(F90) -c $(FFLAGS) cell.f90

constants.o : constants.f90
	$(F90) -c $(FFLAGS) constants.f90

core.o : core.f90 constants.o io.o
	$(F90) -c $(FFLAGS) core.f90

comms.o : comms.F90 constants.o io.o
	$(F90) -c $(FFLAGS) comms.F90

dos.o : dos.f90 dos_utils.o
	$(F90) -c $(FFLAGS) dos.f90

dos_utils.o : dos_utils.f90 algorithms.o cell.o constants.o comms.o electronic.o io.o parameters.o
	$(F90) -c $(FFLAGS) dos_utils.f90

electronic.o : electronic.F90 comms.o constants.o 
	$(F90) -c $(FFLAGS) electronic.F90

io.o : io.F90 constants.o
	$(F90) -c $(FFLAGS) io.F90

jdos.o : jdos.f90 jdos_utils.o
	$(F90) -c $(FFLAGS) jdos.f90
	
jdos_utils.o : jdos_utils.f90 algorithms.o cell.o constants.o comms.o electronic.o io.o parameters.o
	$(F90) -c $(FFLAGS) jdos_utils.f90

optics.o : optics.f90 constants.o io.o
	$(F90) -c $(FFLAGS) optics.f90

parameters.o : parameters.f90  cell.o constants.o io.o
	$(F90) -c $(FFLAGS) parameters.f90

pdos.o : pdos.F90  cell.o comms.o constants.o io.o dos.o electronic.o
	$(F90) -c $(FFLAGS) pdos.F90


# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD


