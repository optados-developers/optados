
include ../make.system



#-------------------------------------------
# T H I N G S   Y O U   S H O U L D
# P R O B A B L Y   L E A V E   T O   T H E 
# D E V E L O P E R S

ifeq ($(BUILD), debug)
  FFLAGS+=$(FFLAGS_DEBUG)
  EXTENSION:=.debug
else
  FFLAGS+=$(FFLAGS_FAST)
endif


ifeq ($(COMMS_ARCH), mpi)
   FFLAGS+=$(FFLAGS_PARALLEL)
   F90:=$(F90_PARALLEL)
   EXTENSION:=.mpi
else
   F90:=$(F90_SERIAL)
endif

EXTENSION:=$(EXTENSION)$(EXE_SUFFIX)

# Put object names here
OBJS=algorithms.o cell.o comms.o constants.o core.o dos.o dos_utils.o electronic.o io.o jdos.o jdos_utils.o optics.o parameters.o pdos.o xmgrace_utils.o

all : optados

optados : optados.f90 $(OBJS)
	$(F90) $(FFLAGS)   optados.f90 $(OBJS) -o optados$(EXTENSION) 

algorithms.o : algorithms.f90 io.o constants.o
	$(F90) -c $(FFLAGS) algorithms.f90

cell.o : cell.f90 comms.o constants.o io.o algorithms.o
	$(F90) -c $(FFLAGS) cell.f90

constants.o : constants.f90
	$(F90) -c $(FFLAGS) constants.f90

core.o : core.f90 parameters.o electronic.o dos_utils.o constants.o io.o xmgrace_utils.o
	$(F90) -c $(FFLAGS) core.f90

comms.o : comms.F90 constants.o 
	$(F90) -c $(FFLAGS) comms.F90

dos.o : dos.f90 dos_utils.o xmgrace_utils.o electronic.o parameters.o io.o
	$(F90) -c $(FFLAGS) dos.f90

dos_utils.o : dos_utils.f90 algorithms.o cell.o constants.o comms.o electronic.o io.o parameters.o
	$(F90) -c $(FFLAGS) dos_utils.f90

electronic.o : electronic.F90 comms.o parameters.o constants.o 
	$(F90) -c $(FFLAGS) electronic.F90

io.o : io.F90 constants.o comms.o
	$(F90) -c $(FFLAGS) io.F90

jdos.o : jdos.f90 jdos_utils.o xmgrace_utils.o electronic.o parameters.o io.o
	$(F90) -c $(FFLAGS) jdos.f90

jdos_utils.o : jdos_utils.f90 algorithms.o cell.o constants.o comms.o electronic.o io.o parameters.o dos_utils.o
	$(F90) -c $(FFLAGS) jdos_utils.f90

optics.o : optics.f90 constants.o io.o electronic.o jdos_utils.o cell.o parameters.o xmgrace_utils.o
	$(F90) -c $(FFLAGS) optics.f90

parameters.o : parameters.f90  cell.o constants.o io.o comms.o
	$(F90) -c $(FFLAGS) parameters.f90

pdos.o : pdos.F90  cell.o comms.o constants.o io.o dos.o electronic.o
	$(F90) -c $(FFLAGS) pdos.F90

xmgrace_utils.o : xmgrace_utils.f90 io.o constants.o
	$(F90) -c $(FFLAGS) xmgrace_utils.f90


# Utility targets
.PHONY: clean veryclean


clean:
	rm -f *.o *.mod *.MOD


veryclean:
	make  clean
	rm -f optados$(EXTENSION) optados$(EXTENSION).debug optados$(EXTENSION).mpi

install: 
	cp optados$(EXTENSION) $(INSTALL_DIR)
