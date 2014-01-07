
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

optados : optados.f90 $(OBJS) ../make.system
	$(F90) $(FFLAGS) $< $(OBJS) -o optados$(EXTENSION) 

algorithms.o : algorithms.f90 io.o constants.o ../make.system
	$(F90) -c $(FFLAGS) $<

cell.o : cell.f90 comms.o constants.o io.o algorithms.o ../make.system
	$(F90) -c $(FFLAGS) $<

constants.o : constants.f90 ../make.system
	$(F90) -c $(FFLAGS) $<

core.o : core.f90 parameters.o electronic.o dos_utils.o constants.o io.o xmgrace_utils.o ../make.system
	$(F90) -c $(FFLAGS) $<

comms.o : comms.F90 constants.o  ../make.system
	$(F90) -c $(FFLAGS) $<

dos.o : dos.f90 dos_utils.o xmgrace_utils.o electronic.o parameters.o io.o ../make.system
	$(F90) -c $(FFLAGS) $<

dos_utils.o : dos_utils.f90 algorithms.o cell.o constants.o comms.o electronic.o io.o parameters.o ../make.system
	$(F90) -c $(FFLAGS) $<

electronic.o : electronic.f90 comms.o parameters.o constants.o  ../make.system
	$(F90) -c $(FFLAGS) $<

io.o : io.f90 constants.o comms.o ../make.system
	$(F90) -c $(FFLAGS) $<

jdos.o : jdos.f90 jdos_utils.o xmgrace_utils.o electronic.o parameters.o io.o ../make.system
	$(F90) -c $(FFLAGS) $<

jdos_utils.o : jdos_utils.f90 algorithms.o cell.o constants.o comms.o electronic.o io.o parameters.o dos_utils.o ../make.system
	$(F90) -c $(FFLAGS) $<

optics.o : optics.f90 constants.o io.o electronic.o jdos_utils.o cell.o parameters.o xmgrace_utils.o ../make.system
	$(F90) -c $(FFLAGS) $<

parameters.o : parameters.f90  cell.o constants.o io.o comms.o ../make.system
	$(F90) -c $(FFLAGS) $<

pdos.o : pdos.F90  cell.o comms.o constants.o io.o dos.o electronic.o ../make.system
	$(F90) -c $(FFLAGS) $<

xmgrace_utils.o : xmgrace_utils.f90 io.o constants.o ../make.system
	$(F90) -c $(FFLAGS) $<


# Utility targets
.PHONY: clean veryclean


clean:
	rm -f *.o *.mod *.MOD


veryclean: clean
	rm -f optados$(EXTENSION) optados$(EXTENSION).debug optados$(EXTENSION).mpi

install: 
	cp optados$(EXTENSION) $(INSTALL_DIR)
