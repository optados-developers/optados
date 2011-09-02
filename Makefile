#-------------------------------------------
# V A R I A B L E S   YO U   S H O U L D   
# S E T

# one of gfortran, g95, ifort, pfg90
SYSTEM := ifort 

# fast / debug
BUILD := fast 

# serial / mpi
COMMS_ARCH := serial 

# Where would you like the executables?
PREFIX:= ~/bin

SYSTEM := $(strip $(SYSTEM))
BUILD  := $(strip $(BUILD))
COMMS_ARCH := $(strip $(COMMS_ARCH))
PREFIX := $(strip $(PREFIX))


#-------------------------------------------
# T H I N G S   Y O U   M I G H T  W A N T  
# T O   F I D D L E   W I T H
ifeq ($(SYSTEM), g95)
   F90_SERIAL= g95
   F90_PARALLEL= openmpif90
   FFLAGS= -fendian=big
   FFLAGS_PARALLEL= -DMPI
   FFLAGS_FAST= -O3
   FFLAGS_DEBUG= -O0 -C -pg -g -Wall
   EXTENSION=.g95
endif

ifeq ($(SYSTEM), gfortran)
   F90_SERIAL= gfortran
   F90_PARALLEL= openmpif90 
   FFLAGS= -fconvert=big-endian
   FFLAGS_PARALLEL=  -DMPI
   FFLAGS_FAST= -O3
   FFLAGS_DEBUG= -O0 -pg -g 
   EXTENSION=.gfortran
endif

ifeq ($(SYSTEM), nag)
   F90_SERIAL= nagfor
   F90_PARALLEL= mpif90
   FFLAGS= -DNAG 
   FFLAGS_PARALLEL=
   FFLAGS_FAST= -O3 -Oassumed -w=all
   FFLAGS_DEBUG= -Ddebug -g -O0 -C=all -gline -w=all 
   EXTENSION=.nag
endif

ifeq ($(SYSTEM), ifort)
   F90_SERIAL= ifort
   F90_PARALLEL= mpifort 
   FFLAGS= -convert big_endian
   FFLAGS_PARALLEL= -DMPI
   FFLAGS_FAST= -O3
   FFLAGS_DEBUG= -O0 -C -pg -g -debug extended -traceback
   EXTENSION=.ifort
endif

ifeq ($(SYSTEM), pathscale)
   F90_SERIAL= pathf95
   F90_PARALLEL= mppathf90 
   FFLAGS= -byteswapio
   FFLAGS_PARALLEL=  -DMPI
   FFLAGS_FAST= -O3 -OPT:Ofast -ffast-math -OPT:recip=ON -OPT:malloc_algorithm=1
   FFLAGS_DEBUG=  -Ddebug -O0 -g -ffortran-bounds-check
   EXTENSION=.path
endif

ifeq ($(SYSTEM), pgf90)
   F90_SERIAL= pgf90
   F90_PARALLEL= mpif90
   FFLAGS= -Mbyteswapio
   FFLAGS_PARALLEL= -DMPI 
   FFLAGS_FAST= -O3
   FFLAGS_DEBUG= -O0 -C -pg -g -Mbounds 
   EXTENSION=.pgf90
endif

ifeq ($(SYSTEM), sun)
   F90_SERIAL= sunf95
   F90_PARALLEL= mpif90
   FFLAGS= 
   FFLAGS_PARALLEL= -DMPI
   FFLAGS_FAST= 
   FFLAGS_DEBUG= -O0 -C -pg -g
   EXTENSION=.sunf95
endif


#-------------------------------------------
# T H I N G S   Y O U   S H O U L D
# P R O B A B L Y   L E A V E   T O   T H E 
# D E V E L O P E R S

ifeq ($(BUILD), debug)
  FFLAGS+=$(FFLAGS_DEBUG)
  EXTENSION:=$(EXTENSION).debug
else
  FFLAGS+=$(FFLAGS_FAST)
endif


ifeq ($(COMMS_ARCH), mpi)
   FFLAGS+=$(FFLAGS_PARALLEL)
   F90=$(F90_PARALLEL)
   EXTENSION:=$(EXTENSION).mpi
else
   F90=$(F90_SERIAL)
endif


EXTENSION:=$(EXTENSION).x86_64


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

docs:
	make -C ./documents all

clean:
	rm -f *.o *.mod *.MOD
	make -C ./documents clean

dist:
	 @(tar cf - \
                ./*.?90 \
                ./examples/*/*.cell \
                ./examples/*/*.param \
                ./examples/*/*.odi \
                ./documents/*.tex \
                ./documents/*.pdf \
		./documents/Makefile \
		./documents/THINGS_TO_DO.txt \
                ./Makefile \
	| gzip -c > \
		./optados.tar.gz)

install: 
	cp optados$(EXTENSION) $(PREFIX)
