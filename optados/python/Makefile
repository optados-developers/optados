MODULE_NAME=OptaPyDOS

F90_FILES=../src/algorithms.f90 ../src/cell.f90 ../src/comms.F90 ../src/constants.f90 ../src/core.f90 ../src/dos.f90 ../src/dos_utils.f90 ../src/electronic.f90 ../src/io.f90 ../src/jdos.f90 ../src/jdos_utils.f90 ../src/optics.f90 ../src/parameters.f90 ../src/pdos.F90 ../src/xmgrace_utils.f90

KIND_MAP= kind_map

python:
	f90wrap -vvvvv -m $(MODULE_NAME) $(F90_FILES)  -k $(KIND_MAP) 
	f2py -c -m _OptaPyDOS -I"../src" ../src/*.o f90wrap*.f90

# Utility targets
.PHONY: clean veryclean

clean:
	rm f90wrap*

veryclean: clean
	rm OptaPyDOS.py  _OptaPyDOS.so
