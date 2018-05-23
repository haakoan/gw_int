[21~BASE = $(shell /bin/pwd)

#defaults for external variables                                                                                   
COMPILER ?= INTEL
#COMPILER ?= INTEL
EXECTBL ?= cheek_t
SRCS = ./src/fortran/datatypes.f90 ./src/fortran/dataread.f90  ./src/fortran/legndrpol.f90 ./src/fortran/tblinf.f90 ./src/fortran/quads.f90 ./src/fortran/sphexp.f90 ./src/fortran/integraltest.f90

SRCS1 = ./src/fortran/datatypes.f90 ./src/fortran/dataread.f90  ./src/fortran/legndrpol.f90 ./src/fortran/tblinf.f90 ./src/fortran/quads.f90 ./src/fortran/sphexp.f90 ./src/fortran/integraltest.f90


HDF5_DIR=${HDF5_HOME}
HDF5_DIR=/afs/mpa-garching.mpg.de/data/haakoan/phd/libs/hdf5_gcc
INC_PATH =  -I./src/fortran/wigxjpf-1.6/inc -I./src/fortran/wigxjpf-1.6/mod -I./src/fortran/fastwigxj-1.0/inc -I./src/fortran/fastwigxj-1.0/mod #-I$(HDF5_DIR)/include/ 
LIB_PATH = -L./src/fortran/wigxjpf-1.6/lib -L./src/fortran/fastwigxj-1.0/lib #-L$(HDF5_DIR)/lib
#LIBS = -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5
LIBS =   -lm -I/usr/lib64/gfortran/modules -I/usr/include -L/usr/lib64 -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lhdf5hl_fortran 


ifeq ($(COMPILER),INTEL)
HDF5_FC = ifort
F90 = h5fc
F90FLAGS = -cpp -lstdc++ -O3 -implicitnone -assume realloc-lhs -xHost -heap-arrays #-openmp
#F90FLAGS = -cpp -lstdc++ -g -debug all -check bounds -fimplicit-none -w  -check all -traceback -implicitnone -ftrapuv -check-uninit
endif

ifeq ($(COMPILER),GCC)
HDF5_FC = gfortran
F90 = gfortran
F90FLAGS = -cpp  -lstdc++ -pedantic -fimplicit-none -w  -fcheck=all -fbacktrace -O3  #-fopenmp
endif



main: 
	icc -c ./src/cpp/lookup.cpp -O3 -std=c++11 
	mv lookup.o ./objs
	h5fc -o $(EXECTBL) ./objs/lookup.o $(SRCS) $(F90FLAGS) $(INC_PATH) $(LIB_PATH) $(LIBS)
	mv *.mod ./objs	
	mv $(EXECTBL) ./bin
	@echo "Compiled by: $$(whoami)"
	@echo "Created on: $$(date)"

clean:
	rm legndrpol.o sph_leg_pol.mod rotator reader.o

#-fconvert=big-endian
