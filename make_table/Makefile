BASE = $(shell /bin/pwd)

#defaults for external variables                                                                                   
COMPILER ?= GCC
#COMPILER ?= INTEL
EXECTBL ?= gen_tbl

SRCS = legndrpol.f90 al_int.f90 gen_main.f90



#HDF5_DIR=${HDF5_HOME}
HDF5_DIR=/afs/mpa-garching.mpg.de/data/haakoan/phd/libs/hdf5_gcc
INC_PATH =-I./wigxjpf-1.6/inc -I./wigxjpf-1.6/mod -I./fastwigxj-1.0/inc -I./fastwigxj-1.0/mod
LIB_PATH =-L./wigxjpf-1.6/lib -L./fastwigxj-1.0/lib
LIBS =  -lfastwigxj -lwigxjpf_quadmath -lwigxjpf  -lm


ifeq ($(COMPILER),INTEL)
HDF5_FC = ifort
F90 = h5fc
F90FLAGS = -O3 -convert big_endian 
endif

ifeq ($(COMPILER),GCC)
HDF5_FC = gfortran
F90 = gfortran
F90FLAGS = -cpp -O3  -pedantic -fimplicit-none -w # -fcheck=all -fbacktrace
endif



main: 
	$(F90) -o $(EXECTBL) $(SRCS) $(F90FLAGS)  $(LIB_PATH) $(INC_PATH) $(LIBS)



clean:
	rm legndrpol.o sph_leg_pol.mod rotator reader.o

#-fconvert=big-endian
