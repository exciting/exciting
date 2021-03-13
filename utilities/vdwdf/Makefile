#-------------------------------------------------------------------------------
# Compiler and flags
#-------------------------------------------------------------------------------

# Note CUBA compiles ifort, but it has not been tested
# Intel
#F90 = ifort
#CC = icc
#F90_OPTS = -O3
#CUBA_FFLAGS = -fPIC

# GCC
F90 = gfortran
CC = gcc
F90_OPTS = -O3 -march=native -ffree-line-length-0
CUBA_FFLAGS = 

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
vpath %.f90 src

SRC_modules = modules.f90

SRC_main = main.f90

SRC_routines = read_densities.f90 \
               read_gradients.f90 \
               read_phi.f90 \
               read_input.f90 \
               integrand.f90 \
               vdWDF_LL.f90 \
               vdWDF_VV.f90 \
               blend103.f90 \
               r3minv.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines)

OBJ = $(SRC:.f90=.o)
LIB_SYS = Cuba-3.2/libcuba.a -lm
EXE = noloco

vdW:cuba $(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_SYS)

cuba:
	cd Cuba-3.2; \
	make clean; \
	./configure --with-maxdim=6 --with-maxcomp=1 CC=$(CC) FFLAGS=$(CUBA_FFLAGS); \
	make lib

clean:
	rm -f *.o *.mod *~ fort.* ifc*

backup:
	tar -czf vdWDF.tgz $(SRC) Makefile make.inc


