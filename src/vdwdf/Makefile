include make.inc

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC_modules = modules.f90

SRC_main = main.f90

SRC_routines = read_xsf.f90 read_phi.f90 read_input.f90 \
               integrand.f90 \
               vdWDF_LL.f90 vdWDF_VV.f90 \
               blend103.f90 r3minv.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines)

OBJ = $(SRC:.f90=.o)
EXE = vdWDF.x

vdW:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_SYS)


clean:
	rm -f *.o *.mod *~ fort.* ifc*

backup:
	tar -czf vdWDF.tgz $(SRC) Makefile make.inc


