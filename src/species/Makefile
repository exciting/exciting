include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix Rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC_modules = modxcifc.f90

SRC_routines = species.f90

SRC_lib = \
 atom.f90 rdirac.f90 rdiracdme.f90 rdiracint.f90 fderiv.f90 spline.f90 \
 polynom.f90 radmesh.f90

SRC_xc = \
  xc_pzca.f90 xc_pwca.f90 xc_pbe.f90 xc_am05.f90 xc_xalpha.f90 xc_wc06.f90 \
  x_wc06.f90 x_pbe.f90 c_pbe.f90 c_pbe_gcor.f90

SRC = $(SRC_modules) $(SRC_routines) $(SRC_lib) $(SRC_xc)
OBJ = $(SRC:.f90=.o)
EXE = species

species:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_SYS)

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)

copy:
	cd ..;cp $(SRC_lib) $(SRC_xc) modxcifc.f90 species