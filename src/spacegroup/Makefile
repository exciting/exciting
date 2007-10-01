include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

SRC_modules = modmain.f90

SRC_main = main.f90

SRC_routines = readinput.f90 gencrystal.f90 sgsymb.f90 seitzgen.f90 \
 gengroup.f90 seitzeq.f90 seitzmul.f90 writegeom.f90 writexsf.f90

SRC_lib = findprim.f90 r3frac.f90 r3taxi.f90 r3mv.f90 r3cross.f90 r3minv.f90 \
 r3mm.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines) $(SRC_lib)

OBJ = $(SRC:.f90=.o)
EXE = spacegroup

spacegroup:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ)

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.OUT *.xsf gmon.out sgg.* $(EXE)

copy:
	cd ..;cp $(SRC_lib) spacegroup

doc:
	../protex -s $(SRC_main) $(SRC_modules) $(SRC_routines) $(SRC_lib) \
         > sgg.tex
	latex sgg;latex sgg;latex sgg
	dvipdf sgg.dvi

