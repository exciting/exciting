
AR = ar

include ../../build/make.inc

#-------------------------------------------------------------------------------
# Suffix Rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .F .f90
.F.o:
	$(F77) $(F77_OPTS) -c $<
	
.f90.o:
	$(F90) $(F90_OPTS) -c $<
#-------------------------------------------------------------------------------
# File dependencies
#-------------------------------------------------------------------------------
SRC77 = Lebedev-Laikov.F
SRC90 = leblaik.f90
OBJ = $(SRC77:.F=.o) $(SRC90:.f90=.o)

leblaiklib:	$(OBJ)
	$(AR) -rc leblaiklib.a $(OBJ)

clean:
	rm -f *.o *.mod *~ *.a ifc* *.gcno

