 
include ../../make.inc
 	
#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .F90
.F90.o:
	$(F90) $(F90_OPTS) -c -I../ $<
#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC = iterativeseceqnfv.F90 residualvector.F90 projectedsecequn.F90 \
  setupprojectedhamilton.F90 calcupdatevector.F90 updateevecfv.F90

OBJ = $(SRC:.F90=.o)
all:libiterative_solver.a
	

libiterative_solver.a: $(OBJ)
	rm libiterative_solver.a; ar -v -q libiterative_solver.a  $(OBJ)

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.OUT gmon.out *.a $(EXE)

