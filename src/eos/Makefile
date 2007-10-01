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
SRC = modeos.f90 eos.f90 amoeba.f90 amotry.f90 fopt.f90 eveos.f90 pveos.f90 \
      readdata.f90 getedata.f90 fitdata.f90 output.f90 constants.f90

OBJ = $(SRC:.f90=.o)
EXE = eos

eos:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ)

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.OUT gmon.out $(EXE)
