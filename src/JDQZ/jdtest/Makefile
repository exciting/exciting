FC =gfortran 
FFLAGS = -u -O

LIBPATH = -L../jdlib -L./
LIBS = -ljdqz liblapack.a libblas.a -llapack -lblas
SUBDIRS = ../jdlib
OBJ = example.o

all: 	$(OBJ)
	$(FC) $(LDFLAGS) -o example $(OBJ) $(LIBPATH) $(LIBS)

clean:
	-rm -f example
	-rm -f `find . -name '*.[ao]' -print`
