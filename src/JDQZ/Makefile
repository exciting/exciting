
FFLAGS = -u -O
include ../../make.inc
FC = $(F77)
LIBPATH = -L../jdlib
LIBS = -ljdqr -llapack -lblas
SUBDIRS = jdlib
TEST    = jdtest

all: 	libs 

jdqr: 	
	@for i in $(TEST); do (cd $$i; $(MAKE)); done

libs:
	@for i in $(SUBDIRS); do (cd $$i; $(MAKE)); done

clean:
	-rm -f jdqz
	-rm -f `find . -name '*.[ao]' -print`
