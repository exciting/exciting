



default: build/make.inc all
 
build/make.inc:
	perl ./setup.pl
	
include build/make.inc

serial:
	cd build/serial; $(MAKE) libs
	cd build/serial; $(MAKE) 

mpi:
	cd build/mpi; $(MAKE) libs
	cd build/mpi; $(MAKE) 

smp:
	cd build/smp; $(MAKE) libs
	cd build/smp; $(MAKE)
debug:
	cd build/debug; $(MAKE) libs
	cd build/debug; $(MAKE)

mpiandsmp:
	cd build/mpiandsmp; $(MAKE) libs
	cd build/mpiandsmp; $(MAKE)

test::
	cd test/; $(MAKE) -i
	
doc::
	$(MAKE) -f build/Make.common doc
 
all:serial mpi  smp mpiandsmp doc
	cp build/make.inc ./
	cd src/eos; $(MAKE)
	cd src/spacegroup; $(MAKE)
	cd src/species; $(MAKE)

clean:

	cd build/serial; $(MAKE) clean cleanlibs
	cd build/mpi; $(MAKE) clean cleanlibs
	cd build/smp; $(MAKE) clean cleanlibs
	cd build/debug; $(MAKE) clean cleanlibs
	cd build/mpiandsmp; $(MAKE) clean cleanlibs
	cd test/build ;$(MAKE) clean cleanlibs
	cd src/eos; $(MAKE) clean
	cd src/spacegroup; $(MAKE) clean
	cd src/species; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*
	rm bin/exciting?*
	rm interfaces/*

