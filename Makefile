



default: build/make.inc all

all: serial mpi  smp mpiandsmp  eos spacegroup species

build/make.inc:
	perl ./setup.pl
	
include build/make.inc

serial:
	cd build/serial; $(MAKE) 

mpi:
	cd build/mpi; $(MAKE) 

smp:
	cd build/smp; $(MAKE)
	
debug:
	cd build/debug; $(MAKE)

mpiandsmp:
	cd build/mpiandsmp; $(MAKE)

test::
	cd test/; $(MAKE) -i

doc: excitingdoc spacegroupdoc
	
excitingdoc::
	$(MAKE) -f build/Make.common doc
	
spacegroupdoc::
	cd src/spacegroup; $(MAKE) doc
 


eos::
	cd src/eos; $(MAKE)
	
spacegroup::
	cd src/spacegroup; $(MAKE)
	
species::
	cd src/species; $(MAKE)

debian:all doc
	cd debian &&  sh makepackage.sh

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
	rm docs/exciting/*
	rm docs/spacegroup/*
	rm -r debian/debian/usr

