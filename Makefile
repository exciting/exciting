
MAKE = make

include make.inc

all:
	cd src; $(MAKE) all
	cd src/eos; $(MAKE)
	cd src/spacegroup; $(MAKE)
	cd src/species; $(MAKE)

clean:
	cd src; $(MAKE) cleanall
	cd src/eos; $(MAKE) clean
	cd src/spacegroup; $(MAKE) clean
	cd src/species; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*

