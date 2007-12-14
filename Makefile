
MAKE = make
serial:
	cd build/serial; $(MAKE) libs
	cd build/serial; $(MAKE) 
parallel:
	cd build/parallel; $(MAKE) libs
	cd build/parallel; $(MAKE) 

all:serial parallel
	cp build/serial/make.inc ./
	cd src/eos; $(MAKE)
	cd src/spacegroup; $(MAKE)
	cd src/species; $(MAKE)

clean:

	cd build/serial; $(MAKE) clean
	cd build/parallel; $(MAKE) clean
	cd src/eos; $(MAKE) clean
	cd src/spacegroup; $(MAKE) clean
	cd src/species; $(MAKE) clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*

