all:
	cd src; make all
	cd src/eos; make
	cd src/spacegroup; make
	cd src/species; make

clean:
	cd src; make cleanall
	cd src/eos; make clean
	cd src/spacegroup; make clean
	cd src/species; make clean
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.exe exdg.*

