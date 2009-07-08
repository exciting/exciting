include ../../build/make.inc
bin:speciesbin

speciesbin:: ../../build/make.inc 
	../../utilities/mkmf -t ./template -f -m Makefile.mkmf -p species \
	./ ../src_lib ../src_xc ../modxcifc.f90 ../atom.f90 \
	../mod_muffin_tin.F90 \
	&& $(MAKE) -f Makefile.mkmf species \
	&& cp species ../../bin/species
	
clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)


#include ../../build/Make.common
