include ../../build/make.inc
bin:speciesbin

speciesbin:: ../../build/make.inc 
	../../utilities/mkmf -t ./template -f -m Makefile.mkmf -p species \
	./ ../src_lib ../src_xc ../modxcifc.f90 ../atom.f90 \
	../zfftifc.f90 \
	../mod_muffin_tin.F90 \
	../mod_atoms.F90 \
	../mod_Gvector.F90 \
    ../mod_SHT.F90 \
    ../mod_potential_and_density.F90 \
    ../src_inputparser/inputmodules.f90 \
    ../src_inputparser/modinputdom.f90 \
	&& $(MAKE) -f Makefile.mkmf species \
	&& cp species ../../bin/species
	
clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)


#include ../../build/Make.common
