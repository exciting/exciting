include ../../build/make.inc

species::fox
	../../build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p species \
	../src_lib ../mod_timing.F90 ../warning.f90 \
	../mod_muffin_tin.F90 ../mod_potential_and_density.F90 ../mod_SHT.F90 ../mod_Gvector.F90 \
	../src_xc ../mod_atoms.F90 ../atom.f90 ../genrmesh.f90 \
	../fftlib/cfftnd.f90 ../zfftifc.f90 \
	../src_inputparser/inputmodules.f90 \
	../src_inputparser/modinputdom.f90 \
	&& $(MAKE) -f Makefile.mkmf species \
	&& mv species ../../bin/species

fox:libFoX_common.a

libFoX_common.a:
	cd ../../src/FoX/ && \
	./configure FC=$(F90) && \
	make
	cp -r ../../src/FoX/objs/* ./

	
clean:
	rm -rf *.o *.mod *~ fort.* ifc* *.gcno lib/ finclude/ species 
