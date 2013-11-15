include ../../build/make.inc

F90_OPTS=$(F90_OPTS) -DSPECIES

GeneratedFiles=  ../../src/src_inputparser/inputmodules.f90 ../../src/src_inputparser/speciesmodules.f90	

species:: $(GeneratedFiles) fox
	../../build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p species \
	../src_lib ../mod_timing.F90 ../warning.f90 ../src_xs/getunit.F90 \
	../mod_muffin_tin.F90 ../mod_potential_and_density.F90 ../mod_SHT.F90 ../mod_Gvector.F90 \
	../src_xc ../mod_atoms.F90 ../atom.f90 ../genrmesh.f90 \
	../fftlib/cfftnd.f90 ../zfftifc.f90 \
	../src_inputparser/inputmodules.f90 \
	../src_inputparser/modinputdom.f90 \
	&& $(MAKE) -f Makefile.mkmf species \
	&& mv species ../../bin/species

fox:lib/libFoX_common.a

lib/libFoX_common.a:
	cd ../../src/FoX/ && \
	./configure FC=$(F90) && \
	make
	cp -r ../../src/FoX/objs/* ./

../../xml/excitinginput.xsd: ../../xml/schema/*.xsd  ../../xml/schema/schemaexpand.xsl
	xsltproc ../../xml/schema/schemaexpand.xsl ../../xml/schema/input.xsd > ../../xml/excitinginput.xsd

../../src/src_inputparser/inputmodules.f90: ../../xml/excitinginput.xsd ../../xml/schematofortran.xsl
	xsltproc ../../xml/schematofortran.xsl   ../../xml/excitinginput.xsd > ../../src/src_inputparser/inputmodules.f90

../../src/src_inputparser/speciesmodules.f90:   ../../xml/species.xsd ../../xml/schematofortran.xsl
	xsltproc  ../../xml/schematofortran.xsl  ../../xml/species.xsd > ../../src/src_inputparser/speciesmodules.f90

clean:
	rm -rf *.o *.mod *~ fort.* ifc* *.gcno lib/ finclude/ species 
