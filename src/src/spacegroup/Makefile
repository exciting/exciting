include ../../build/make.inc

bin:spacegroupbin

spacegroupbin:: ../../build/make.inc sginputmodule.f90
	../../build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p spacegroup \
	./ ../src_lib/r3*   \
	sginputmodule.f90   \
	../src_inputparser/modinputdom.f90   \
	&& $(MAKE) -f Makefile.mkmf spacegroup \
	&& cp spacegroup ../../bin/spacegroup

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)

doc::
	cd ../../docs/spacegroup/
	perl ../../build/utilities/scripts/genprotexdoc.pl ../../src/spacegroup/
	mv doc.pdf spacegroup.pdf

sginputmodule.f90:../../xml/schematofortran.xsl ../../xml/schema/symmetries.xsd
	xsltproc --stringparam tool "spacegroup" ../../xml/schematofortran.xsl ../../xml/schema/symmetries.xsd >sginputmodule.f90
