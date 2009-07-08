include ../../build/make.inc

bin:spacegroupbin

spacegroupbin:: ../../build/make.inc 
	../../utilities/mkmf -t ./template -f -m Makefile.mkmf -p spacegroup \
	./ ../src_lib/r3*   \
	 ../src_inputparser/inputmodules.f90    \
	 ../src_inputparser/modinputdom.f90   \
	&& $(MAKE) -f Makefile.mkmf spacegroup \
	&& cp spacegroup ../../bin/spacegroup
	
clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)
	
doc::
	cd ../../docs/spacegroup/
	perl ../../utilities/scripts/genprotexdoc.pl ../../src/spacegroup/
	mv doc.pdf spacegroup.pdf