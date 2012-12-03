include ../../build/make.inc

bin:stateconvertbin

stateconvertbin:: ../../build/make.inc 
	../../build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p stateconvert \
	./ ../portstate.F90   \
	 ../ioarray.F90   \
	 ../mod_misc.F90   \
	 ../versions_gt.F90   \
	&& $(MAKE) -f Makefile.mkmf stateconvert \
	&& cp stateconvert ../../bin/stateconvert

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)

doc::
	cd ../../docs/stateconvert/
	perl ../../build/utilities/scripts/genprotexdoc.pl ../../src/stateconvert/
	mv doc.pdf stateconvert.pdf
