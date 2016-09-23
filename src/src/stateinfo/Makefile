include ../../build/make.inc

bin:stateinfobin

stateinfobin:: ../../build/make.inc 
	../../build/utilities/mkmf -t ./template -f -m Makefile.mkmf -p stateinfo \
	./ ../portstate.F90   \
	 ../ioarray.F90   \
	 ../mod_misc.F90   \
	 ../versions_gt.F90   \
	&& $(MAKE) -f Makefile.mkmf stateinfo \
	&& cp stateinfo ../../bin/stateinfo

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out $(EXE)

doc::
	cd ../../docs/stateinfo/
	perl ../../build/utilities/scripts/genprotexdoc.pl ../../src/stateinfo/
	mv doc.pdf stateinfo.pdf
