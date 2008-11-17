include ../build/make.inc

.IGNORE:test cleantests test01 test02 test03 test04 testmpi

ifeq ($(BUILDMPI),true)
test: cleantests test01 test02 test03 test04 testmpi
	perl perl/summary.pl
endif 
ifeq ($(BUILDMPI),false)
test: cleantests test01 test02 test03 test04
	perl perl/summary.pl
endif 
	
test01::
	-cd build;$(MAKE) libs  test01
	-cd ./test01/run ; ../../build/testprogram1
	
test02::
	-cd build;$(MAKE) libs  test02
	-cd ./test02/runlapack ; ../../build/testprogram2
	-cd ./test02/runarp; ../../build/testprogram2
	-cd ./test02/runarpnpacked; ../../build/testprogram2
	-cd ./test01 ; perl assert.pl

test04::
	-cd build;$(MAKE) libs  test04
	-cd ./test04/run ; ../../build/testprogram4

exciting::
	cd ../build/debug/ && $(MAKE) libs && $(MAKE)

test03::exciting cleantests
	-cd test03/run && ../../../build/debug/exciting
	-cd test03/runmixer2 && ../../../build/debug/exciting
	-cd test03 && perl report.pl


testmpi::
	-cd ../build/mpi/; $(MAKE) libs; $(MAKE) 
	-cd test02/runlapackmpi ;sh startmpirun 
	-cd test02 ; perl assertmpi.pl

cleantests::
	-rm -f  test*/*.xml
	-rm -f  test*/run*/*OUT
