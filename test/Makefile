include ../build/make.inc

.IGNORE:test cleantests test01 test02 test03 test04 testmpi

ifeq ($(BUILDMPI),true)
test: cleantests  test02 test03 test04 testmpi
	perl perl/summary.pl
else
test: cleantests  test02 test03 test04
	perl perl/summary.pl
endif	


test02::
	-cd build;$(MAKE) libs  test02
	-cd ./test02/runlapack ; ../../build/testprogram2
	-cd ./test02/runarp; ../../build/testprogram2
	-cd ./test02/runarpnpacked; ../../build/testprogram2
	-cd ./test02 ; perl assert.pl

test04::
	-cd build;$(MAKE) libs  test04
	-cd ./test04/run ; ../../build/testprogram4

exciting::
	cd ../build/debug/ && $(MAKE) libs && $(MAKE)

test03::exciting cleantests
	-cd test03/run && ../../../build/debug/exciting
	-cd test03/runmixer2 && ../../../build/debug/exciting
	-cd test03/runmixer3 && ../../../build/debug/exciting
	-cd test03 && perl report.pl

test05::exciting cleantests
	-cd test05/runconverge && ../../../build/debug/exciting
	-cp test05/runconverge/STATE.OUT   test05/runrestart
	-cd test05/runrestart && ../../../build/debug/exciting

test06::exciting
	-cd build;$(MAKE)   test06
	-cd ./test06/run ; ../../build/testxml
	
testxmlin::exciting	
	-cd ./test06/run2 ; ../../../build/debug/exciting
	
testxmlind::exciting
	-cd ./test06/run2 ; gdb ../../../build/debug/exciting
	
testxmlAld::exciting 
	-cd ./test06/runAl ; gdb -x debug ../../../build/debug/exciting
	
testspin::exciting 
	-cd ./test06/spinpol ; gdb -x debug ../../../build/debug/exciting

testrelax::exciting 
	-cd ./test06/struct ; gdb -x debug ../../../build/debug/exciting


testmpi::
	-cd ../build/mpi/; $(MAKE) libs; $(MAKE) 
	-cd test02/runlapackmpi ;sh startmpirun 
	-cd test02 ; perl assertmpi.pl

cleantests::
	-rm -f  test*/*.xml
	-rm -f  test*/run*/*OUT



