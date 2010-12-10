include ../build/make.inc

.IGNORE:test cleantests test01 test02 test03 test04 test05 testmpi

ifeq ($(BUILDMPI),true)
test: cleantests  test02 test03 test04 test05 testmpi
	perl perl/summary.pl
else
test: cleantests  test02 test03 test04 test05
	perl perl/summary.pl
endif	


test02::exciting
	-cd ./test02/runlapack ;../../../build/debug/exciting
	-cd ./test02/runarp; ../../../build/debug/exciting
	-cd ./test02/runarpnpacked;../../../build/debug/exciting
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

test05::exciting  
	-cd test05/runlibxc && ../../../build/debug/exciting
	-cd test05/runelectrstr && ../../../build/debug/exciting
	-cd test05 && bash assert.sh  

testmpi::
	-cd ../build/mpi/; $(MAKE) libs; $(MAKE) 
	-cd test02/runlapackmpi ;sh startmpirun 
	-cd test02 ; perl assertmpi.pl

cleantests::
	-rm -f  test*/*.xml
	-rm -f  test*/run*/*OUT
    -rm -f  test*/run*/dos.xml
    -rm -f  test*/run*/bandstructure.xml
    -rm -f  test*/run*/fermisurface.xml
    	



