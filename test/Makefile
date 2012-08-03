include ../build/make.inc


ifeq ($(TestMPI),true)
TESTDIRS =test02 test03 test04 test05 test07 test08 test09  testmpi
else
TESTDIRS =test02 test03 test04 test05 test07 test08 test09
endif	

test: 
	$(MAKE) cleantests 
	$(MAKE) $(TESTDIRS)

summary:: test
	perl perl/makedlist.pl
	xsltproc report/collectTests.xsl directories.xml >report/reports.xml
	xsltproc report/reports2html.xsl report/reports.xml >report/report.html
	cp -f report/reports.xml report/oldreports.xml

exciting::
	cd ../build/serial/ && $(MAKE) libs && $(MAKE)
	cd ../build/debug/ && $(MAKE) libs && $(MAKE)

	
test02::exciting
	cd ./test02/runlapack ;../../../build/debug/exciting
	cd ./test02/runarp; ../../../build/debug/exciting
	cd ./test02/runarpnpacked;../../../build/debug/exciting
	cd ./test02 ; perl assert.pl

test03::exciting  
	cd test03/run && ../../../build/debug/exciting
	cd test03/runmixer2 && ../../../build/debug/exciting
	cd test03/runmixer3 && ../../../build/debug/exciting
	cd test03 && perl report.pl

test04::
	cd build;$(MAKE) libs  test04
	cd ./test04/run ; ../../build/testprogram4


test05::exciting  
	cd test05/runlibxc && ../../../build/debug/exciting
	cd test05/runelectrstr && ../../../build/debug/exciting
	cd test05 && bash assert.sh  

test07a::exciting 	
	cd test07/runBSE && ../../../build/serial/exciting
test07b::exciting 
	cd test07/runtddft && ../../../build/serial/exciting
	
test07:: test07a test07b	 
	cd test07 && bash assert.sh  

test08::exciting  
	cd test08/rungw && ../../../build/debug/exciting
	cd test08 && perl assert.pl 
	 
test09::
	cd test09 && bash validate.sh
	cd test09 && perl assert.pl  



testmpi::
	cd ../build/mpi/; $(MAKE) libs; $(MAKE) 
	cd test02/runlapackmpi ;sh startmpirun test08
	cd test02 ; perl assertmpi.pl
	
cleantests::
	-rm -f  test*/*.xml
	-rm -f  test*/run*/*OUT
	-rm -f  test*/run*/dos.xml
	-rm -f  test*/run*/bandstructure.xml
	-rm -f  test*/run*/fermisurface.xml
	-rm -f  test*/run*/LOSS*
	-rm -f  test*/run*/EPSILON*


