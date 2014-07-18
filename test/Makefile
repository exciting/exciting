include ../build/make.inc


ifeq ($(TestMPI),true)
TESTDIRS =test01 test02 test03 test04 test05 test06 test07 test08 test09  testmpi
else
TESTDIRS =test01 test02 test03 test04 test05 test06 test07 test08 test09
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

	
test01::exciting
	cd ./test01/runHe ;../../../build/debug/exciting
	cd ./test01/runAr; ../../../build/debug/exciting input.xml      ; cp info.xml info-default.xml ; cp atoms.xml atoms-default.xml
	cd ./test01/runAr; ../../../build/debug/exciting input-nr.xml   ; cp info.xml info-nr.xml      ; cp atoms.xml atoms-nr.xml
	cd ./test01/runAr; ../../../build/debug/exciting input-zora.xml ; cp info.xml info-zora.xml    ; cp atoms.xml atoms-zora.xml
	cd ./test01/runAr; ../../../build/debug/exciting input-iora.xml ; cp info.xml info-iora.xml    ; cp atoms.xml atoms-iora.xml
	cd ./test01 ; bash assert.sh 

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

test06::exciting
	cd ./test06/runPbTiO3; ../../../build/debug/exciting input-nonsym.xml      ; cp info.xml info-nonsym.xml
	cd ./test06/runPbTiO3; ../../../build/debug/exciting input-sym.xml  ; cp info.xml info-sym.xml      
	cd ./test06 ; bash assert.sh 

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
	-rm -f  test*/run*/atoms*
	-rm -f  test*/run*/info*xml
	-rm -f  test*/run*/geometry*xml
	-rm -f  test*/run*/FERMISURF*
	-rm -f  test*/run*/RHO2D*
	-rm -f  test*/run*/err*
	-rm -f  test*/err
	-rm -f  test*/run*/result*
	-rm -f  test*/result*
