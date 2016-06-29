include ../build/make.inc

ifeq ($(BUILDMPI),true)
EXE = ../bin/excitingmpi
else
EXE = ../bin/excitingser
endif

TESTDIRS = test01 test02 test03 test05 test06 test07 test08 test09 test10

test: 
	$(MAKE) cleantests
	$(MAKE) $(TESTDIRS)

summary:: test
	perl perl/makedlist.pl
	xsltproc report/collectTests.xsl directories.xml > report/reports.xml
	xsltproc report/reports2html.xsl report/reports.xml > report/report.html
	cp -f report/reports.xml report/oldreports.xml

onlysummary::
	perl perl/makedlist.pl
	xsltproc report/collectTests.xsl directories.xml > report/reports.xml
	xsltproc report/reports2html.xsl report/reports.xml > report/report.html
	cp -f report/reports.xml report/oldreports.xml

exciting::
	cd $(BINDIR) && $(MAKE) libs && $(MAKE)

test01::
	cd ./test01/runHe; ../../$(EXE)
	cd ./test01/runAr; ../../$(EXE) input.xml; cp info.xml info-default.xml; cp atoms.xml atoms-default.xml
	cd ./test01/runAr; ../../$(EXE) input-nr.xml; cp info.xml info-nr.xml; cp atoms.xml atoms-nr.xml
	cd ./test01/runAr; ../../$(EXE) input-zora.xml; cp info.xml info-zora.xml; cp atoms.xml atoms-zora.xml
	cd ./test01/runAr; ../../$(EXE) input-iora.xml; cp info.xml info-iora.xml; cp atoms.xml atoms-iora.xml
	cd ./test01 ; bash assert.sh 

test02::
	cd ./test02/runlapack; ../../$(EXE)
	cd ./test02/runarp; ../../$(EXE)
	cd ./test02/runarpnpacked; ../../$(EXE)
	cd ./test02; perl assert.pl

test03::  
	cd test03/run && ../../$(EXE)
	cd test03/runmixer2 && ../../$(EXE)
	cd test03/runmixer3 && ../../$(EXE)
	cd test03 && perl report.pl

test04::
	cd build; $(MAKE) test04
	cd ./test04/run ; ../../build/testprogram4

test05::  
	cd test05/runlibxc && ../../$(EXE)
	cd test05 && bash assert.sh

test06::
	cd test06/runPbTiO3; ../../$(EXE) input-nonsym.xml; cp info.xml info-nonsym.xml
	cd test06/runPbTiO3; ../../$(EXE) input-sym.xml; cp info.xml info-sym.xml      
	cd test06; bash assert.sh 

test07::	 
	cd test07/runBSE && ../../$(EXE)
	cd test07/runtddftBSE && ../../$(EXE)
	cd test07 && python assert.py

test08::
	cd test08/rungw && ../../$(EXE)
	cd test08 && perl assert.pl 

test09::
	cd test09 && bash validate.sh
	cd test09 && perl assert.pl  

test10::
	cd test10/run_properties && ../../$(EXE)
	cd test10 && perl assert.pl

test11::
	cd test11/run-pbe0-hf  && ../../$(EXE)
	cd test11/run-pbe0-oep && ../../$(EXE)
	cd test11 && perl assert.pl 

cleantests::
	-rm -f  test*/*.xml
	-rm -f  test*/run*/*OUT
	-rm -f  test*/run*/*TXT
	-rm -f  test*/run*/dos.xml
	-rm -f  test*/run*/*_scf.xml
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
	-rm -f  test*/run*/*3D.xml
	-rm -f  test*/run*/*.dat
	-rm -f  test*/run*/*.DAT
	-rm -f  test*/run*/*.xsf
	-rm -f  test*/run*/*.cube
