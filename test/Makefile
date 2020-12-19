include ../build/make.inc

ifeq ($(BUILDMPISMP),true)
EXE = ../bin/excitingmpismp
else
EXE = ../bin/excitingser
endif

TESTDIRS = test01 test02 test03 test04 test05 test06 test07 test08 test09 test10 test11 test12 test13 test14 test15

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
	cd ./test01/runAr; ../../$(EXE) input.xml; cp info.xml info-default.xml 2>/dev/null || :; cp atoms.xml atoms-default.xml 2>/dev/null || :
	cd ./test01/runAr; ../../$(EXE) input-nr.xml; cp info.xml info-nr.xml 2>/dev/null || :; cp atoms.xml atoms-nr.xml 2>/dev/null || :
	cd ./test01/runAr; ../../$(EXE) input-zora.xml; cp info.xml info-zora.xml 2>/dev/null || :; cp atoms.xml atoms-zora.xml 2>/dev/null || :
	cd ./test01/runAr; ../../$(EXE) input-iora.xml; cp info.xml info-iora.xml 2>/dev/null || :; cp atoms.xml atoms-iora.xml 2>/dev/null || :
	cd ./test01 ; bash assert.sh

test02::
	cd ./test02/runlapack; ../../$(EXE)
	cd ./test02/runarp; ../../$(EXE)
	cd ./test02; perl assert.pl

test03::
	cd test03/run && ../../$(EXE)
	cd test03/runmixer2 && ../../$(EXE)
	cd test03/runmixer3 && ../../$(EXE)
	cd test03 && perl report.pl

test04::
	cd test04/runso && ../../$(EXE)
	cd test04 && perl assert.pl

test05::
	cd test05/runlibxc && ../../$(EXE)
	cd test05 && bash assert.sh

test06::
	cd test06/runPbTiO3; ../../$(EXE) input.xml
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
	cd test11/run-pbe0 && ../../$(EXE)
	cd test11/run-hse && ../../$(EXE)
	cd test11 && perl assert.pl

test12::
	cd test12/runcollinear && ../../$(EXE)
	cd test12/runnoncollinear && ../../$(EXE)
	cd test12 && perl assert.pl

test13::
	cd test13/runldos && ../../$(EXE)
	cd test13 && perl assert.pl

test14::
	cd test14/runwannier && export MKL_NUM_THREADS=1 && ../../$(EXE)
	cd test14 && python assert.py
	
test15::
	cd test15/run && ../../$(EXE)
	cd test15 && perl assert.pl


cleantests::
	-rm -f  test*/*.xml
	-rm -f  test*/run*/*OUT
	-rm -f  test*/run*/*TXT
	-rm -f  test*/run*/dos*.xml
	-rm -f  test*/run*/*_scf.xml
	-rm -f  test*/run*/*_ini.xml
	-rm -f  test*/run*/bandstructure*.xml
	-rm -f  test*/run*/fermisurface.xml
	-rm -rf  test*/run*/LOSS
	-rm -rf  test*/run*/EPSILON
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
	-rm -f  test*/run*/spintext.xml
	-rm -f  directories.xml
	-rm -f  test*/run*/*.out
