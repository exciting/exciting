include ../build/make.inc

ifeq ($(BUILDMPI),true)
test:test01 test02 test03 test04 testmpi
endif 
ifeq ($(BUILDMPI),false)
test:test01 test02 test03 test04
endif 
	perl perl/summary.pl
test01::
	cd build;make libs  test01

test02::
	cd build;make libs  test02

test04::
	cd build;make libs  test04


exciting::
	cd ../build/debug/ && make libs && make

test03::exciting
	cd test03/run && ../../../build/debug/exciting
	cd test03/runmixer2 && ../../../build/debug/exciting
	cd test03 && perl report.pl


testmpi::
	cd ../build/mpi/; make libs; make 
	cd test02/runlapackmpi ;sh startmpirun 
	cd test02 ; perl assertmpi.pl


