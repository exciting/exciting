test:test01 test02 test03

test01::
	cd build;make libs  test01

test02::
	cd build;make libs  test02

exciting::
	cd ../build/debug/ && make libs && make

test03::exciting
	cd test03/run && ../../../build/debug/exciting
	cd test03 && perl report.pl
