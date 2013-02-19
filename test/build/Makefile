all:test01 test02 test04
include ../../build/make.inc
include ../../build/Make.common

Makefile.test1::$(GeneraedFiles)
	perl ../../build/utilities/versionstamp.pl
	../../build/utilities/mkmf   -t template -f -m Makefile.test1 -p testprogram1  \
	$(SOURCEDIRS) $(LIBXCMOD) ../test01/src/ \
	../src/testframework/


Makefile.test2::$(GeneraedFiles)
	perl ../../build/utilities/versionstamp.pl
	../../build/utilities/mkmf -t template -f -m Makefile.test2 -p testprogram2  \
	$(SOURCEDIRS) $(LIBXCMOD) ../test02/src \
	../src/testframework


Makefile.test4::$(GeneraedFiles)
	perl ../../build/utilities/versionstamp.pl
	 ../../build/utilities/mkmf -t template -f -m Makefile.test4 -p testprogram4  \
	$(SOURCEDIRS) $(LIBXCMOD) ../test04/src \
	../src/testframework


testprogram1:Makefile.test1 libs
	make -f Makefile.test1 testprogram1 

testprogram2:Makefile.test2 libs
	make -f Makefile.test2 testprogram2 
	
testprogram4:Makefile.test4 libs
	make -f Makefile.test4 testprogram4
 
	
test01::testprogram1
	

test02::testprogram2


test04::testprogram4
	
test06:testxml
	
	
testxml:Makefile.testxml libs 
	make -f Makefile.testxml testxml
	
	
	
Makefile.testxml::  ../../src/src_inputparser/inputmodules.f90 ../../src/src_inputparser/speciesmodules.f90 
	perl ../../build/utilities/versionstamp.pl
	../../build/utilities/mkmf  -t template -f -m Makefile.testxml -p testxml  \
	../../src/src_inputparser/ \
	$(SOURCEDIRS) ../test06/src/ \

	
#for the rules to build the libraries include the other makefile


clean:
	rm -f test* 
