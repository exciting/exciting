all:test01 test02 test04

include ../../build/Make.common

Makefile.test1::
	perl ../../utilities/versionstamp.pl
	../../utilities/mkmf   -t template -f -m Makefile.test1 -p testprogram1  \
	$(SOURCEDIRS) ../test01/src/ \
	../src/testframework/


Makefile.test2::
	perl ../../utilities/versionstamp.pl
	../../utilities/mkmf -t template -f -m Makefile.test2 -p testprogram2  \
	$(SOURCEDIRS) ../test02/src \
	../src/testframework


Makefile.test4::
	perl ../../utilities/versionstamp.pl
	 ../../utilities/mkmf -t template -f -m Makefile.test4 -p testprogram4  \
	$(SOURCEDIRS) ../test04/src \
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
	

#for the rules to build the libraries include the other makefile



