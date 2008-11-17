all:test

Makefile.test1::
	perl ../../utilities/versionstamp.pl
	../../utilities/mkmf   -t template -f -m Makefile.test1 -p testprogram1  \
	../../src  \
	../src \
	../../src/src_xs \
	../../src/src_sym \
	../../src/src_fermisurfdx \
	../../src/src_iterative_solver \
	../../src/src_iterative_solver/src_diis \
	../../src/src_mixing \
	../src/testplans/test1.f90 \
	../src/testframework/


Makefile.test2::
	perl ../../utilities/versionstamp.pl
	../../utilities/mkmf -t template -f -m Makefile.test2 -p testprogram2  \
	../../src  \
	../src \
	../../src/src_xs \
	../../src/src_sym \
	../../src/src_fermisurfdx \
	../../src/src_iterative_solver \
	../../src/src_iterative_solver/src_diis \
	../../src/src_mixing \
	../test02/src \
	../src/testframework


Makefile.test4::
	perl ../../utilities/versionstamp.pl
	 ../../utilities/mkmf -t template -f -m Makefile.test4 -p testprogram4  \
	../../src  \
	../src \
	../../src/src_xs \
	../../src/src_sym \
        ../../src/src_fermisurfdx \
        ../../src/src_iterative_solver \
	../../src/src_iterative_solver/src_diis \
	../../src/src_mixing \
	../test04/src \
	../src/testframework




testprogram1::Makefile.test1
	make -f Makefile.test1 testprogram1 

testprogram2::Makefile.test2
	make -f Makefile.test2 testprogram2 
	
testprogram4::Makefile.test4
	make -f Makefile.test4 testprogram4
 
	
test01::testprogram1
	

test02::testprogram2


test04::testprogram4
	

#for the rules to build the libraries include the other makefile


include ../../build/Make.common
