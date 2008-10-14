all:test

Makefile.test1::
	sh ../../utilities/versionstamp.sh
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
	sh ../../utilities/versionstamp.sh
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


Makefile.test2mpi::
	sh ../../utilities/versionstamp.sh
	../../utilities/mkmf -t templatempi -f -c -DMPI-DMPISEC-DMPIRHO\
	 -m Makefile.test2mpi -p testprogram2mpi \
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

testprogram1::Makefile.test1
	make -f Makefile.test1 testprogram1 

testprogram2::Makefile.test2
	make -f Makefile.test2 testprogram2 
	
testprogram2mpi::Makefile.test2mpi
	make -f Makefile.test2mpi testprogram2mpi 
	
test01::testprogram1
	cd ../test01/run ; ../../build/testprogram1

test02::testprogram2
	cd ../test02/runlapack ; ../../build/testprogram2
	cd ../test02/runarp; ../../build/testprogram2
	cd ../test02 ; perl assert.pl

testmpi::
	cd ../../build/mpi/; make libs exciting
	cd ../test02/runlapackmpi ;sh startmpirun 
	cd ../test02 ; perl assertmpi.pl

libs:libbzint.a	libblas.a liblapack.a fftlib.a  libarpack.a 

libblas.a: 
	cp -f ../../build/make.inc ../..
	cd ../../src/BLAS && make clean && make
	cp ../../src/BLAS/libblas.a .
	cd ../../src/BLAS && make clean

liblapack.a: 
	cp -f ../../build/make.inc ../..
	cd ../../src/LAPACK && make clean && make
	cp ../../src/LAPACK/liblapack.a .
	cd ../../src/LAPACK && make clean

fftlib.a: 
	cp -f ../../build/make.inc ../..
	cd ../../src/fftlib && make clean && make
	cp ../../src/fftlib/fftlib.a .
	cd ../../src/fftlib && make clean

libbzint.a: 
	cp -f ../../build/make.inc ../..
	cd ../../src/libbzint && make clean && make
	cp ../../src/libbzint/libbzint.a .
	cd ../../src/libbzint && make clean


libarpack.a: 
	cp -f ../../build/make.inc ../..
	cd ../../src/ARPACK && $(MAKE) clean && make lib
	cp ../../src/ARPACK/libarpack.a .
	cd ../../src/ARPACK && $(MAKE) clean



clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a

### build libbzint, do not clean source tree of library and build exciting
bzex:
	cp -f ../../build/make.inc ../..
	cd ../../src/src_libbzint && make
	cp ../../src/src_libbzint/libbzint.a .
	rm -f exciting
	make serial
