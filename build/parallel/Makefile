include make.inc


parallel: make.inc
	../../utilities/mkmf -t template -f -m Makefile.parallel -p \
		excitingmpi \
        ../../src  \
        ../../src/src_iterative_solver \
	../../src/src_fermisurfdx \
	../../src/src_xs \
        && make -f Makefile.parallel excitingmpi \
        && cp -f excitingmpi ..

libs:	libblas.a liblapack.a fftlib.a libbzint.a libarpack.a

libblas.a: make.inc
	cp -f make.inc ../..
	cd ../../src/BLAS && make clean && make
	cp ../../src/BLAS/libblas.a .
	cd ../../src/BLAS && make clean

liblapack.a: make.inc
	cp -f make.inc ../..
	cd ../../src/LAPACK && make clean && make
	cp ../../src/LAPACK/liblapack.a .
	cd ../../src/LAPACK && make clean

fftlib.a: make.inc
	cp -f make.inc ../..
	cd ../../src/fftlib && make clean && make
	cp ../../src/fftlib/fftlib.a .
	cd ../../src/fftlib && make clean

libbzint.a: make.inc
	cp -f make.inc ../..
	cd ../../src/src_libbzint && make clean && make
	cp ../../src/src_libbzint/libbzint.a .
	cd ../../src/src_libbzint && make clean

libarpack.a: make.inc
	cp -f make.inc ../..
	cd ../../src/ARPACK && make clean && make lib
	cp ../../src/ARPACK/libarpack.a .
	cd ../../src/ARPACK && make clean

clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a
