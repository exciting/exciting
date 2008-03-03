include make.inc

MAKE=make

serial:make.inc
	../../utilities/mkmf -t template -f -m Makefile.serial -p excitings \
	../../src  \
	../../src/src_xs \
	../../src/src_iterative_solver \
	../../src/src_iterative_solver/src_jacobidavidson \
	../../src/src_iterative_solver/src_diis \
	&& make -f Makefile.serial excitings \
	&& cp -f excitings ..

libs:	libblas.a liblapack.a fftlib.a libbzint.a libarpack.a libjdqz.a

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
	cd ../../src/ARPACK && $(MAKE) clean && make lib
	cp ../../src/ARPACK/libarpack.a .
	cd ../../src/ARPACK && $(MAKE) clean
	
libjdqz.a: make.inc
	cp -f make.inc ../..
	cd ../../src/JDQZ && $(MAKE) clean && make 
	cp ../../src/JDQZ/jdlib/libjdqz.a ./
	cd ../../src/JDQZ && $(MAKE)  clean


clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a

### build libbzint, do not clean source tree of library and build exciting
bzex:
	cp -f make.inc ../..
	cd ../../src/src_libbzint && make
	cp ../../src/src_libbzint/libbzint.a .
	rm -f exciting
	make serial
