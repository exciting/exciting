include make.inc


parallel:
	../../utilities/mkmf -t template -f -m Makefile.parallel -p excitingmpi\
        ../../src  \
        && make -f Makefile.parallel excitingmpi \
        && cp -f excitingmpi ..

libs:	blas lapack fftlib 

blas:
	cp -f make.inc ../..
	cd ../../src/BLAS && make clean && make
	cp ../../src/BLAS/libblas.a .
	cd ../../src/BLAS && make clean

lapack:
	cp -f make.inc ../..
	cd ../../src/LAPACK && make clean && make
	cp ../../src/LAPACK/liblapack.a .
	cd ../../src/LAPACK && make clean

fftlib:
	cp -f make.inc ../..
	cd ../../src/fftlib && make clean && make
	cp ../../src/fftlib/fftlib.a .
	cd ../../src/fftlib && make clean



clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a
