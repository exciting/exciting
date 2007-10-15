include make.inc


parallel:
	../../utilities/mkmf -t template -f -m Makefile.parallel -p excitingp \
        ../../src  \
	../../src/src_xs \
        && make -f Makefile.parallel excitingp \
        && cp -f excitingp ..

libs:	blas lapack fftlib 

blas:
	cp -f make.inc ../..
	cd ../../src/BLAS && make clean && make
	cp ../../src/BLAS/blas.a .
	cd ../../src/BLAS && make clean

lapack:
	cp -f make.inc ../..
	cd ../../src/LAPACK && make clean && make
	cp ../../src/LAPACK/lapack.a .
	cd ../../src/LAPACK && make clean

fftlib:
	cp -f make.inc ../../..
	cd ../../src/fftlib && make clean && make
	cp ../../src/fftlib/fftlib.a .
	cd ../../src/fftlib && make clean

libbzint:
	cp -f make.inc ../..
	cd ../../src/src_libbzint && make clean && make
	cp ../../src/src_libbzint/libbzint.a .
	cd ../../src/src_libbzint && make clean

clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a

