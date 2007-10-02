
include make.inc

serial:
	mkmf -t template -f -m Makefile.serial -p excitings \
        ../.. ../../src_tddft && make -f Makefile.serial excitings \
        && cp -f excitings ..

libs:	blas lapack fftlib libbzint

blas:
	cp -f make.inc ../../..
	cd ../../BLAS && make clean && make
	cp ../../BLAS/blas.a .
	cd ../../BLAS && make clean

lapack:
	cp -f make.inc ../../..
	cd ../../LAPACK && make clean && make
	cp ../../LAPACK/lapack.a .
	cd ../../LAPACK && make clean

fftlib:
	cp -f make.inc ../../..
	cd ../../fftlib && make clean && make
	cp ../../fftlib/fftlib.a .
	cd ../../fftlib && make clean

libbzint:
	cp -f make.inc ../../..
	cd ../../src_libbzint && make clean && make
	cp ../../src_libbzint/libbzint.a .
	cd ../../src_libbzint && make clean

clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a
