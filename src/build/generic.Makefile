
# Copyright (C) 2007 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
# C. Ambrosch-Draxl.
# This file is distributed under the terms of the GNU General Public License.
# See the file COPYING for license details.

# Generic Makefile
# To be used from within the serial/ and parallel/ directory.
# It utilized the "mkmf" perl script to generate Makefile
# taking into account dependencies due to Fortran modules.

# default name of executable
EXE_name = exciting

# include compiler/linker/archiver flags
INCLUDEFILE = make.inc
TEMPLATEFILE = ../mkmf.template
include $(INCLUDEFILE)

# directories
DIR_root	= ../../../
DIR_main 	= ../../
DIR_blas 	= ../../BLAS/
DIR_lapack 	= ../../LAPACK/
DIR_fftlib 	= ../../fftlib/
DIR_libbzint	= ../../src_libbzint/
DIR_xs 		= ../../src_xs/

# temporary Makefile name
MKMF_MAKEFILE = Makefile.mkmf

# library file names
LIB_blas 	= blas.a
LIB_lapack 	= lapack.a
LIB_fftlib 	= fftlib.a
LIB_libbzint 	= libbzint.a

exciting:
	mkmf -t $(TEMPLATEFILE) -f -m $(MKMF_MAKEFILE) -p $(EXE_name) \
        $(DIR_main) $(DIR_xs) && make -f $(MKMF_MAKEFILE) $(EXE_name) 

#\
#&& mv -f $(EXE_name) ..

libs:	blas lapack fftlib libbzint

blas:
	cp -f $(INCLUDEFILE) $(DIR_root)
	cd $(DIR_blas) && make clean && make
	cp $(DIR_blas)/$(LIB_blas) .
	cd $(DIR_blas) && make clean

lapack:
	cp -f $(INCLUDEFILE) $(DIR_root)
	cd $(DIR_lapack) && make clean && make
	cp $(DIR_lapack)/$(LIB_lapack) .
	cd $(DIR_lapack) && make clean

fftlib:
	cp -f $(INCLUDEFILE) $(DIR_root)
	cd $(DIR_fftlib) && make clean && make
	cp $(DIR_fftlib)/$(LIB_fftlib) .
	cd $(DIR_fftlib) && make clean

libbzint:
	cp -f $(INCLUDEFILE) $(DIR_root)
	cd $(DIR_libbzint) && make clean && make
	cp $(DIR_libbzint)/$(LIB_libbzint) .
	cd $(DIR_libbzint) && make clean && rm -f $(LIB_libbzint)

clean:
	rm -f *.o *.mod

cleanlibs:
	rm -f *.a
