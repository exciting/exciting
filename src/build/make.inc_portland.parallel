
# ---------------------------- PGI compiler -----------------------------------
F90		= pgf95
F90_OPTS	= -DMPI -fastsse -O4
F77 		= pgf95
F77_OPTS 	= -fastsse -O4
LIB_SYS 	= /installadmin/mpich/1.2.6/pgi/lib/libmpich.a
### the path below might not be consistent with the 6.2 version of PGI
##LIB_SYS 	= /software/pgi/linux86-64/6.0/lib/libmpich.a

### BLAS and LAPACK
# included libraries
LIB_LPK = lapack.a blas.a

# system libraries
#LIB_LPK = -llapack -lblas

### FFT library
LIB_FFT = fftlib.a

### Library for Brillouin zone integration
LIB_BZINT = libbzint.a
# ---------------------------- PGI compiler -----------------------------------
