# ---------------------------- PGI compiler -----------------------------------
F90 = pgf95
F90_OPTS = -DMPI -fastsse -O4
F77 = pgf95
F77_OPTS = -fastsse -O4
LIB_SYS = /installadmin/mpich/1.2.6/pgi/lib
### the path below might not be consistent with the 6.2 version of PGI
LIB_SYS = /software/pgi/linux86-64/6.0/lib/libmpich.a
# use custom blas and lapack
LIB_LPK = lapack.a blas.a
#LIB_LPK = -llapack -lblas
#LIB_LPK = -L /software/intel/mkl721/lib/em64t/ -lmkl_lapack -lmkl_em64t -lguide -lpthread
LIB_FFT = fftlib.a libbzint.a
# ---------------------------- PGI compiler -----------------------------------
