F90 =mpif90
F90_OPTS = -g -DMPI -DMPISEC -DMPIRHO
F77 =ifort 
F77_OPTS = 
LIB_SYS = 
LIB_LPK = -llapack -lblas  libbzint.a libarpack.a
LIB_FFT = fftlib.a
