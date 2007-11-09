F90 =ifort
F90_OPTS = -g  #-DDEBUG
F77 =ifort 
F77_OPTS =  -g -O3
LIB_SYS = 
LIB_LPK = liblapack.a libblas.a libarpack.a -llapack -lblas
LIB_FFT = fftlib.a
LIB_BZINT=libbzint.a
