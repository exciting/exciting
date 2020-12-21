# Example make.inc for Intel 2019

# Compiler, flags and preprocessing variables 
F90 = ifort
F90_OPTS = -O3 -ip -unroll -scalar_rep -I"${MKLROOT}/include"
CPP_ON_OPTS = -cpp -DXS -DISO -DLIBXC -DIFORT -DFFTW
F77 = $(F90)
F77_OPTS = -O3
# Required for libXC compilation 
FCCPP = cpp

F90_DEBUGOPTS = -g -O0 -DUSE_ASSERT -debug all -implicitnone -warn unused \
   -fp-stack-check -heap-arrays -ftrapuv -check pointers \
   -check bounds -check all -check noarg_temp_created -traceback -I"${MKLROOT}/include"

F77_DEBUGOPTS = $(F90_DEBUGOPTS)

#Libraries
LIB_ARP =libarpack.a
export USE_SYS_LAPACK=true
LIB_LPK = -mkl=parallel
LIB_FFT = fftlib.a
LIB_BZINT = libbzint.a
LIBS = $(LIB_ARP) $(LIB_LPK) $(LIB_FFT) $(LIB_BZINT)

# SMP and MPI compilers, flags and preprocessing variables 
MPIF90 = mpiifort
MPIF90_OPTS = -DMPI -DSCAL
MPI_LIBS = -L${MKLROOT}/lib/intel64 -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64


SMPF90_OPTS = -qopenmp -DUSEOMP
SMPF77_OPTS = $(SMPF90_OPTS)
SMP_LIBS = 

BUILDMPI = true
BUILDSMP = true
BUILDMPISMP = true
