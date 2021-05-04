#!/usr/bin/env bash

#--------------------------------------------------------
# Build Stack using GCC 8.3.0 (Debian Buster default GCC)
#
# Run with `source linux-gcc-8-3-0.sh` 
# 
# After installing, the following:
#  ```
#  spack env activate exciting_gcc-8-3-0
#  ```
# can be added to the ~/.bashrc such that the virtual 
# environment loads with the shell 
#
# Installs:
#  openblas@0.3.10 (one could switch this for MKL)
#  openmpi@4.0.4
#  scalapack@2.1.0
#  hdf5@1.12.0
#--------------------------------------------------------
spack env create exciting_gcc-8-3-0
spack env activate exciting_gcc-8-3-0
spack install hdf5@1.12.0+fortran+hl+mpi ^openmpi@4.0.4 %gcc@8.3.0
spack install openblas@0.3.10 %gcc@8.3.0
spack install netlib-scalapack@2.1.0 ^openmpi@4.0.4 %gcc@8.3.0