#!/usr/bin/env bash

#------------------------------------------------------------------------
# Build Stack using Intel 2019, on Debian Buster. 
# 
# Run with `source sol-intel-2019.sh`
# else the module command may not correctly execute.
#
# After installing, the following:
#  ```
#  spack env activate exciting_intel-2019
#  ```
# can be added to the ~/.bashrc such that the virtual 
# environment loads with the shell 
#
# Note: 
# m4 has a compilation bug with Intel 2019. 
# See, for example:  https://github.com/spack/spack/issues/17877
# Hence this script modifies m4's package.py before starting a virtual env
#
# Installs:
#  hdf5@1.12.0
#------------------------------------------------------------------------

# Remove the CFLAGS to enable Intel compilation of m4
# \x27 =  single quote in ASCII
M4_PACKAGE="${SPACK_ROOT}/var/spack/repos/builtin/packages/m4"
cd $M4_PACKAGE
sed -i.bak 's/args.append(\x27CFLAGS=-no-gcc\x27)/pass/g' package.py 
cd - 

spack env create exciting_intel-2019
spack env activate exciting_intel-2019
module load intel/2019
spack compiler find    
spack install hdf5@1.12.0+fortran+hl+mpi ^intel-mpi@2019.0.117 %intel@19.0.0.117