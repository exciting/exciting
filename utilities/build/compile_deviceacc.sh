#!/bin/bash
#
# This bash script compiles deviceacc
# library using the make.inc build system
# WARNING: Only CPU support, for GPU support use CMake to build exciting


# Checking if executed with the correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Error: Expected 1 argument, but received $#"
    exit 1
fi

# Getting the compiler from the variable
fortran_compiler=$1

# Getting the folders
compilation_folder=$PWD
cd ../../
exciting_root=`echo $PWD`
deviceacc_folder="$exciting_root/external/deviceacc"

# Change to the deviceacc folder
cd $deviceacc_folder
# Remove old install
if [ -d "build_deviceacc" ]; then
    rm -rf "build_deviceacc"
fi
if [ -d "deviceacc" ]; then
    rm -rf "deviceacc"
fi

mkdir build_deviceacc
cd build_deviceacc
FC=$fortran_compiler cmake -DCPUBACKEND=ON ..
make -j4 -l4
make install

# Return to the original folder
cd $compilation_folder

# Providing access to deviceacc modules and library
ln -s $deviceacc_folder/deviceacc/include/* .
ln -s $deviceacc_folder/deviceacc/lib/* .

