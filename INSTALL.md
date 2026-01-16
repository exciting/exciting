- [Compiling exciting](#compiling-exciting)
  - [Requirements](#requirements)
  - [Compiling (CMake)](#compiling-cmake)
    - [Compilation Steps](#compilation-steps)
      - [Notes](#notes)
    - [Example Configurations](#example-configurations)
      - [Intel Machines](#intel-machines)
      - [AMD-Based Machines](#amd-based-machines)
    - [CMake Options](#cmake-options)
    - [Mac OS](#mac-os)
  - [Documentation](#documentation)
  - [Using Apptainer containers](#using-apptainer-containers)
    - [Definition file](#definition-file)
      - [Intel oneAPI based containers](#intel-oneapi-based-containers)
      - [Containers based on GCC, OpenMPI, OpenBLAS](#containers-based-on-gcc-openmpi-and-openblas)
    - [Building the container image](#building-the-container-image)
    - [Compiling exciting after the image has been built](#compiling-exciting-after-the-image-has-been-built)
    - [Running the test suite](#running-the-test-suite)
  - [Compiling exciting with SIRIUS](#compiling-exciting-with-sirius)
  - [SIRIUS Gotchas](#sirius-gotchas)
  - [fastBSE](#fastbse)
  - [Known Issues](#known-issues)
    - [Intel MPI](#intel-mpi)
    - [FORD Fails to Find Graphviz](#ford-fails-to-find-graphviz)

# Compiling exciting

Requirements
------------------
exciting requires `xsltproc` to preprocess its XML schema into code.
Additionally, the code requires the installation of FFTW3 (such as oneMKL, AOCL-FFTW, FFTW3, Cray-FFTW, etc.) and a BLAS/LAPACK implementation (such as oneMKL, BLIS+libFLAME, OpenBLAS, LibSci, etc.) to compile. Make sure that you are using a multithreading-aware (thread-safe) BLAS/LAPACK library when enabling OpenMP.
**PLEASE** ensure you have these libraries and binaries installed before proceeding.

exciting comes with the following external libraries required to compile the code:

* [FoX XML](https://github.com/andreww/fox) library for parsing the input (2012 version).
		
* [LIBXC V7](https://libxc.gitlab.io/) library of DFT exchange and correlation functionals. It is also possible to link against an external Libxc installation, provided it is version 5.0.0 or newer.
	
* [BSPLINE-FORTRAN](https://github.com/jacobwilliams/bspline-fortran) Multidimensional B-Spline interpolation of data on a regular grid.

Compilation for fully parallel execution requires an MPI library, such as Open MPI, MPICH or Intel MPI library, and optionally a version of ScaLAPACK. These can be installed with package managers such as 
APT, Conda, Spack or EasyBuild, via a containerized environment, or built manually from source. 

Test suite dependencies are specified in [test/README](test/README).  


Compiling (CMake)
------------------

The `exciting` code is built with CMake. The following compilers are supported:
- Intel Classic (ifort): 2021.0.3, 2021.13.1
- Intel LLVM (ifx): 2025.0.0 
- GNU: 11, 12, 14, 15
- Cray: 18.0.1, 19.0.1
- LLVM-Flang-based compilers (must support Fortran 2018 standard)

### Compilation Steps

To compile `exciting`, run the following commands from the `exciting` root directory (**Note that the following is not a working example, for these refer to the next section**):
```shell
  mkdir build
  cd build
  FC=[SERIAL_FORTRAN_COMPILER] CC=[SERIAL_C_COMPILER] CXX=[SERIAL_C++_COMPILER] ../external/cmake-3.31.3-linux-x86_64/bin/cmake [OPTIONS] ..
  make -j N -l N exciting_NAME
  make install
```
#### Notes:
- `exciting_NAME` is determined by CMake during configuration, depending on the selected options:
  - `exciting_serial` (`-DOMP=OFF -DMPI=OFF`)
  - `exciting_smp` (`-DOMP=ON -DMPI=OFF`)
  - `exciting_mpismp` (`-DOMP=ON -DMPI=ON`)
- We provide a bundled version of CMake located at `external/cmake-3.31.3-linux-x86_64/bin/cmake`.
  **Note:** This version is only valid for `x86_64` systems.
- The `[SERIAL_FORTRAN_COMPILER]`, `[SERIAL_C_COMPILER]`, and `[SERIAL_C++_COMPILER]` must be adjusted to the compilers to use.
- The `[OPTIONS]` section must be adjusted based on your processor, compilers, and/or required features.

### Example Configurations

#### Intel Machines:
- **Classic Intel Compilers:**   
**N.B.:** The classic Intel compilers are deprecated. For instance, `icc` and `icpc`
were removed from all Intel `oneAPI` packages starting with the `2024.0` release, 
and `ifort` was removed starting with the `2025.0` release.   
Unless BSE calculations are planned, it is recommended to use the Intel 
LLVM-based compilers instead.
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpiifort -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DMKL=ON ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **Intel LLVM Compilers with a processor supported by `-ax` (e.g., Intel(R) Xeon(R) Platinum 8480L, code name SAPPHIRERAPIDS):**   
**N.B.:** Please keep in mind that BSE results are not fully reproducible when 
using Intel LLVM compilers.
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpiifx -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DMKL=ON -DINTEL_CODE_NAME=SAPPHIRERAPIDS ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **Intel LLVM Compilers with an unsupported processor name:**
```shell
mkdir build
cd build
FC=ifx CC=icx CXX=icpx ../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpiifx -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DMKL=ON ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
#### AMD-Based Machines:
- **With OpenBLAS and ScaLAPACK:**
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DOPENBLAS=ON -DSCALAPACK=ON ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **With AOCL-FFTW and BLIS + libFLAME:**
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DAMDLINALG=ON ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **With OpenBLAS+SIRIUS+HDF5**
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DOPENBLAS=ON -DSIRIUS=ON -DHDF5=ON ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **GPU accelerated: AMD GPU MI250**
AMD GPU acceleration requires the MAGMA library, assuming it is installed in _MAGMA\_INSTALL\_DIR_
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DOPENBLAS=ON -DAMD=ON -DAMDTARGET=gfx90a -DMAGMA_ROOT=MAGMA_INSTALL_DIR ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **GPU accelerated with USM -  MI300A**
Note that MI300A is an APU and thus the RAM is physically shared between the GPU and CPU. _USM_ option is mandatory in those cases.
AMD GPU acceleration requires of MAGMA library. Assuming it is installed in _MAGMA_INSTALL_DIR_
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DOPENBLAS=ON -DUSM=ON -DAMD=ON -DAMDTARGET=gfx942 -DMAGMA_ROOT=MAGMA_INSTALL_DIR ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```
- **GPU accelerated: Intel PVC with Intel Sapphire Rapids**
```shell
mkdir build
cd build
../external/cmake-3.31.3-linux-x86_64/bin/cmake -DCMAKE_Fortran_COMPILER=mpiifx -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DMKL=ON -DINTEL=ON -DINTEL_CODE_NAME=SAPPHIRERAPIDS ..
make -j`nproc` -l`nproc` exciting_mpismp
make install
```

### CMake Options

CMake installation can be customized using the following options (**Notice that in CMake options are written as -DOPTION=OPTION_VALUE**):

- **_MPI_**: Controls MPI support (default: ON).
- **_OMP_**: Controls OpenMP support (default: ON).
- **_HDF5_**: Enables HDF5 support (default: OFF).
- **_IDIEL_**: Enables IDieL library for the anisotropic averaging of the RPA dielectric matrix (default: ON).
- **_MKL_**: Uses MKL for linear algebra and FFT (default: OFF).
- **_OPENBLAS_**: Uses OpenBLAS for linear algebra (default: OFF).
- **_AMDLINALG_**: Uses AMD linear algebra libraries (BLIS and FLAME) (default: OFF).
- **_CRAYLIBSCI_**: Uses Cray LibSci for linear algebra and, if required, ScaLAPACK (default: OFF).
- **_OTHERLINALG_**: Uses another linear algebra library (not officially supported) (default: OFF).
- **_LINALGLIB_**: If `OTHERLINALG` is ON, specify the full path to the desired linear algebra libraries (default: None).
- **_FFTW3_ROOT_**: For non-standard compilation, provide the install directory of FFTW3 (default: None).
- **_SCALAPACK_**: Enables ScaLAPACK support (default: OFF).
- **_SCALAPACK_ROOT_**: For non-standard compilation, provide the install directory of ScaLAPACK (default: None).
- **_USE_INTERNAL_LIBXC_**: Uses the bundled libXC version (default: ON).
- **_LIBXC_ROOT_**: For non-standard compilation, provide the install directory of libXC (default: None).
- **_SIRIUS_**: Compiles EXCITING with SIRIUS (default: OFF).
- **_NVIDIA_**: Enables GPU support for NVIDIA GPUs (default: OFF).
- **_NVIDIAARCH_**: Sets the NVIDIA architecture (default: 89).
- **_AMD_**: Enables GPU support for AMD GPUs (default: OFF).
- **_AMDTARGET_**: Specifies the AMD GPU target (e.g., gfx90a) (default: None).
- **_AMD_HIPSETVALIDDEVICE_SUPPORTED_**: Set to ON if `hipSetValidDevices` is supported (after ROCm 6.2.0) (default: OFF).
- **_INTEL_**: Enables GPU support for Intel GPUs (default: OFF).
- **_USM_**: Enables an optimized compilation for CPU-GPU systems with USM; i.e. physically shared RAM (default: OFF).
- **_CPUBACKEND_**: Enables a CPU-only build (default: ON).
- **_MAGMA_ROOT_**: If `AMD` or `NVIDIA` are ON, and MAGMA is not in a default path this provides the path to its install directory (default: None).
- **_INTEL_CODE_NAME_**: For Intel processors, this can be modified to match the processor name, allowing `ifx` to generate optimized code paths. If not set, defaults to the host machine Intel architecture. **Do not modify for non-Intel machines.** (default: Host).
- **_DOCUMENTATION_**: Controls whether documentation is built (default: None).
- **_UNIT_TESTS_**: Enables unit tests via make test command (default: ON).
- **_REGRESSION_TESTS_**: Enables regression tests via make test. This requires Python 3 (default: OFF).
- **_BUILD_EXCITING_**: Builds EXCITING (default: ON).
- **_ELPA_**: Enables ELPA support (default: OFF).
- **_ELPA_ROOT_**: Root path to the ELPA installation when `ELPA` is ON (default: None).

On top of the specific options for _exciting_, the following CMake default options can be modified to tune your installation:
- **_CMAKE_INSTALL_PREFIX_**: controls the installation directory (default: **install** in the exciting root directory).
- **_CMAKE_BUILD_TYPE_**: controls the installation type, it can be: Debug or Release (default: Release).
- **_CMAKE_Fortran_FLAGS_**: string to add extra compiler options for the Fortran compiler. Note that this affects the whole compilation, and that the default flags are generating high performant executable for the _Release_ build.
### Mac OS

exciting can be compiled on mac OS, but it is complicated by nonstandard installation
locations for libraries, which can vary between OS versions and package managers.

Because production calculations typically require HPC resources, compilation on mac is
primarily intended for developers.

For mac, exciting has been tested on:

* OS Monterey using GCC 12.2.0, openBLAS 0.3.21 and openmpi 4.1.4

Documentation
------------------

exciting's documentation can be built with FORD:

```shell
  cd build
  make GenerateExcitingDocs
```
Notice that this requires the option `-DDOCUMENTATION=ON` in the configuration step.

FORD is available as a python package and can be installed with pip. To install FORD, type:

```shell
  pip3 install ford
```

FORD generates html-based documentation, including graphical dependency analysis, which can be 
viewed by opening docs/exciting_ford/index.html in a web browser. More details of FORD can be found 
on its [GitHub page](https://github.com/Fortran-FOSS-Programmers/ford), and additional details regarding
installation of dependencies can be found under 'Known Issues', below. 


## Using Apptainer containers
This section illustrates how to install `exciting` using an [Apptainer](https://apptainer.org/) container.

**Prerequisite:** 
Ensure that `Apptainer` is installed on your system.

### Definition file

#### Intel oneAPI-Based Containers

For Intel compilers, use the `oneapi-hpckit` Docker image.  
Example definition file (`my-container.def`) for the Intel LLVM compiler (`ifx`):

```shell
Bootstrap: docker
From: intel/oneapi-hpckit:2025.1.0-0-devel-ubuntu24.04

%post
    # Update package manager and install xsltproc
    apt-get update && apt-get install -y xsltproc python3-dev python3-pytest python3-numpy
```

For the classic Intel compiler (`ifort`):
```shell
Bootstrap: docker
From: intel/oneapi-hpckit:2024.0.1-devel-ubuntu22.04

%post
    wget -qO - https://repositories.intel.com/gpu/intel-graphics.key | \
      gpg --yes --dearmor --output /usr/share/keyrings/intel-graphics-archive-keyring.gpg && \
    apt-get update && apt-get install -y xsltproc python3-dev python3-pytest python3-numpy cmake
```

#### Containers Based on GCC, OpenMPI, and OpenBLAS
To build a container using GCC, OpenMPI, and OpenBLAS, use the following definition file (`my-container.def`) based on the `xinzhewu/hpc-gcc12-openmpi5-openblas-scalapack:latest` Docker image:

```shell
Bootstrap: docker
From: xinzhewu/hpc-gcc12-openmpi5-openblas-scalapack:latest

%post
    # Update package manager and install xsltproc
    apt-get update && apt-get install -y xsltproc python3-dev python3-pytest python3-numpy
    # Install FFTW
    wget http://www.fftw.org/fftw-3.3.10.tar.gz
    tar -xzf fftw-3.3.10.tar.gz && cd fftw-3.3.10 && mkdir build && cd build 
    cmake .. -DENABLE_THREADS=ON -DENABLE_OPENMP=ON -DENABLE_AVX2=ON && make && make install
    cd ../ && mkdir build_float && cd build_float
    cmake .. -DENABLE_THREADS=ON -DENABLE_OPENMP=ON -DENABLE_AVX2=ON -DENABLE_FLOAT=ON && make && make install
    cd ../.. && rm -rf fftw-3.3.10 fftw-3.3.10.tar.gz
```
**N.B.:** If the target system does not support `AVX2`, remove `-DENABLE_AVX2=ON` from CMake commands above.

### Building the container image
To build the container, run:
```bash
apptainer build my-container.sif my-container.def
```
Sometimes, the default temporary working space (typically `/tmp`) may not have
enough space to accommodate the full uncompressed container image with other 
temporary files generated during the build process. This can cause the build to fail 
with an error indicating insufficient disk space. 

To work around this, you can specify a custom temporary directory with sufficient 
space:
```bash
mkdir -p $HOME/temp/
apptainer build --tmpdir $HOME/temp/ my-container.sif my-container.def
```
After the build is complete, you may remove the temporary directory if it is no 
longer needed: 
```bash
rm -rf $HOME/temp/
```

### Compiling exciting after the image has been built
To compile `exciting`, follow the steps below. These instructions assume: 
  1. The container file `my-container.sif` is located in your `home` directory.
  2. The `exciting` source code is located in `$HOME/exciting`, and
  3. The build and install directories are `$HOME/build-farm/build-$branch` and 
  `$HOME/build-farm/install-$branch`, respectively, where `$branch` can be used 
  to represent the `git` branch being compiled.
If your setup differs, please adjust the paths accordingly. In addition, use the following table to set environment variables for the chosen compiler and container:

| Variable | `ifx` | `ifort` | `gfortran`|
|--|--|--|--|
| `export FC=`| `ifx`| `ifort`| `gfortran`|
| `export MKL=`| `ON`| `ON`| `OFF`|
| `export fftw3=`| `""`| `""`| `"-DFFTW3_ROOT=/usr/local/"`|
| `export env=`| `""`| `""`| `"--env LD_LIBRARY_PATH=$install_dir/lib/"`|


```bash
export container="$HOME/my-container.sif"
export exciting_dir="$HOME/exciting"
export branch="mybranch"
export build_dir="$HOME/build-farm/build-$branch"
export install_dir="$HOME/build-farm/install-$branch"
export cmake_custom=$exciting_dir/external/cmake-3.31.3-linux-x86_64/bin/cmake 
export FC="(adjust accordingly)"
export MKL="(adjust accordingly)"
export fftw3="(adjust accordingly)"
cd $exciting_dir
apptainer exec $container $cmake_custom \
  -DMKL=$MKL -DMPI=ON -DOMP=ON \
  -S $exciting_dir \
  -B $build_dir \
  -DCMAKE_INSTALL_PREFIX=$install_dir \
  $fftw3 --trace --fresh
apptainer exec $container $cmake_custom --build $build_dir -j 4
apptainer exec $container $cmake_custom --install $build_dir
```

Finally, if the installation completes successfully, you can execute `exciting` using:
```bash
export env="(adjust accordingly)"
apptainer exec $env $container mpirun -np 4 $install_dir/bin/exciting_mpismp
```
where the `-np` argument specifies the number of MPI processes to launch.   
**N.B.:** This approach launches `mpirun` from within the container. This is acceptable for single-node runs, but not suitable for multi-node executions.

### Running the test suite
```bash
export env="(adjust accordingly)"
python3 -m venv $HOME/.excitingvenv
source $HOME/.excitingvenv/bin/activate
pip install pyyaml
pip install -e $exciting_dir/tools/exciting_tools/
cd $exciting_dir/test
python runtest.py -np 2 -omp 2 -mpirun "apptainer exec $env $container mpirun" -bp $install_dir/bin -e exciting_mpismp
```

Compiling exciting with SIRIUS
-----------------------------------

[SIRIUS](https://github.com/electronic-structure/SIRIUS) is a domain specific library for electronic 
structure calculations. It implements pseudopotential plane wave (PP-PW) and full potential linearized 
augmented plane wave (FP-LAPW) methods, and is designed for GPU acceleration of popular community codes 
such as Exciting, Elk and Quantum ESPRESSO.

Compiling exciting with SIRIUS is complex. To simplify the process of building dependencies, SIRIUS can be
completely installed with the python package manager [spack](https://spack.readthedocs.io/en/latest/getting_started.html).

As of exciting Sodium, only a CPU build chain using GCC on Ubuntu Focal is regularly tested in exciting's 
CI. This is provided in [build/utilities/docker/Dockerfile_ci_sirius](build/utilities/docker/Dockerfile_ci_sirius). 

```shell
# Install sirius dependencies with spack
# ---------------------------------------
# Install spack: https://spack.readthedocs.io/en/latest/index.html

# Find any preinstalled compilers or dependencies
# Note, it's preferable to use a preinstalled compiler where possible, to minimise total build time
spack compiler find
spack external find

# Define the spack specification for sirius as a string
# See: https://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies
# Remove `target` or adjust accordingly, if not on x86_64 hardware
export SPEC="sirius@develop build_type=Release +scalapack +fortran ^mpich@3.3.2 ^intel-oneapi-mkl+cluster ^spfft target=x86_64"

# Install sirius's dependencies
spack install $SPEC

# Load SIRIUS
spack load sirius

# Compile it with SIRIUS options (-DSIRIUS=ON -DUSE_INTERNAL_LIBXC=OFF). Notice that we need to compile exciting with the same libXC version than SIRIUS. 
```

SIRIUS Gotchas
------------------

* `spack` command is not added to the `$PATH` by default. To source `spack` in bash, type 
  `source spack/share/spack/setup-env.sh`

* sirius is under active development, and its dependencies change regularly. Please use `spack info sirius`
  for the current status

* spack does not query a server when `spack info` is called. To ensure the info is up-to-date, run
  `git pull` in the root of the spack directory.

fastBSE
------------------
To use fastBSE, you need to build exciting with HDF5 and FFTW3. While FFTW3 is 
used by default, HDF5 support must be activated during the cmake step by 
setting the build option `-DHDF5=ON`. Do not forget to install the HDF5 
library or to load the appropriate module.


**GCC**

To install the libraries in a Debian-based distribution
```shell
  sudo apt-get install -y libhdf5-serial-dev libhdf5-mpi-dev libfftw3-dev
```

Then, to compile `exciting`
```bash
  mkdir build
  cd build
  FC=gfortran CC=gcc CXX=g++ ../external/cmake-3.31.3-linux-x86_64/bin/cmake -DHDF5=ON ..
  make -j N -l N exciting_mpismp
  make install
```

**INTEL**

If you intend to use Intel, you must manually compile HDF5. We recommend using HDF5 1.12.0. The source code can be found on the [HDF5 website](https://www.hdfgroup.org/solutions/hdf5/). 
We recommend configuring the build as follows:

```shell
  configure --enable-build-mode=production --enable-parallel --enable-fortran --enable-fortran2003
```

Alternatively, one can use spack with the following options `cxx=true fortran=true mpi=true %oneapi@2025.0.0`. Notice that you need to replace the compiler
version by the correct tag. See $HOME/.spack/linux/compilers.yaml for the compiler tags. In case, your Intel compilers do not appear, make sure they are loaded
and execute `spack compiler find`.

For details on the options, refer to the HDF5 documentation. 
After a successful installation, update the PATH variable with the location of the HDF5 binaries:

```shell
  export PATH="path/to/hdf5-build/hdf5/bin:$PATH"
```

Please note that this code is not executable as-is; it provides instructions for setting up the fastBSE environment. 
Make sure to adjust the paths and options according to your system and requirements.

Known Issues
------------------

### Intel MPI

Intel MPI 2019 and 2021.3 are known to exhibit memory leaks. This can be troublesome for memory-intensive calculations
and cause exciting to crash. See the threads associated with [Qbox](https://groups.google.com/g/cp2k/c/BJ9c21ey0Ls) and
[CP2K](https://github.com/cp2k/cp2k/issues/1830) for further discussion.

### FORD Fails to Find Graphviz

On Debian Buster, we have noticed that FORD fails to find graphviz, despite pip installing
it as part of FORD's dependencies. This can be fixed by installing graphviz with APT. Optionally
one can also install python-dev and LXML, which facilitate faster generation of the search database:

```shell
  sudo apt install python3-all-dev graphviz
  pip3 install lxml
```

