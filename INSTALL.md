# Compiling exciting


Requirements
------------------
exciting requires `xsltproc` to preprocess its XML schema into code.
Additionally, the code requires the installation of FFTW3 (such as oneMKL, AOCL-FFTW, FFTW3, Cray-FFTW, etc.) and a BLAS/LAPACK implementation (such as oneMKL, BLIS+libFLAME, OpenBLAS, LibSci, etc.) to compile. Be aware, of using an multithreading-aware version of BLAS/LAPACK libraries.
**PLEASE** ensure you have these libraries and binaries installed before proceeding.

exciting comes with the following external libraries required to compile the code:

* [FoX XML](https://github.com/andreww/fox) library for parsing the input (2012 version).
		
* [LIBXC V7](https://libxc.gitlab.io/) library of DFT exchange and correlation functionals. We also allow to use external libXC, given that they are version 5.0.0 or higher.
	
* [BSPLINE-FORTRAN](https://github.com/jacobwilliams/bspline-fortran) Multidimensional B-Spline interpolation of data on a regular grid.

Compilation for fully parallel execution requires an MPI library, such as open MPI, MPICH or Intel MPI library and optionally, a version of Scalapack. These can be installed with a package manager such as 
APT, Conda, Spack or EasyBuild, or built manually from source. Some limited spack recipes for 
installing external libraries, compiled with GCC and Intel, are provided in the [repository](build/utilities/spack). 

Test suite dependencies are specified in [test/README](test/README).  

exciting can be built using CMake. 

Compiling (CMake)
------------------

The code can be compiled using CMake. Be aware that only Intel (classic and LLVM), GNU, Cray and LLVM-Flang-based compilers are supported. Note that LLVM-Flang-based compilers need to support the Fortran2018 standard.

To compile exciting using CMake run in exciting root directory:
```shell
  mkdir build
  cd build
  ../external/cmake-3.31.3-linux-x86_64/bin/cmake [OPTIONS] ..
  make -j N -l N exciting_NAME
  make install
```
Where `exciting_NAME` is determined by CMake on configuration step depending on the options: `exciting_serial` (-DOMP=OFF -DMPI=OFF), `exciting_smp` (-DOMP=ON -DMPI=OFF), or `exciting_mpismp` (-DOMP=ON -DMPI=ON). Notice that we provide a bundled version of cmake with exciting `external/external/cmake-3.31.3-linux-x86_64/bin/cmake` for the build. Be aware that this version is only valid for `x86_64` systems.

CMake installation can be tuned with the following options:
* _MPI_: controls the MPI support (default: ON).
* _OMP_: controls the OpenMP support (default: ON)
* _HDF5_: enables the HDF5 support (default: OFF)
* _MKL_: use MKL for linear algebra and FFT (default: OFF)
* _OPENBLAS_: use OpenBLAS for linear algebra, as oposed to normal call this will look for accelerated versions (default: OFF)
* _AMDLINALG_: use AMD linear algebra libraries, i.e. BLIS and FLAME (default: OFF).
* _OTHERLINALG_: use another linear algebra library, as AMD (BLIS + libFLAME) or Cray libsci. We do not officially support exciting reliability with these libs; indeed some libsci versions are known to produce garbage results in exciting for multithreading (default: OFF)
* _LINALGLIB_: If _OTHERLINALG_ is ON, this option should contain the full path to the desired linear algebra libraries (default: None)
* _FFTW3\_ROOT_ : For non-standard compilation provide a PATH to the install directory of FFTW3 (default: None). 
* _SCALAPACK_ : Build exciting with scalapack support (default: OFF).
* _SCALAPACK\_ROOT_ : For non-standard compilation provide a PATH to the install directory of scalapack (default: None).
* _CRAY\_SCALAPACK_: Indicates we are using the scalapack functions from Cray LibSci library (default: OFF).
* _USE\_INTERNAL\_LIBXC_ : Use the libXC version bundled with exciting (default: ON).
* _LIBXC\_ROOT_ : For non-standard compilation provide a PATH to the install directory of libxc (default: None).
* _SIRIUS_ : compile exciting with SIRIUS (default: OFF). 
* _NVIDIA_: activate GPU support for NVIDIA GPUs (default: OFF)
* _NVIDIAARCH_: set the proper NVIDIA architecture (default: 89)
* _AMD_: activate GPU support for AMD GPUs (default: OFF)
* _AMDTARGET_: provides info about the AMD GPU; e.g. gfx90a (default: None)
* _AMD\_HIPSETVALIDDEVICE\_SUPPORTED_: set ON if hipSetValidDevices is supported; only after ROCm 6.2.0 (default: OFF)
* _INTEL_: activate GPU support for INTEL GPUs (default: OFF)
* _CPUBACKEND_: CPU-only build (default: ON)
* _MAGMA\_DIR_: If _AMD_ or _NVIDIA_ are ON, this option provides the path to MAGMA's install directory (default: None)
* _INTEL\_CODE\_NAME_: For Intel processors this can be modified to mach the processor name, so ifx can generate code paths specific for that processor. If not set the build suit will select generic Intel subset based on the presence of AVX512 and/or AVX2 instructions. Do not modify for non-Intel machines. (default: None)
* _DOCUMENTATION_: Controls if the documentation is built or not. (default: None).
* _UNIT\_TESTS_: Activate the unit tests for the build system. Requires python3 in the system (default: ON).
* _REGRESSION\_TESTS_: Activate the regression tests for the build system  (default: ON).
* _BUILD\_EXCITING_: Build exciting (default: ON).

### CMake build examples

Notice that the following examples assume the required modules/libraries to be loaded, and that the used `cmake` is 
the one bundled with `exciting`, i.e. `external/cmake-3.31.3-linux-x86_64/bin/cmake`. In all cases, the command is run after 
executing `mkdir -p build; cd build`.

* Intel machine (Classic Intel Compilers): `FC=ifort CC=icc CXX=icpc cmake -DMKL=ON .. `
* Intel machine with processor name supported by `-ax` option (e.g. Intel(R) Xeon(R) Platinum 8480L is SAPPHIRERAPIDS) 
  (Intel LLVM Compilers): `FC=ifx CC=icx CXX=icpx cmake -DMKL=ON -DINTEL_CODE_NAME=SAPPHIRERAPIDS .. `
* Intel machine with processor name not supported by `-ax` option (Intel LLVM Compilers): `FC=ifx CC=icx CXX=icpx cmake -DMKL=ON ..`
* AMD-based machines with `openBLAS` and `scaLAPACK`: `FC=gfortran CC=gcc CXX=gcc cmake -DOPENBLAS=ON -DSCALAPACK=ON ..`
* AMD-based machine with `AOCL-FFTW` and `BLIS`+`libFLAME`: `FC=gfortran CC=gcc CXX=gcc cmake -DAMDLINALG=ON ..`

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
on its [Github page](https://github.com/Fortran-FOSS-Programmers/ford), and additional details regarding
installation of dependencies can be found under 'Known Issues', below. 


Compiling with exciting with SIRIUS
-----------------------------------

[SIRIUS](https://github.com/electronic-structure/SIRIUS) is a domain specific library for electronic 
structure calculations. It implements pseudopotential plane wave (PP-PW) and full potential linearized 
augmented plane wave (FP-LAPW) methods, and is designed for GPU acceleration of popular community codes 
such as Exciting, Elk and Quantum ESPRESSO.

Compiling exciting with SIRIUS is complex. To simplify the procecss of building dependencies, SIRIUS can be
completely installed with the python package manager [spack](https://spack.readthedocs.io/en/latest/getting_started.html).

As of exciting Neon, only a CPU build chain using GCC on Ubuntu Focal is regularly tested in exciting's 
CI. This is provided in [build/utilities/docker/Dockerfile_ci_sirius](build/utilities/docker/Dockerfile_ci_sirius). 

```shell
# Install sirius dependencies with spack
# ---------------------------------------
# Install spack: https://spack.readthedocs.io/en/latest/index.html

# Find any preinstalled compilers or dependencies
# Note, it's preferable to use a preinstalled compiler were possible, to minimise total build time
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

SIRIUS Gotchas
------------------

* `spack` command is not added to the `$PATH` by default. To source `spack` in bash, type 
  `source spack/share/spack/setup-env.sh`

* sirius is under active development, and its dependencies change regularly. Please use `spack info sirius`
  for the current status

* spack does not query  a server when `spack info` is called. To ensure the info is up-to-date, run
  `git pull` in the root of the spack directory.

* On some architectures, sirius will install to `$SIRIUS_ROOT/lib64`, not `$SIRIUS_ROOT/lib`. This will be
  clear at the linking step, where exciting fails to find sirius, and requires one to manually edit the
  `make.inc` file to point to the correct directory.

fastBSE
------------------
To use fastBSE, you need to link it with HDF5 and FFTW3.

**GCC**

To install the libraries in a Debian-based distribution
```shell
  sudo apt-get install -y libhdf5-serial-dev libhdf5-mpi-dev libfftw3-dev
```

Then, to compile `exciting`
```bash
  mkdir build
  cd build
  FC=gfortran CC=gcc CXX=gcc ../external/cmake-3.31.3-linux-x86_64/bin/cmake -DHDF5=ON ..
  make -j N -l N exciting_mpismp
  make install
```

**INTEL**

If you intend to use Intel, you must manually compile HDF5. We recommend using HDF5 1.12.0. The source code can be found on the HDF5 website. 
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


Compiler Support
------------------

exciting requires an Fortran2018-compliant compiler. exciting is known to compile with:

* Intel ifort: 2021

* Intel ifx: 2025.0.0
  
* GNU gfortran: 10, 12, 14, 15
  
* Crayftn: 18.0.1

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

