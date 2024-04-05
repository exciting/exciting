# Compiling exciting


Requirements
------------------
exciting requires `xsltproc` to preprocess its XML schema into code. **PLEASE** ensure you have
this installed before proceeding.

exciting comes with all the external libraries required to compile the code in serial.
These include:

* [FoX XML](http://www.uszla.me.uk/FoX/DoX/Licensing.html) library for parsing the input

* [BLAS & LAPACK](http://www.netlib.org/lapack/LICENSE) for linear algebra operations 
		
* [LIBXC V2](https://www.tddft.org/programs/libxc/) library of DFT exchange and correlation functionals 
	
* [BSPLINE-FORTRAN](https://github.com/jacobwilliams/bspline-fortran) Multidimensional B-Spline interpolation of data on a regular grid

Compilation for fully parallel execution requires an MPI library, such as open MPI, MPICH or Intel MPI library and optionally, a version of Scalapack. These can be installed with a package manager such as 
APT, Conda, Spack or EasyBuild, or built manually from source. Some limited spack recipes for 
installing external libraries, compiled with GCC and Intel, are provided in the [repository](build/utilities/spack). 

Test suite dependencies are specified in [test/README](test/README). 


Compiling
------------------

In order to compile the code, one is required to copy a `make.inc` file from the examples:

```shell
  cp build/platforms/make.inc.ifort build/make.inc
```   

which copies the `make.inc` for Intel. Examples are provided for the Intel and GCC compilers.
The code can then be built by typing:

```shell
  make 
```

in the exciting root directory. By default, the *mpiandsmp* and *smp* versions of exciting
are built, where *mpiandsmp* refers to a binary with both MPI and openMP parallelism enabled.
Serial and pure MPI builds are also possible by typing:

```shell
  make serial
```

and

```shell
  make mpi
```

respectively. Debug builds are available for the *serial* and *mpiandsmp* builds, by typing:

```shell
  make debug
```

or

```shell
  make debugmpiandsmp
```

respectively. The test suite can be run by typing:

```shell
  make test
```

For a full list of build options, one may consult the [Makefile](./Makefile) in the root directory. 
Auxiliary programs include "spacegroup" for producing crystal geometries
from spacegroup data, "species" for generating species files, and "eos" for
fitting equations of state to energy-volume data.

Please note that in order to use: 

```shell
  make clean
```

a `make.inc` file must be present in the build directory. 


CMake Support
------------------

Whilst CMake is present in exciting, this is in the developmental phase and CANNOT
be used to compile the code.


OS Support 
------------------

### Mac OS

exciting can be compiled on mac OS, but it is complicated by nonstandard installation
locations for libraries, which can vary between OS versions and package managers.

Because production calculations typically require HPC resources, compilation on mac is
primarily intended for developers, and as such, we only provide schematic instructions
for modifying the `make.inc`. Details are provided in [build/platforms/os_mac](build/platforms/os_mac)

For mac, exciting has been tested on:

* OS Catalina using GCC 9 and 10, openBLAS 0.3.12 and openmpi 4.1.0 

* OS Monterey using GCC 12.2.0, openBLAS 0.3.21 and openmpi 4.1.4


Documentation
------------------

exciting's input documentation and the majority of its source documentation can be built with:

```shell
  make doc 
```

PDF documents are built in their respective `doc/` subfolders. All documentation is parsed
through LaTex and requires the user to have a Tex distribution installed. 

exciting's source code documentation is in the process of migrating to FORD. At present, only
select excited-state modules are documented. These are specified in [docs/ford_settings.md](docs/ford_settings.md) 
under the `src_dir` tag. 

FORD is available as a python package and can be installed with pip. To install FORD, type:

```shell
  pip3 install ford
```

To build exciting's documentation using FORD, type:

```shell
  make ford
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

As of exciting Flurine, only a CPU build chain using GCC on Ubuntu Focal is regularly tested in exciting's 
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
export SPEC="sirius@develop %gcc@9.4.0:12.2.0 build_type=Release +scalapack +fortran ^mpich@3.3.2 ^intel-oneapi-mkl+cluster ^spfft target=x86_64"

# Install sirius's dependencies
spack install --fresh --only=dependencies $SPEC

# Define these directories where you want to:
#  Install sirius and all corresponding dependencies
export SIRIUS_ROOT=spack-sirius-install
#  Install sirius environment details (yaml file)
export SIRIUS_ENV=spack-sirius-env

# Create a spack environment for sirius
spack env create --with-view $SIRIUS_ROOT $SIRIUS_ENV
spack -e $SIRIUS_ENV add $SPEC

# Install sirius into the defined env
spack -e $SIRIUS_ENV install $SPEC
```

exciting can then be built, linking to SIRIUS, as follows from the project root:

```shell
# Copy the appropriate make.inc
cp build/platforms/make.inc.sirius build/make.inc

# Ensure `SIRIUS_ROOT` is edited in the `build/make.inc`
# to be consistent with the user's own installation path:
# emacs build/make.inc
# SIRIUS_ROOT=<ADD ME>

spack build-env ${SPEC} -- make mpiandsmp -j 4
```

and tested with:

```shell
cd test/
# Ensure the bottom of `test/defaults_config.yml` is set to:
# group_execution:
#   NONE: False
#   SLOW_TESTS: False
#   GW: False
#   LIBXC: False
#   SIRIUS: True
python3 runtest.py -e exciting_mpismp -t sirius
```


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
  cd $EXCITINGROOT
  cp build/platforms/make.inc.gfortran10+.hdf5.fftw3 build/make.inc
  make smp # Choose the required build type as explained above.
```
Please note that the recommended `make.inc` file supports **GCC10+**. If you use an older version drop the flags `-fallow-argument-mismatch` and `-fallow-invalid-boz`.

**INTEL**

If you intend to use Intel, you must manually compile HDF5. We recommend using HDF5 1.12.0. The source code can be found on the HDF5 website. 
We recommend configuring the build as follows:

```shell
  configure --enable-build-mode=production --enable-parallel --enable-fortran --enable-fortran2003
```

For details on the options, refer to the HDF5 documentation. 
After a successful installation, update the PATH variable with the location of the HDF5 binaries:

```shell
  export PATH="path/to/hdf5-build/hdf5/bin:$PATH"
```

Intel provides an FFTW3 API to their own FFT library in the MKL library. We recommend using this API.

Then follow these steps:

```shell
  cp build/platforms/make.inc.intel.hdf5.fftw3-mkl-api build/make.inc
  make smp # Choose the required build type as explained above.
```

Please note that this code is not executable as-is; it provides instructions for setting up the fastBSE environment. 
Make sure to adjust the paths and options according to your system and requirements.


Compiler Support
------------------

exciting requires an F2008-compliant compiler. exciting is known to compile with:

* Intel ifort: 2019, 2021
  
* GCC gfortran: 8, 9, 10, 12
  

Compliant but not tested:

* Intel: 2022

* GCC: 11


Known Issues
------------------

### Intel MPI

Intel MPI 2019 and 2021.3 are known to exhibit memory leaks. This can be troublesome for memory-intensive calculations
and cause exciting to crash. See the threads associated with [Qbox](https://groups.google.com/g/cp2k/c/BJ9c21ey0Ls) and
[CP2K](https://github.com/cp2k/cp2k/issues/1830) for further discussion.


### Intel Compiler 2018

Compilation is known to fail at `src_lib/mod_manopt.f90` with the error:

```
   catastrophic error: **Internal compiler error: segmentation violation signal raised** 
```

This is likely caused by a bug in the 2018 version of the ifort compiler.

### GCC 10 and Above

If compiling with GCC 10, one should append the F90_OPTS, F77_OPTS and F90_DEBUGOPTS variables in
make.inc with `-fallow-argument-mismatch`. See:

  https://gcc.gnu.org/gcc-10/porting_to.html

for more details.

### FORD Fails to Find Graphviz

On Debian Buster, we have noticed that FORD fails to find graphviz, despite pip installing
it as part of FORD's dependencies. This can be fixed by installing graphviz with APT. Optionally
one can also install python-dev and LXML, which facilitate faster generation of the search database:

```shell
  sudo apt install python3-all-dev graphviz
  pip3 install lxml
```
