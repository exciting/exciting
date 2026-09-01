In this tutorial you will lear how to compile of **`exciting`** based on your needs. Further technical details can be found in the `INSTALL.md` file in the `exciting` source directory.

## Table of Contents

1. [Basic Compilation](#basic-compilation)  
   - 1.1. [Requirements](#requirements)  
   - 1.2. [Basic compilation using GNU compilers (gfortran)](#basic-compilation-using-gnu-compilers-gfortran)  
   - 1.3. [Basic compilation using Intel compilers (ifx)](#basic-compilation-using-intel-compilers-ifx)  

2. [How to Tune Your Compile Options for Better Performance](#how-to-tune-your-compile-options-for-better-performance)  
   - 2.1. [Selecting a performant library](#selecting-a-performant-library)  
     - 2.1.1. [Linear Algebra Libraries](#linear-algebra-libraries)  
     - 2.1.2. [FFT Libraries](#fft-libraries)  
   - 2.2. [Passing Extra Options to the Fortran Compiler](#passing-extra-options-to-the-fortran-compiler)  
   - 2.3. [Extra Options for Further Tuning with the Modern Intel Compiler (IFX) for Intel Processors](#extra-options-for-further-tuning-with-the-modern-intel-compiler-ifx-for-intel-processors)  

3. [Advanced: Other Production Options](#advanced-other-production-options)  
   - 3.1. [HDF5](#hdf5)  
   - 3.2. [ScaLAPACK](#scalapack)  
   - 3.3. [SIRIUS](#sirius)  

4. [Advanced: GPU Compilation](#advanced-gpu-compilation)  
   - 4.1. [Compiling for AMD and NVIDIA GPUs using `Crayftn`](#compiling-for-amd-and-nvidia-gpus-using-crayftn)  
   - 4.2. [Compiling for Intel GPUs using `ifx`](#compiling-for-intel-gpus-using-ifx)  
     - 4.2.1. [Unified Shared Memory Support](#unified-shared-memory-support)  

## 1. Basic Compilation

**`exciting`** uses **CMake** (version 3.25 or later) for building and installation.
You can download a precompiled version of CMake suitable for your system from [here](https://cmake.org/download/).
For **x86-64 Linux**, bundled CMake 3.31.3 binaries are provided with **exciting** under `external/cmake-3.31.3-linux-x86_64/`.

### 1.1 Requirements

**`exciting`** requires `xsltproc` to preprocess its XML schema into code.

Additionally, the code depends on the installation of:

- **FFTW3** implementations (such as `oneMKL`, `AOCL-FFTW`, `FFTW3`, `Cray-FFTW`, etc.)
- A **BLAS/LAPACK** implementation (such as `oneMKL`, `BLIS` + `libFLAME`, `OpenBLAS`, `LibSci`, etc.)

> **Please ensure you have these libraries and binaries installed before proceeding.**

> **Note:** For performant compilation, the choice of these libraries is essential.

The **exciting** package requires the following external libraries, which are either included or downloaded during compilation:

- [FoX XML](https://github.com/andreww/fox) — Library for parsing input (2012 version).
- [LIBXC V7](https://libxc.gitlab.io/) — Library of DFT exchange and correlation functionals. Version 7.0.0 is downloaded from https://gitlab.com/libxc/libxc/ during the build; external libXC versions 5.0.0 or higher are also supported.
- [BSPLINE-FORTRAN](https://github.com/jacobwilliams/bspline-fortran) — Multidimensional B-spline interpolation on a regular grid.

> **Note:** For performant compilations supporting fully parallel execution, an MPI library is required, such as Open MPI, MPICH, or Intel MPI. Optionally, a version of ScaLAPACK can be used.
> These dependencies can be installed via package managers like APT, Conda, Spack, EasyBuild, containerized environments, or built manually from source.

Finally, C, C++, and Fortran compilers are required. In particular, a Fortran compiler supporting the Fortran 2018 standard is needed. Currently, the code supports:

- Intel Classic (ifort): 2021.13.0
- Intel LLVM (ifx): 2025.0.0
- GNU: versions 12, 14, 15
- Cray: 18.0.1, 19.0.1
- LLVM-Flang-based compilers (**must support Fortran 2018 standard**)

### 1.2 Basic Compilation Using GNU Compilers (`gfortran`)

This subsection shows how to perform a basic and **non-performant** compilation of **exciting** using `GNU` compilers.

The default build/compilation process attempts to automatically create a performant binary with CPU-level parallelization enabled. However, several options exist to further tune the build for improved performance.

In this example, we will compile **exciting** with limited parallelization (only multithreading), using the default libraries which may not be optimized for performance.

```bash
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=gcc -DMPI=OFF ..
make -j
make install
```

This will install **`exciting`** to: `$EXCITINGROOT/install/bin/exciting`

### 1.3 Basic Compilation Using Intel Compilers (`ifx`)

This subsection shows how to perform a basic and **non-performant** compilation of **exciting** using modern `Intel` compilers.

The default build/compilation process attempts to automatically create a performant binary with CPU-level parallelization enabled. However, several options exist to further tune the build for improved performance.

In this example, we will compile **exciting** with limited parallelization (only multithreading), using the `oneAPI` libraries which should offer a relatively good performance for Intel machines.

```bash
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DMPI=OFF -DMKL=ON ..
make -j
make install
```

This will install **`exciting`** to: `$EXCITINGROOT/install/bin/exciting`

> **Note** recall to load the appropriate modules, e.g. `module load intel-oneapi/2025.1` as in most systems Intel oneAPI is not loaded by default.

## 2 How to tune your compile options for a better performance

As mentioned in the previous section, `cmake` will by default attempt to compile a performant version of **exciting** tailored to the architecture of the compilation machine.

However, there are cases where manual tuning can further improve performance; either by selecting specific high-performance libraries or by adding additional compilation options.

This tuning is especially important for HPC systems, where the login node or build environment may differ architecturally from the compute nodes on which **exciting** will run.

### 2.1 Selecting a performant library:

The `exciting` code, like most Density Functional Theory (DFT) codes, relies heavily on linear algebra and fast Fourier transforms (FFT), which are provided by external libraries. Therefore, the overall performance of the code depends significantly on the efficiency and parallelization schemes of these libraries.

For optimal compilation performance, we recommend using multithreaded versions of these libraries and, whenever possible, vendor-specific libraries tailored for the target architecture.

#### 2.1.1 Linear Algebra Libraries

The choice of linear algebra libraries can be controlled through the [`BLAVENDOR`](https://cmake.org/cmake/help/latest/module/FindBLAS.html) CMake option. Additionally, **exciting** provides specific tested options that automatically select appropriately parallelized library versions based on other compilation settings.

You can tune the linear algebra backend using the following CMake options (**Note that options in `cmake` are written as follows `-DOPTION=OPTION_VALUE`):

- **_MKL_**: Uses Intel MKL for linear algebra (default: OFF).
  _Recommended for Intel machines._
- **_OPENBLAS_**: Uses OpenBLAS for linear algebra (default: OFF).
- **_AMDLINALG_**: Uses AMD linear algebra libraries (BLIS and FLAME) (default: OFF).
  _Recommended for AMD machines._
- **_CRAYLIBSCI_**: Uses Cray LibSci for linear algebra and, if required, ScaLAPACK (default: OFF).
- **_OTHERLINALG_**: Uses another linear algebra library (not officially supported) (default: OFF).
- **_LINALGLIB_**: If `OTHERLINALG` is ON, specify the full path to the desired linear algebra libraries (default: None).

#### 2.1.2 FFT Libraries

`exciting` requires an FFTW3 implementation. Similar to linear algebra libraries, we recommend using multithreaded and vendor-optimized FFT libraries when available.

Recommended FFT implementations based on target architecture:

- **AOCL-FFTW**: Recommended for AMD machines
- **Cray-FFTW**: Recommended for HPE HPC systems
- **Intel MKL FFT**: Recommended for Intel machines

The FFT library used can be controlled via the module system or the CMake option:

- **_FFTW3_ROOT_**: For non-standard compilation, specify the installation directory of FFTW3 (default: None).

> If you're using **Intel MKL**, the build system will use the FFTW3 implementation included in this library.

#### 2.2 Passing Extra Options to the Fortran Compiler

If needed, extra options can be passed to the Fortran compiler via the CMake option:

```bash
cmake [...] -DCMAKE_Fortran_FLAGS="extra options"
```
> [...] indicates other possible options

This allows, for instance, changing the target architecture for which the compilation is performed, which by default is the architecture of the machine where the build takes place. 

For instance, to lower the optimization level from the default:

```bash
cmake [...] -DCMAKE_Fortran_FLAGS="-O1"
```

##### 2.3 Extra Options for Further Tuning with the Modern Intel Compiler (IFX) for Intel Processors

Obtaining a performant binary using `IFX` is more nuanced than with classical Intel compilers. In particular, the optimization option `-O3` does not enable Intel-specific optimizations as it does in the classical compilers, but only generic optimizations.

To enable Intel-specific optimizations, `IFX` requires the option:
```bash
ifx [...] -x<code>
```
where `<code>` is the [Intel Processor ID](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2025-0/x-qx.html), and `[...]` indicates other possible compiler options.

By default, the exciting build system sets <code> to Host, which generates optimizations for the CPU architecture of the machine where the build/compilation takes place.

If this does not match the architecture of the machine where the code will run, you can modify it with the CMake option:

```bash
cmake [..] -DINTEL_CODE_NAME=ID
```
where `ID` is the Intel Processor ID of the target machine.

> **Note:** Do not use binaries compiled with the Intel LLVM Fortran compiler (IFX) for BSE (Bethe-Salpeter Equation) calculations. IFX is known to produce lower numerical precision in `exciting` for BSE, which can lead to unreliable results—especially in advanced cases such as dichroism or spin-polarized BSE. In those cases, we recommend using other compilers (e.g., IFORT, GFortran) for better accuracy and stability.

## 3. Advanced: Other Production Options

In addition to the previously discussed options, `exciting` can be compiled with support for additional libraries. This section covers those options.

### HDF5
[HDF5](https://www.hdfgroup.org/solutions/hdf5/) is a data model and file format designed to store and organize large, complex scientific data efficiently and flexibly.

This library is mandatory for both [`fastBSE`](https://exciting-code.org/uploads/exciting/tutorial_notebooks/optical_spectra_from_fastBSE.html) and [RIXS](!!!!TODO: link Elias RIXS tutorial).

To enable it, use the `cmake` `HDF5` option as follows:

```bash
cmake [...] -DHDF5=ON ..
```
> [...] indicates other possible options

### ScaLAPACK

[ScaLAPACK](https://www.netlib.org/scalapack/) is a library of high-performance linear algebra routines designed for parallel distributed-memory machines.

[BSE](https://exciting-code.org/uploads/exciting/tutorial_notebooks/excited_states_from_bse.html) calculations involve the diagonalization of relatively large matrices, so using ScaLAPACK is strongly recommended for efficient and scalable performance.

To do so, use the `cmake` `SCALAPACK` option as follows:

```bash
cmake [...] -DSCALAPACK=ON ..
```
> **Note** in case your ScaLAPACK installation is not in a standard path, you can specify the installation directory as `cmake [...] -DSCALAPACK=ON -DSCALAPACK_ROOT="pat/to/ScaLAPACK/installation/directory"`

> If you're using **Intel MKL** or **Cray LibSci**, and set `SCALAPACK=ON`, the build system will use the ScaLAPACK implementation included in those libraries.

> [...] indicates other possible options

### SIRIUS

[SIRIUS](https://github.com/electronic-structure/SIRIUS) is a domain-specific library for electronic structure calculations. For ground-state calculations, it offers higher levels of parallelization than vanilla `exciting`, making it suitable for extremely large systems.

> **Note:** **SIRIUS** can only be used for ground-state calculations. The resulting ground-state cannot be used as a starting point for any other calculation beyond ground-state (e.g., GW, BSE etc).

To enable **SIRIUS**, use the following `cmake` options:

```bash
cmake [..] -DSIRIUS=ON -DUSE_INTERNAL_LIBXC=OFF ..
```

> [...] indicates other possible options

The `USE_INTERNAL_LIBXC=OFF` option forces `exciting` to link against the same `libXC` version that **SIRIUS** was compiled with. Otherwise, runtime errors may occur.

> **Note** in case your SIRIUS installation is not in a standard path, you can specify the installation directory as `cmake [...] -DSIRIUS=ON  -DSIRIUS_ROOT="path/to/SIRIUS/installation/directory" -DUSE_INTERNAL_LIBXC`


## 4. Advanced: GPU Compilation

`exciting` can benefit from GPU offloading for the following tasks:
- [G₀W₀](https://exciting-code.org/uploads/exciting/tutorial_notebooks/electronic_band_structure_from_gw.html) calculations
- Ground-state calculations when using [hybrid](https://exciting-code.org/uploads/exciting/tutorial_notebooks/hybrid_functional_calculations.html) functionals

`exciting` uses a vendor-independent approach based on OpenMP offload (version 5.2 or higher) for loops and memory management, combined with calls to wrappers for vendor-specific routines handling linear algebra and fast Fourier transforms.

> **Note:** If GPU compilation is enabled, `cmake` will check for OpenMP offload support. If this test fails, compilation will be aborted.

### 4.1 Compiling for AMD and NVIDIA GPUs using `Crayftn`

To compile GPU-accelerated `exciting` for NVIDIA and AMD GPUs with the Cray Fortran compiler, you need version 18.0.1 or higher. In both cases, you must link against [`MAGMA`](https://icl.utk.edu/magma/) version 2.7.2 or higher, compiled against a `BLAS/LAPACK` library also compiled with the Cray compiler.
**We recommend using OpenBLAS, as LibSci has shown instability in some cases with multithreading.**

> **Compiling OpenBLAS with the Cray compiler:** It is important to include the `-hnopattern` flag for the Fortran compiler, disable tests during build, and use `crayftn` and `craycc` directly instead of the Cray wrappers for this specific case.

> **Note:** Before running the configuration/compilation commands on an HPE system, you need to load the appropriate accelerator module, for example:
> `module load craype-accel-amd-gfx90a`

For AMD cards, use the following `cmake` options:

```bash
export FC=ftn
export CC=cc
export CXX=cc
cmake [...] -DAMD=ON -DAMDTARGET=<AMD-target-id> -DMAGMA_ROOT=<path/to/MAGMA/install/directory> ..
```

where `AMDTARGET` is the [AMD-target-id](https://llvm.org/docs/AMDGPUUsage.html#target-triples), and `[...]` indicates other CMake options.

For NVIDIA cards, use:
```bash
export FC=ftn
export CC=cc
export CXX=cc
cmake [...] -DNVIDIA=ON -DNVIDIAARC=<NVIDIA-architecture> -DMAGMA_ROOT=<path/to/MAGMA/install/directory> ..
```
where `NVIDIAARC` is the NVIDIA architecture version without the `sm_` prefix (e.g., `89` for compute capability 8.9).

### 4.2.1 Unified Shared Memory Support

`exciting` supports offloading to architectures with unified shared memory (USM) via `Crayftn`, particularly AMD Accelerated Processing Units (APUs), where the CPU and GPU share the same physical memory. This eliminates the need for explicit memory transfer or management.

To enable this, add the `-DUSM=ON` option to your `cmake` configuration step:

```bash
cmake [...] -DUSM=ON ..
```

> [...] indicates other possible options

This instructs the compiler to add the appropriate flags and compile versions of the functions compatible with unified shared memory.

### 4.2 Compiling for Intel GPUs using `ifx`

To compile GPU-accelerated `exciting` for Intel GPUs, you need the `ifx` compiler version `2025.0` or higher. Additionally, enabling `MKL` support is mandatory. Use the following `cmake` options:
```bash
cmake [...] -DMKL=ON -DINTEL=ON ..
```
> [...] indicates other possible options

