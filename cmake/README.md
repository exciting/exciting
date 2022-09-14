# CMake for exciting

CMake for exciting is a *work-in-progress*, and cannot be used to build exciting.

## Outline of TODOs:

* Check for essential packages: `xsltproc`
* Building external libraries (start by building those in `external/`, extend to cloning if not present)
  * FoX
  * Lbfgsb.3.0
  * leblaiklib
  * libXC - see `cmake/libxc.cmake`
* Generate version stamp, calling or replacing `versionstamp.pl`
* Custom command to preprocess `inputmodules.f90` and `speciesmodules.f90` from XML schema
* Find Ford and `make docs`
* Build exciting as a lib
* Link lib to exciting binary
* Migrate to a unit-test lib
* Migrate regression-testing to CMake
* Introduce unit tests in CMake/ctest

## For Developers

Create a build directory in exciting's root (noting build is already taken):

```bash
mkdir build_dir && cd build_dir
# Configure command
cmake ../
# Build command
make -j
```
