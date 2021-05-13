# Mac OS Catalina Compilation

## GCC 9 with openBLAS Macports Installation

### Macport Dependencies

  To obtain the necessary compiler and libraries:

  ```
  sudo port install gcc9 openmpi-gcc9 OpenBLAS
  ```

### Modifying the make.inc

  The include path to the OMP header should be added to the INC variable:

  ```
  INC = -I/opt/local/lib/gcc9/gcc/x86_64-apple-darwin19/9.3.0/finclude
  ```

  and linking should be set to OpenBLAS:

  ```
  LIB_LPK = -L./ -llapack -lopenblas
  ```

  where ` -L./` has been retained in the library search path for the bundled libraries.


## GCC 10 with openBLAS Macports Installation

### Macport Dependencies

  To obtain the necessary compiler and libraries:

  ```
  sudo port install gcc10 openmpi-gcc10 OpenBLAS
  ```

## Modifying the make.inc

  GCC 10 is stricter with regards to argument mismatches. To switch the error to a warning,
  append the fortran compilation flags F90_OPTS and F90_DEBUGOPTS, with:

  ```
  -fallow-argument-mismatch
  ```

  The include path to the OMP header should be added to the INC variable:

  ```
  INC = -I/opt/local/lib/gcc10/gcc/x86_64-apple-darwin19/10.2.0/finclude/
  ```

  and linking should be set to OpenBLAS:

  ```
  LIB_LPK = -L./ -llapack -lopenblas
  ```

  where ` -L./` has been retained in the library search path for the bundled libraries.


## GCC 9 with Homebrew Installation

Installation of exciting has not yet been attempted with external dependencies installed using
homebrew. If you have a working solution, please open a pull/merge request. 
