# For Developers

Documentation and points of discussion for GreenX developers.

## Portability

* The minimax library has been tested on:
  - macOS Catalina, with GCC9 and openBLAS
  - Ubuntu-20.04, with GCC9 and blas/lapack
  
* There has been no testing with Intel or MKL. This should be done.

## Testing 

The test framework is currently set up such that one:
* Writes a .f90 test driver, which calls some GreenX API
* Writes a .py file to: 
    * Supply inputs to the .f90 driver
    * Execute the driver binary
    * Parse the outputs of the test
    * Make assertions
    
This has some requirements, and some technical open questions associated with it.
Two better approaches would be to:  

a). Replace a fortran unit testing framework. In most cases, this is more appropriate.
     This would remove the need for all I/O and path-setting.
     
b). Call the respective fortran library directly from python, removing the need
    for the .f90 driver, and removing all the complications associated with
    consistent binary/file paths.

### Issues with the Test Set Up

#### Test Binary Name

* The fortran binary name is specified in the python test.
  - If one enforces the convention that the .f90 and .py files must have the same
    prefix, the binary name can automatically be retrieved by python
    `os.path.basename(__file__).py`

#### Choice of Test Framework

* One should test the library by directly calling from python, however Alex 
  (thinks) memory allocation in the library prevents this.
  - Refactor the library, such that the caller handles the memory. 

* Using a fortran test framework OR calling the library from python removes
  all of the issues with the test set-up.