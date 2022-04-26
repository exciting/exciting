# Writing Unit Tests in exciting

## Running Existing Unit Tests

exciting has its own unit-testing framework which comprises a `unit_test_type` for performing assertions, and a 
series of test drivers for running tests. Unit tests can be run by typing:

```bash 
bin/exciting_serial -run-unit-tests=all
```

or 

```bash 
mpirun -np 2 bin/exciting_mpismp -run-unit-tests=all
```

for serial and parallel binaries, respectively (remember to also set threads appropriately). To run a specific unit test, 
or mutliple unit tests, one can type:

```bash
bin/exciting_serial -run-unit-tests= name_of_test, name_of_another_test
```

where the test names must be comma-separated. The command line parsing has been written such that it is insensitive to 
whitespaces between options (-run-unit-tests) and sub-options (name_of_test). 

If an assertion fails, the default behaviour is for all tests to continue running. To change this, use the 
`kill-on-failure` option:

```bash
bin/exciting_serial -run-unit-tests=all -kill-on-failure
```

## Adding to Unit Test Sub-Options

It is proposed that the command line sub-options approximately correspond to directories in `src`. For example, we might 
expect:

```bash
bin/exciting_serial -run-unit-tests= gw, math, hybrids
```

runs all tests in the directories: `math`, `src_gw` and `src_hyrbids`. If a sub-option is not recognised, that implies
no tests exist for that directory. The developer is therefore required to add some infrastructure to enable it. 

To add the sub-option to run unit tests in the `hybrids` directory, for example, one needs to first extend the 
`unit_tests_type` type and `subroutine set_unit_test` in `src/input/cmd_line_args.f90`:

```fortran
   !> Each logical is initialised to .false.
   type unit_tests_type
      logical :: all = .false.    !! All unit tests
      logical :: math = .false.   !! tests in the math directory
      logical :: gw = .false.     !! GW
      logical :: hybrids = .false. !! For tests in the src_hybrids directory  
   contains
      procedure :: init => set_unit_tests
   end type unit_tests_type
```

```fortran
  subroutine set_unit_test(run, test_name, mpi_env)

      select case (test_name)
      case ('all')
         run%all = .true.
      case ('gw')
         run%gw = .true.
      case ('math')
         run%math = .true.
      ! Add the new module for testing here 
      case ('hyrbids')
         run%hybrids = .true.
```

where `set_unit_test` essentially maps the sub-option string to the corresponding logical. Note that `select case` can 
take more than one argument, i.e. `(‘tddft’, ‘TDDFT’)`although the suggestion is to stick to lowercase.

One is then required to add a call to the test driver for the `src_hybrids` directory to 
`src/testframework/unit_test_drivers.f90`:

```fortran
module unit_test_drivers

contains 

   subroutine unit_test_driver(kill_on_failure)
     ! ... some code
   
        if (run%math .or. run%all) then
            call math_test_driver(mpiglobal, kill_on_failure) 
        end if
       
        if (run%gw .or. run%all) then
            call gw_test_driver(mpiglobal, kill_on_failure)
        end if
    
       ! Test driver call goes here
        if (run%hyrbids .or. run%all) then
            call hybrids_test_driver(mpiglobal, kill_on_failure)
        end if
   
   end subroutine unit_test_driver

end module unit_test_drivers
```

Finally, one writes `hybrids_test_driver.f90` in `src_hybrids`:

```fortran
!> Module for managing calls to all test drivers within src_hyrbids 
module hyrbids_test_drivers
  use modmpi, only: mpiinfo
  ! Load test drivers here. One test driver per module of the directory
  use some_module, only: some_module_test_driver
  use another_module, only: another_module_test_driver

  private
  public :: hyrbids_test_driver

contains

  !> Calls test drivers for each module of the directory
  subroutine hybrids_test_driver(mpiglobal, kill_on_failure)
    type(mpiinfo), intent(in) :: mpiglobal
    logical, optional :: kill_on_failure 
    ! Call test drivers for each module here
    call calc_vnlmat_test_driver(mpiglobal, kill_on_failure)
    call expint_test_driver(mpiglobal, kill_on_failure)
    !... etc
  end subroutine hybrids_test_driver

end module 
```

## Writing a Test Driver for a Module

It is proposed that each directory of the exciting source has one module that collates the calls to all tests within
the directory. Such as `hyrbids_test_drivers`, shown above. Each module of the directory will subsequently have its
own test module.

For example, `calc_vnlmat.f90` will have the test module `calc_vnlmat_tests.f90`, which contains `calc_vnlmat_test_driver`. 
`expint.f90` will have the test module `expint_tests.f90`, which contains `expint_test_driver`. Both of these test drivers
are called by `hybrids_test_driver`.

In general a `module_tests.f90` will contain a test_driver and a series of tests. Schematically, for `calc_vnlmat.f90`
this would look like:

```fortran

module calc_vnlmat_tests

   public :: calc_vnlmat_test_driver
contains 

   subroutine calc_vnlmat_test_driver(mpiglobal, kill_on_failure)
      type(mpiinfo), intent(in) :: mpiglobal
      logical, intent(in), optional :: kill_on_failure

      !> Test report object
      type(unit_test_type) :: test
      !> Number of assertions
      integer, parameter :: n_assertions = 19

      call test_report%init(n_assertions, mpiglobal)
      call test_vnlmat(test)
      
      call test_report%finalise()
   end subroutine calc_vnlmat_test_driver
   
   
   subroutine test_vnlmat(test)
     class(unit_test_type), intent(inout) :: test
     ! Define some reference data for the expected H
     H_ref = []
     ! Perform some sensible assertions
     call calc_vnlmat(some, args, H)
     test%assert(all_close(H, H_ref), message='Expect H = H_ref')
   end subroutine test_vnlmat

end module
```
 
 For each subroutine of the module, one can make a series of assertions on the data it returns. In many cases writing
 unit tests might be difficult, and we may need to consider writing routines to facilitate testing, such as ones that
 return a basis, or a set of crystal positions, for example. Additionally, we can also perform 
 [mocking](https://stackoverflow.com/questions/2665812/what-is-mocking) of input arguments. 
 
 Unfortunately, the number of assertions needs to be manually specified, and therefore updated as new assertions are
 added. Without a linked list implementation, this is a constraint of fortran and the current test framework 
 implementation. The test suite will warn you if the number is too small. 
 
 The `test` object and its `assert` method have been written to behave like the test object in the [Zofu unit testing
 framework](https://github.com/acroucher/zofu). The longer-term plan is to switch out the current test framework in 
 favour of Zofu, once the build system has been overhauled.  
 
 Threaded, and in principle MPI code, are also unit-testable. 
