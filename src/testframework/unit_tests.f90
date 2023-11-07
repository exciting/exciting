!> Unit tests type, indicating which unit tests to run.
!> One notes that a change of build system to CMake, where each
!> unit test could be defined as its own main program, would
!> remove the need for this module.
module unit_tests
   use modmpi, only: mpiinfo, terminate_mpi_env
   implicit none
   private
   public :: unit_tests_type

   !> Type containing logicals for all unit test modules in exciting.
   !> Used to indicate which unit tests to run.
   !> Unit test modules are defined as one per
   !>
   !> All unit tests initialised to .false.
   type unit_tests_type
      !> All unit tests
      logical :: all = .false.

      logical :: advanced = .false.
      logical :: eigensystem = .false.
      logical :: fermisurfdx = .false.
      logical :: gw = .false.
      logical :: hybrids = .false.
      logical :: LDAU = .false.
      logical :: mixing = .false.
      logical :: mpi = .false.
      logical :: optics = .false.
      logical :: phonon = .false.
      logical :: raman = .false.
      logical :: rdmft = .false.
      logical :: stm = .false.
      logical :: sym = .false.
      logical :: vdw = .false.
      logical :: wannier = .false.
      logical :: xc = .false.
      logical :: xs = .false.

      logical :: xstring = .false.
      logical :: lapack = .false.
      logical :: math = .false.
      logical :: simplified_input = .false.
      logical :: structure = .false.
      logical :: testframework = .false.
      logical :: file_io = .false.
      logical :: matrix_elements = .false.
      logical :: xgrid = .false.
      logical :: xhdf5 = .false.

   contains
      procedure :: init => set_unit_tests
   end type unit_tests_type

contains

   !> Initialise object of type 'unit_tests_type'.
   !>
   !> Given a unit_tests_str of the form:
   !>
   !>    "gw,math,dft,tddft"
   !> or
   !>    "all"
   !>
   !> sequentially get each unit test name and use
   !> it to set its corresponding logical in an
   !> instance of unit_tests_type
   !>
   !> Expects:
   !>  Test names to be comma-delimited
   !>  Doesn't care about whitespaces
   !
   subroutine set_unit_tests(this, unit_tests_str, mpi_env)
      !> unit tests
      class(unit_tests_type), intent(inout) :: this
      !> String containing all unit tests to run
      character(len=*), intent(in) :: unit_tests_str
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env

      !> Delimiter for unit test names in unit_tests_str
      character(len=1), parameter :: delimiter = ","
      !> Indices
      integer :: i, start

      start = 1
      do i = 1, len(unit_tests_str)
         if (unit_tests_str(i:i) == delimiter) then
            call set_unit_test(this, unit_tests_str(start:i - 1), mpi_env)
            start = i + 1
         end if
      end do

      ! Last test name does not end in a comma
      ! hence is missed in the do loop
      call set_unit_test(this, unit_tests_str(start:), mpi_env)

   end subroutine set_unit_tests

   !> Set the logical in `run` that corresponds to the string `test_name`.
   subroutine set_unit_test(run, test_name, mpi_env)
      !> Unit tests to run
      type(unit_tests_type), intent(inout) :: run
      !> Unit test name
      character(len=*), intent(in) :: test_name
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env

      select case (test_name)
      case ('all')
         run%all = .true.
         ! src_
      case ('advanced')
         run%advanced = .true.
      case ('eigensystem')
         run%eigensystem = .true.
      case ('fermisurfdx')
         run%fermisurfdx = .true.
      case ('gw')
         run%gw = .true.
      case ('hybrids')
         run%hybrids = .true.
      case ('LDAU')
         run%LDAU = .true.
      case ('mixing')
         run%mixing = .true.
      case ('mpi')
         run%mpi = .true.
      case ('optics')
         run%optics = .true.
      case ('phonon')
         run%phonon = .true.
      case ('raman')
         run%raman = .true.
      case ('rdmft')
         run%rdmft = .true.
      case ('stm')
         run%stm = .true.
      case ('sym')
         run%sym = .true.
      case ('vdw')
         run%vdw = .true.
      case ('wannier')
         run%wannier = .true.
      case ('xc')
         run%xc = .true.
      case ('xs')
         run%xs = .true.


         ! Modern directory conventions
      case ('xstring')
         run%xstring = .true.
      case ('lapack')
         run%lapack = .true.
      case ('math')
         run%math = .true.
      case('simplified_input')
         run%simplified_input = .true.
      case ('structure')
         run%structure = .true.
      case ('testframework')
         run%testframework = .true.
      case ('file_io')
         run%file_io = .true.
      case ('matrix_elements')
         run%matrix_elements = .true.
      case ('xgrid')
         run%xgrid = .true.
      case ('xhdf5')
         run%xhdf5 = .true.

      case default
         call terminate_mpi_env(mpi_env, 'Unrecognised unit test name: '&
            & //trim(test_name)//'. Please add to set_unit_test and unit_test_names_type')
      end select

   end subroutine set_unit_test

end module unit_tests
