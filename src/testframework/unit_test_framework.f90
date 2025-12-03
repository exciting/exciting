!> Module providing a lightweight unit testing object.
!> The assertion routine signature has deliberately been written
!> to copy that of the [Zofu unit testing framework](https://github.com/acroucher/zofu),
!> which should ultimately supersede this module.
module unit_test_framework
   use precision
   use modmpi, only: mpiinfo, terminate_mpi_env, mpiglobal

   implicit none
   private
   public :: unit_test_type

   !> Message size
   integer(i32), parameter :: message_size = 512

   !> Unit Test object
   !> Collect the results and error messages from unit tests
   type unit_test_type
      !> Collect an optional message for each assertion
      character(len=message_size), allocatable, private :: messages(:)
      !> Stores the indices of failed assertions
      integer, allocatable, private :: failures(:)
      !> Internal state index
      integer, private :: i
      !> MPI environment
      type(mpiinfo), private :: mpi_env
   contains
      procedure :: init
      procedure :: assert
      procedure :: report
      procedure :: finalise
   end type unit_test_type

contains

   !> Initialise the unit test object
   !>
   !> This routine copies the MPI environment instance
   !> such that the methods have access to it without having
   !> to modify the API of assert (and therefore break compatibility
   !> with zofu). this%mpi_env will only be used to terminate
   !> the code.
   subroutine init(this, mpi_env)
      !> Test report object
      class(unit_test_type), intent(inout) :: this
      !> MPI environment
      type(mpiinfo), intent(in), optional :: mpi_env

      allocate (this%messages(0))
      allocate (this%failures(0))

      this%messages = ""
      this%failures = 0
      this%i = 1
      if(present(mpi_env)) then
         this%mpi_env = mpi_env
      else
         this%mpi_env = mpiglobal
      end if
   end subroutine init

   !> Log failed assertions
   subroutine assert(this, condition, message)
      !> Test report object
      class(unit_test_type), intent(inout) :: this
      !> Assertion status
      logical, intent(in) :: condition
      !> Optional assertion error message
      character(len=*), intent(in) :: message


      this%messages = [this%messages, trim(adjustl(message))]

      if (.not. condition) then
         this%failures = [this%failures, this%i]
      end if

      this%i = this%i + 1
   end subroutine assert

   !> Report results
   !>
   !> @todo Note, the failure reporting is going to get messy for multiple processes.
   !> This should be addressed as soon as parallel tests are written.
   subroutine report(this, name, kill_on_failure)

      !> Unit test object
      class(unit_test_type), intent(in) :: this
      !> Name of the test
      character(len=*), intent(in) :: name
      !> Kill the program before the test driver finishes, if an assertion fails
      logical, intent(in), optional :: kill_on_failure

      !> Copy of the mpi environment
      type(mpiinfo) :: mpi_env
      integer :: i, n_assertions, n_failures
      character(len=10) :: process_id

      n_assertions = size(this%messages)
      n_failures = size(this%failures)

      write(process_id, '(I0)') this%mpi_env%rank
      if (this%mpi_env%is_root) then
         write (*, *) 'Report test: '//trim(adjustl(name))
      endif

      ! I3 assumes less than 1000 assertions
      write (*, '(X, A, I0, A, I0)') 'Assertions passed on process '//trim(process_id)//':', &
         n_assertions - n_failures, ' out of ', n_assertions

      if (n_failures > 0) then
         write (*, *) 'Failing tests on process '//trim(process_id)//':'
         do i = 1, n_failures
            write (*, *) 'Error message: ', trim(adjustl(this%messages(this%failures(i))))
         end do

         if (present(kill_on_failure)) then
            if (kill_on_failure) then
               ! Avoids making 'this' intent(inout)
               mpi_env = this%mpi_env
               call terminate_mpi_env(mpi_env, "unit test failed")
            end if
         end if
      end if
      write(*,*) ' '

   end subroutine report

   !> Clean up unit test object by deallocating
   !> any allocatable data.
   subroutine finalise(this)
      !> Test report object
      class(unit_test_type), intent(inout) :: this

      this%i = 0
      if(allocated(this%messages)) deallocate (this%messages)
      if(allocated(this%failures)) deallocate (this%failures)
   end subroutine finalise

end module unit_test_framework
