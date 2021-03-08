!> Module providing a lightweight unit testing object. 
!> The assertion routine signature has deliberately been written
!> to copy that of the [Zofu unit testing framework](https://github.com/acroucher/zofu),
!> which should ultimately superceed this module. 
module unit_test_framework
  implicit none
  private
  public :: unit_test_type

  !> Unit Test object
  !> Collect the results and error messages from unit tests
  type unit_test_type
    !> Collect an optional message for each assertion 
    character(len=256), allocatable, private :: messages(:)
    !> Stores the indices of failed assertions 
    integer, allocatable, private :: failures(:)
    !> Number of failed assertions
    integer, private :: n_failures
    !> Internal state index
    integer, private :: i         
  contains
    procedure :: init
    procedure :: assert
    procedure :: report
    procedure :: finalise
  end type unit_test_type

contains
  
  !> Initialise the unit test object
  subroutine init(this, n_assertions)
    !> Test object   
    class(unit_test_type), intent(inout) :: this
    !> Total number of assertions
    integer, intent(in) :: n_assertions
    allocate(this%messages(n_assertions))
    ! Allocated with an upper bound
    allocate(this%failures(n_assertions))  
    this%messages = ""
    this%failures = 0
    this%n_failures = 0                  
    this%i = 1
  end subroutine init

  !> Log failed assertions
  subroutine assert(this, condition, message)
    !> Test object
    class(unit_test_type), intent(inout) :: this
    !> Assertion status
    logical, intent(in) :: condition
    !> Optional assertion error message 
    character(len=*), intent(in), optional :: message
    
    if (present(message)) then
       this%messages(this%i) = trim(adjustl(message))
    endif
       
    if (.not. condition) then
       this%n_failures = this%n_failures + 1
       this%failures(this%n_failures) = this%i
    end if
    
    this%i = this%i + 1
  end subroutine assert

  !> Report results
  subroutine report(this, name, kill_on_failure)
    use modmpi, only: terminate

    !> Unit test object
    class(unit_test_type), intent(in) :: this
    !> Name of the test
    character(len=*) :: name 
    !> Kill the program before the test driver finishes,
    !> if an assertion fails 
    logical, optional :: kill_on_failure 
    integer :: i, n_assertions, n_failures

    n_assertions = size(this%messages)
    n_failures = this%n_failures
    
    write(*,*) 'Report test: '//trim(adjustl(name))
    ! I3 assumes less than 1000 assertions
    write(*,'(X, A, I3, A, I3)') 'Assertions passed: ', n_assertions - n_failures, ' out of ', n_assertions

    if (n_failures > 0) then
      write(*,*) 'Failing tests:'
      do i = 1, n_failures
         write(*,*) 'Error message: ', trim(adjustl(this%messages(this%failures(i))))
      end do

      if (present(kill_on_failure)) then
         if (kill_on_failure) call terminate("feature test failed")
      endif
   end if

 end subroutine report

  !> Clean up unit test object
 subroutine finalise(this)
   !> Test object
   class(unit_test_type), intent(inout) :: this
   deallocate(this%messages)
   deallocate(this%failures)
 end subroutine finalise
    
end module unit_test_framework
