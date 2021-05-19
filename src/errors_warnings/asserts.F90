!> Module implementing assertions, emulating the C library macro `void assert`
module asserts
    use iso_fortran_env, only: error_unit

    use trace, only: trace_back
    implicit none
    private

    !> Error code for Exciting to return to the invoking environment
    Integer, Parameter :: error_code_logical = -101

    !>  A generic interface for assertions
    Interface assert
        Module Procedure assert_true
    End Interface assert

    ! Pubically-visible interface 
    public :: assert

contains

    !> Terminate following a failed assertion. 
    !>
    !> Cannot use the terminate subroutine in the \p modmpi 
    !> because we (eventually) want the MPI routines to use asserts,
    !> which would lead to a circular dependency. 
    !>
    !> @todo Note, this kills ALL mpi processes. One may not want to do this 
    !> if the communicator has been split
    subroutine terminate()
#ifdef MPI
        use mpi, only: MPI_COMM_WORLD, MPI_ABORT
        !> Dummy MPI error integer
        integer :: ierr 
        Call MPI_ABORT(MPI_COMM_WORLD, error_code_logical, ierr) 
#else
        stop error_code_logical
#endif        
    end subroutine terminate

    
    !> Assert if a logical condition is true. 
    !>
    !> If an assertion, passed as `logical_condition`, is false,
    !> assert_true first provides a starck trace of exciting's call 
    !> stack, then terminates.  
    !>
    !> One notes that due to a lack of a intrinsic backtrace function
    !> this behaviour will differ between GCC and Intel. Indeed, Intel
    !> will kill execution following the stack trace, whereas GCC 
    !> allows it to continue. 
    !>
    !> If not compiled in USE_ASSERT (debug) mode, the compiler should
    !> optimise by removing the routine, which will be empty.
    !>
    subroutine assert_true(logical_condition, message)
        !> Condition to test
        logical, intent(in) :: logical_condition
        !> Optional message 
        character(len=*), intent(in), optional :: message
#ifdef USE_ASSERT
        if (.not. logical_condition) then
            if (present(message)) then
                write (error_unit, '(/,1x,a)') trim(adjustl(message))
            endif
            call trace_back()
            call terminate()
        end if
#endif
    end subroutine assert_true

end module asserts