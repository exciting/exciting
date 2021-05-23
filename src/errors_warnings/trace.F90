!> Provide backtrace functionality. 
module trace
#ifdef __INTEL_COMPILER
    use ifcore, only: tracebackqq       
#endif 
    implicit none
    private 

    public :: trace_back

contains

    !> Performs a backtrace, ending termination. 
    !>
    !> There is no intrinsic backtrace function, so one is required 
    !> to implement vendour-specific routines. 
    !>
    !> Behaviour unfortunately varies between the GCC and Intel routines.
    !> GCC backtrace allows execution to continue, whereas Intel
    !> tracebackqq terminates execution. As such, it is recommended to
    !> use this routine in conjunction with MPI_ABORT.  
    !>
    !> The preprocessor variables used for GCC and Intel are macros
    !> defined by the compilers, not the developers of exciting.
    !>
    !> For example, one can check  check macros available to GCC by
    !> typing: 
    !>   touch dummy.f90
    !>   gfortran -cpp -E -dM dummy.f90
    !> 
    !> References given by the [fortran wiki](http://fortranwiki.org/fortran/show/Predefined+preprocessor+macros)
    !> 
    !> If exciting has not been compiled with GCC or Intel, do nothing.
    !> If exciting has been compiled in production mode, do nothing as
    !> the debug symbols have not been used. 
    subroutine trace_back(error_message)
        character(len=*), optional, intent(in) :: error_message
#ifdef USE_ASSERT
        character(len=256) :: message 

        if (present(error_message)) then
            message = trim(error_message)
        else 
            message = "Trace back:"
        endif

#ifdef __GFORTRAN__
        write(*,*) trim(message)
        call backtrace()    
#elif __INTEL_COMPILER
        call tracebackqq(string=message) 
#else 
        write(*,*) 'trace_back not overloaded for this compiler' 
        continue
#endif
#endif 

    end subroutine trace_back

end module
