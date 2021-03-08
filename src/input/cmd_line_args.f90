!> Parse arguments from the command line and store in an object 
! In the future, we might consider adopting one of the libraries
! listed here, http://fortranwiki.org/fortran/show/Command-line+arguments,
! depending on the our needs. 

module cmd_line_args
  implicit none
  private
  public :: cmd_line_args_type

  !> A type to hold all arguments that can be passed via the
  !> command line to exciting
  type cmd_line_args_type
     !> Logical indicating whether to run unit tests 
     logical :: run_unit_tests = .false. 
   contains
     !> Initialise by parsing arguments from the command line
     procedure :: parse => parse_command_line_args
  end type cmd_line_args_type


contains

  !> Parse arguments from the command line
  !
  ! For reference, if one needs to convert a string to real:
  ! READ(arg1, *) num1
  ! where num1 is the real
  subroutine parse_command_line_args(this)
    !> Instance of command line options 
    class(cmd_line_args_type), intent(inout) :: this

    character(len=100) :: arg
    integer :: i 

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        select case (arg)
            case ('-run_unit_tests')
               this%run_unit_tests = .true. 

            case default
                write(*,*), 'unrecognised command-line option ignored: ', arg
        end select
    end do

  end subroutine parse_command_line_args


end module cmd_line_args
