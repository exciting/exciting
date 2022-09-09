!> Parse arguments from the command line and store in an object. 
!>
!> In the future, we might consider adopting one of the libraries
!> listed on the [fortran wiki](http://fortranwiki.org/fortran/show/Command-line+arguments),
!> depending on the our needs.
module cmd_line_args
   use modmpi, only: mpiinfo, terminate_mpi_env
#ifdef MPI
   use mpi
#endif 
   implicit none
   private

   public :: cmd_line_args_type, get_unit_tests_string

   !> Arbitrary NULL value for solver threads.
   integer, parameter, public :: null_solver_threads = 1003629

   !> A type to hold all settings/options that can be passed via the
   !> command line to exciting
   type cmd_line_args_type
      !> Indicate whether to run unit tests
      logical :: run_unit_tests = .false.
      !> Kill the tests if an assertion fails
      logical :: kill_on_failure = .false.
      !> Number of threads used by the eigensolver, in the ground state.
      !> Option to store it as GW only gives fully consistent numbers
      !> between executions if the init SCF call is performed with a
      !> number of threads consistent with the ground state (one `assumes` 
      !> specifically for the eigensolver)
      integer :: ground_state_solver_threads
   contains
      !> Initialise by parsing arguments from the command line
      procedure :: parse => parse_command_line_args
   end type cmd_line_args_type

contains

   !> Parse arguments from the command line.
   !>
   !> Separate arguments are defined by the whitespace delimiter.
   !
   ! For reference, if one needs to convert a string to real:
   ! READ(arg1, *) num1
   ! where num1 is the real
   subroutine parse_command_line_args(this, mpi_env)
      !> Instance of command line options
      class(cmd_line_args_type), intent(inout) :: this
      !> Mpi environment
      type(mpiinfo), intent(inout) :: mpi_env

      !> Command line argument, defined by the whitespace delimiter
      character(len=100) :: arg
      integer :: i, j, end

      this%ground_state_solver_threads = null_solver_threads

      ! Parse command line arguments from root, only 
      if (mpi_env%is_root) then

         do i = 1, command_argument_count()
            call get_command_argument(i, arg)
   
            ! If setting=options is specified, apply select case to the setting only
            j = INDEX(arg, '=')
            if (j /= 0) then
               end = j - 1
            else
               end = len(trim(arg))
            end if
   
            select case(arg(:end))
            case('-run-unit-tests')
               this%run_unit_tests = .true.
            
            case('-kill-on-failure')
               this%kill_on_failure = .true.

            case('-gw-scf-threads')
               ! Get value following the matched option
               call get_command_argument(i + 1, arg)
               ! Not likely to have more than 999 threads
               read(arg,'(i3)') this%ground_state_solver_threads 
               cycle
            end select
   
         end do
   
      endif

      ! Broadcast object data to other processes
      ! If command-line arguments grow, one should bcast the whole object instance (this) instead 
#ifdef MPI
      call mpi_bcast(this%run_unit_tests, 1, MPI_LOGICAL, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      call mpi_bcast(this%kill_on_failure, 1, MPI_LOGICAL, mpi_env%root, mpi_env%comm, mpi_env%ierr)
      call mpi_bcast(this%ground_state_solver_threads, 1, MPI_LOGICAL, mpi_env%root, mpi_env%comm, mpi_env%ierr)
#endif 

   end subroutine parse_command_line_args


   !> Get the unit tests substring from the command line arguments.
   !>
   !> Specifically, extract "test, names" from the substring
   !>       "-run-unit-tests= test, names"
   !> if present in the command line arguments.
   !>
   !> Expects command line options to start with '-'
   !>
   !> If '-run-unit-tests' is present but no list of tests is specified, this will
   !> be treated as run all tests, and 'all' will be returned.
   !>
   !> If neither '-run-unit-tests' nor '-run-unit-tests=' are present,
   !> a character string of length 0 is returned.
   !
   function get_unit_tests_string(mpi_env) result(unit_tests_str)
      !> MPI environment
      type(mpiinfo), intent(inout) :: mpi_env

      !> Command line string argument, defined by whitespace delimiters
      character(len=100) :: arg
      !> All command line args as a single string, with no whitespace delimiters
      character(len=500) :: input
      !>  Integer returned by SCAN and INDEX if a substring is not present
      integer, parameter :: NULL = 0
      !> Indices in input character
      integer :: i, istart, iend, index_next_setting
      !> character string listing unit tests to run
      character(len=200) :: unit_tests_str

      ! Parse command line from root process, only
      if (mpi_env%is_root) then
         ! Parse command line arguments as a string with no delimiters
         call get_command_argument(1, input)
         do i = 2, command_argument_count()
            call get_command_argument(i, arg)
            input = trim(input)//trim(arg)
         end do
      endif

#ifdef MPI
      call mpi_bcast(input, len(input), MPI_CHARACTER, mpi_env%root, mpi_env%comm, mpi_env%ierr)
#endif

      ! Find where the '-run-unit-tests=' substring begins
      i = INDEX(input, '-run-unit-tests=')

      ! Limiting cases
      if (i == NULL) then
         if (INDEX(input, '-run-unit-tests') /= NULL) then
            unit_tests_str = 'all'
         else
            unit_tests_str = ''
         end if
         return
      end if

      ! Find where the list of unit tests starts
      istart = i + INDEX(input(i:), '=')

      ! Find where the list of unit tests ends
      index_next_setting = SCAN(input(istart:), '-')

      if (index_next_setting == NULL) then
         iend = len(trim(input))
      else
         iend = istart + index_next_setting - 2
      end if

      unit_tests_str = input(istart:iend)
   end function get_unit_tests_string

end module cmd_line_args
