!> Module for unit tests for for the functions in grid_utils.
module modmpi_test
   use precision, only: dp
   use constants, only: zone
   use modmpi, only: mpiinfo, find_2d_grid, distribute_loop
   use unit_test_framework, only: unit_test_type
   use math_utils, only: all_close

   implicit none

   private
   public :: modmpi_test_driver

contains

   !> Run tests for modmpi
   subroutine modmpi_test_driver(mpiglobal, kill_on_failure)
      !> mpi environment
      type(mpiinfo), intent(in) :: mpiglobal
      !> Kill the program upon failure of an assertion
      logical, intent(in), optional :: kill_on_failure

      !> Test report object
      type(unit_test_type) :: test_report
      !> Number of assertions
      integer, parameter :: n_assertions = 100

      call test_report%init(n_assertions, mpiglobal)

      ! Run unit tests
      call test_find_2d_grid(test_report)
      call test_distribute_loop(test_report)

      if (present(kill_on_failure)) then
         call test_report%report('modmpi', kill_on_failure)
      else
         call test_report%report('modmpi')
      end if

      call test_report%finalise()
   end subroutine modmpi_test_driver


   !> Test generation of 2D grid
   subroutine test_find_2d_grid(test_report)
      !> Test report object
      type(unit_test_type), intent(inout) :: test_report

      integer :: n_groups
      integer :: n_processes
      integer :: rows_per_group
      integer :: cols_per_group

      n_groups = 3
      n_processes = 30
      call find_2d_grid(n_processes, n_groups, rows_per_group, cols_per_group)
      call test_report%assert( rows_per_group == 5,&
                      & 'Test MPI Grid general. Expected number of rows: 5')
      call test_report%assert( cols_per_group == 2,&
                      & 'Test MPI Grid general. Expected number of cols: 2')


      n_groups = 30
      n_processes = 30
      call find_2d_grid(n_processes, n_groups, rows_per_group, cols_per_group)
      call test_report%assert( rows_per_group == 1, &
                        & 'Test MPI Grid n_groups = n_processes. Expected number of rows: 1')
      call test_report%assert( cols_per_group == 1, &
                        & 'Test MPI Grid n_groups = n_processes. Expected number of cols: 1')

      n_groups = 1
      n_processes = 30
      call find_2d_grid(n_processes, n_groups, rows_per_group, cols_per_group)
      call test_report%assert( rows_per_group == 6, &
                        & 'Test MPI Grid n_groups = 1. Expected number of rows: 6')
      call test_report%assert( cols_per_group == 5, &
                        & 'Test MPI Grid n_groups = 1. Expected number of cols: 5')

      n_groups = 3
      n_processes = 27
      call find_2d_grid(n_processes, n_groups, rows_per_group, cols_per_group)
      call test_report%assert( rows_per_group == 3, &
                        & 'Test MPI Grid sqrt(ranks_per_group) is integer. Expected number of rows: 3')
      call test_report%assert( cols_per_group == 3, &
                        & 'Test MPI Grid sqrt(ranks_per_group) is integer. Expected number of cols: 3')                 


      n_groups = 0
      n_processes = 27
      call find_2d_grid(n_processes, n_groups, rows_per_group, cols_per_group)
      call test_report%assert( rows_per_group == 1, &
                        & 'Test MPI Grid n_groups < 1. Expected number of rows: 1')
      call test_report%assert( cols_per_group == 1, &
                        & 'Test MPI Grid n_groups < 1. Expected number of cols: 1')                        
      
   end subroutine test_find_2d_grid


   subroutine test_distribute_loop(test_report)
      !> Test report object
      type(unit_test_type), intent(inout) :: test_report
      !> Mock mpi environment, such that tests can run in serial
      type(mpiinfo) :: mpi_mock_env
      !> Local test variables
      integer :: n_elements, first, last, i

      ! n_elements evenly divisible by the number of processes
      n_elements = 10
      mpi_mock_env%procs = 2
      ! [1, 2, 3, 4, 5 | 6, 7, 8, 9, 10]

      mpi_mock_env%rank = 0
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 1, message = "Expect first = 0 for rank 0")
      call test_report%assert(last == 5, message = "Expect last = 5 for rank 0")

      mpi_mock_env%rank = 1
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 6, message = "Expect first = 6 for rank 1")
      call test_report%assert(last == 10, message = "Expect last = 10 for rank 1")


      ! n_elements not evenly divisible by the number of processes
      ! Excess elements distributed amongst the lowest ranks
      n_elements = 10
      mpi_mock_env%procs = 3

      ! Expect distribution [1, 2, 3, 4| 5, 6, 7| 8, 9, 10]
      mpi_mock_env%rank = 0
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 1, message = "Expect first = 1 for rank 0")
      call test_report%assert(last == 4, message = "Expect last = 4 for rank 0")

      mpi_mock_env%rank = 1
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 5, message = "Expect first = 5 for rank 1")
      call test_report%assert(last == 7, message = "Expect last = 7 for rank 1")

      mpi_mock_env%rank = 2
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 8, message = "Expect first = 8 for rank 2")
      call test_report%assert(last == 10, message = "Expect last = 10 for rank 2")


      !rank receiving (first, last) = [101, 110] from 200 elements 
      n_elements = 200
      mpi_mock_env%procs = 20

      mpi_mock_env%rank = 10
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 101, &
            message = "n_elements divisible by nproc. Expect first = 101 for rank 10")
      call test_report%assert(last == 110, &
            message = "n_elements divisible by nproc. Expect first = 110 for rank 10")
                     

      ! Test for number of processes > n_elements.   
      mpi_mock_env%procs = 10 
      n_elements = 9 
      call test_report%assert(mpi_mock_env%procs > n_elements, message = "Expect n_processes > n_elements")

      ! Processes 0 to (last_rank - 1) get one element each i = [ 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9]
      ! where last_rank = mpi_mock_env%procs -1 
      do i = 1, n_elements
         mpi_mock_env%rank = i - 1
         call distribute_loop(mpi_mock_env, n_elements, first, last)
         call test_report%assert(first == i, message = "Processes 0-8 get one element") 
         call test_report%assert(first == last, message = "first == last")
      enddo

      ! Process 9 should be idle. Returned limits causes `do first, last` to not get evaluated
      mpi_mock_env%rank = 9
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 0, message = "Process 9 should be idle. Expect first = 0") 
      call test_report%assert(last == -1, message = "Process 9 should be idle. Expect first = -1")

      
      ! n_elements = 0
      n_elements = 0
      mpi_mock_env%procs = 4
      mpi_mock_env%rank = 0
      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 0, message = "If n_elements == 0, expect first = 0") 
      call test_report%assert(last == -1, message = "If n_elements == 0, expect first = -1")


      ! Serial behaviour
      mpi_mock_env%procs = 1 
      mpi_mock_env%rank = 0
      n_elements = 9

      call distribute_loop(mpi_mock_env, n_elements, first, last)
      call test_report%assert(first == 1, message = "No distribution. Expect first = 1") 
      call test_report%assert(last ==  9, message = "no distribution. Expect first = 9")

   end subroutine test_distribute_loop

   
end module modmpi_test
