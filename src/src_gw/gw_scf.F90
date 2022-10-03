!> Ensure the SCF routine is called from GW with a number of threads
!> equal to that used to compute the ground state density.
!>
!> ## Issue with Performing the Same Approach with openBLAS.
!>
!> OpenBLAS does not appear to provide functions that allow getting or setting of
!> active threads. One can read about it on these Github issues regarding,
!> [openBLAS](https://github.com/xianyi/OpenBLAS/issues/760) and [Julia](https://github.com/JuliaLang/julia/pull/21105)
!>
!> Specifically, it provides:
!> ```fortran
!>   ! Get the physical number of cpus on the machine
!>   call openblas_get_num_procs().
!>   ! Get total number of threads available to the system, rather than the number used by openBLAS
!>   call openblas_get_num_threads().
!>   ! Presumably set the number of available threads?
!>   call openblas_set_num_threads()
!> ```
!>
!> Alex tested the behaviour with the above routines, and was unable to achieve what is implemented with MKL.
!>
module gw_scf

#ifdef MKL
   ! Defined in  MKL's include directory
   use mkl_service, only: mkl_get_max_threads
#endif
   use modmpi, only: mpiinfo
   use cmd_line_args, only: null_solver_threads, cmd_line_args_type

   implicit none
   private
   public :: set_gs_solver_threads, thread_consistent_scf

contains

#ifdef MKL
   !> Set the threads used by the linear algebra lib, in a ground state calculation.
   !>
   !> Intel reference, as a [tiny URL](tinyurl.com/2p8f478s).
   function set_gs_solver_threads(mpiglobal) result(gs_solver_threads)
      !> Instance of the MPI env
      type(mpiinfo), intent(inout) :: mpiglobal
      !> Ground state solver threads
      integer :: gs_solver_threads
      !> Command-line arguments
      type(cmd_line_args_type) :: args

      ! Get number of threads used for the ground state eigensolver from cmd line
      call args%parse(mpiglobal)

      gs_solver_threads = args%ground_state_solver_threads

      ! If not provided as a cmd-line argument, one takes the available threads
      if (gs_solver_threads == null_solver_threads) then
         gs_solver_threads = mkl_get_max_threads()
      end if

   end function

   !> Run the non-self-consistent call to SCF using the number of
   !> threads for the calculation that generated
   !> the input density (in I assume, the eigensolver).
   subroutine thread_consistent_scf(gs_mkl_threads)
      !> N threads used in the solver of the ground state SCF calculation
      integer, intent(in) :: gs_mkl_threads
      !> MKL threads specified for GW calculation
      integer :: max_mkl_threads
      !> SCF verbosity
      integer, parameter :: verbosity = -2

      ! Linear albegra lib threads requested for GW calculation from MKL
      max_mkl_threads = mkl_get_max_threads()

      ! Compute density with the number of threads used for the ground state calculation
      ! Required for fully-consistent results between GW executions
      call mkl_set_num_threads(gs_mkl_threads)
      call mkl_set_dynamic(0)

      call scf_cycle(verbosity)

      ! Perform the GW linear algebra with all available threads
      call mkl_set_num_threads(max_mkl_threads)
      call mkl_set_dynamic(1)

   end subroutine

#else

   !> Dummy routine for setting the threads used by the linear algebra lib,
   !> in a ground state calculation.
   !>
   !> Inform the user that the command-line argument will be unused, if provided.
   function set_gs_solver_threads(mpiglobal) result(gs_solver_threads)
      !> Instance of the MPI env
      type(mpiinfo), intent(inout) :: mpiglobal
      !> Ground state solver threads
      integer :: gs_solver_threads
      !> Command-line arguments
      type(cmd_line_args_type) :: args

      ! Get number of threads used for the ground state eigensolver from cmd line
      call args%parse(mpiglobal)

      gs_solver_threads = args%ground_state_solver_threads

      if (args%ground_state_solver_threads /= null_solver_threads) then
         write(*, *) 'Threads for GW to use in SCF call have been provided &
                      but MKL is not being used. Expect this value to be unused.'
      end if
   end function

   !> Resort to undressed SCF call.
   !>
   !> User must manually ensure that the ground state calculation
   !> and the GW calculation are performed with the same number of
   !> linear algebra lib threads if one desires consistent QP energies.
   subroutine thread_consistent_scf(threads)
      !> Dummy arg to maintain consistent API
      integer, intent(in), optional :: threads
      !> SCF verbosity
      integer, parameter :: verbosity = -2
      ! Avoid warning forr unused arg
      if (present(threads)) continue
      call scf_cycle(verbosity)
   end subroutine
#endif

end module
