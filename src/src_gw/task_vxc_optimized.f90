!> Module designed for the task vxc_optimized,
!> i.e. the one that prints the optimized
!> XC potential for a QSGW run
module task_optimized_vxc
  use asserts, only: assert
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo_boxmessage
  use modmpi, only: terminate_if_false, mpiglobal, distribute_loop, barrier
  use mod_kqpts, only: kpoints_sets
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32, str_128
  use constants, only: zzero
  use mod_qsgw, only: compute_optimized_vxc, write_optimized_vxc_to_a_file
  use mod_gw_degeneracies, only: ibgw_including_degeneracy, nbgw_including_degeneracy
  use mod_band_to_lapw_transform, only: transform_from_band_representation_to_lapwlo
  use mod_selfenergy, only: read_selfec_from_files, read_selfex_from_files, generate_frequency_grid_for_correlation_self_energy
  use gw_info, only: write_to_gwinfo
  use modgw,   only: kset
  use mod_selfconsistent_gw, only: gw_first_iteration

  implicit none 

  private

  character(len=*), parameter :: task_name = "optimizedVxc"

  !> Interface to the parameters defined in the input file
  type task_optimized_vxc_parameters
    private
    type(kpoints_sets) :: k_points
    logical :: offdiagonal
  contains
    procedure :: parse_input, sanity_checks
  end type

  public :: execute_task_optimized_vxc

contains
!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_kpt )
  class(task_optimized_vxc_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in):: gw_inp
  !> maximum number of k-points
  integer(i32), intent(in) :: n_kpt
  

  call this%sanity_checks( gw_inp )
  call this%k_points%parse_input( gw_inp%taskGroup%optimizedVxc%kpointsarray, n_kpt )
  this%offdiagonal = gw_inp%taskGroup%optimizedVxc%offdiagonal

end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  class(task_optimized_vxc_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%optimizedVxc), &
    'Element optimized_vxc must be present when executing '//'"'//task_name//'"' )

end subroutine


!> Execute a task optimized_vxc calculation
subroutine execute_task_optimized_vxc( n_kpoints_max, file_format )
  integer(i32), intent(in) :: n_kpoints_max
  character(len=*), intent(in) :: file_format
  integer(i32) :: ik
  integer(i32) :: i, i_start, i_end
  type(task_optimized_vxc_parameters) :: input_parameters
  character(len=str_128) :: string
  complex(dp), allocatable :: vxc_opt(:,:,:), vxc_opt_lapwlo(:,:)
  type(mpiinfo) :: mpi_environment_kpoints
  
  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  
  call input_parameters%parse_input( input%gw, n_kpoints_max )
  call input_parameters%k_points%obtain_list_of_indexes()
  call mpi_environment_kpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_kpoints, size(input_parameters%k_points%list_of_indexes), &
    i_start, i_end )

  ! Handle the processes without jobs
  if (i_start == 0 .and. i_end == -1) return

  ! Iterate over the k-points assigned to this MPI rank
  do i = i_start, i_end
    ! Read the diagonal elements of the self-energy
    call read_selfec_from_files( input_parameters%k_points%list_of_indexes(i:i), file_format )
    call read_selfex_from_files( input_parameters%k_points%list_of_indexes(i:i), file_format )
    if ( trim(input%gw%selfenergy%method) == "ac" ) then
      ! Analytical continuation of the correlation self-energy from the complex to the real frequency axis
      call generate_frequency_grid_for_correlation_self_energy( input%gw )
      call calcselfc_ac()
    end if
    ik = input_parameters%k_points%list_of_indexes(i)
    if( mpiglobal%rank == 0 ) then
      write( string, * ) '('//task_name//'): k-point cycle, ik (irreducible) = ', ik
      call write_to_gwinfo( string )
    end if
    allocate(vxc_opt(ibgw_including_degeneracy:nbgw_including_degeneracy, &
                   ibgw_including_degeneracy:nbgw_including_degeneracy, &
                   ik:ik), source=zzero)
    ! Compute the optimized potential in band representation
    ! If selected use the offdiagonal elements of the self-energy 
    ! in addition to the diagonal ones to compute the optimized QSGW
    ! potential.
    call compute_optimized_vxc(ik, input_parameters%offdiagonal, vxc_opt, file_format)
    ! Transform to LAPW+LO representation
    call transform_from_band_representation_to_lapwlo(kset%ikp2ik(ik), ik, &
                                                      gw_first_iteration(), &
                                                      vxc_opt(:,:,ik), &
                                                      vxc_opt_lapwlo, file_format)
    ! ! Save the potential to a file
    call write_optimized_vxc_to_a_file(vxc_opt_lapwlo, ik, file_format)
    deallocate(vxc_opt)
  end do

end subroutine execute_task_optimized_vxc

end module task_optimized_vxc
