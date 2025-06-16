module task_sigmac
  use asserts, only: assert
  use constants, only: zzero, real_zero
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use math_utils, only: all_zero
  use modinput, only: input, gw_type
  use modgw, only: freq, ibgw, nbgw
  use modmpi, only: mpiglobal, terminate_if_false, distribute_loop
  use mod_coulomb_potential, only: read_barcev_vmat_from_file, calculate_sqrt_bare_coulomb, delete_coulomb_potential, barc
  use mod_dielectric_function, only: read_inverse_epsilon_from_file, epsilon
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_product_basis, only: read_sgi_from_file
  use mod_selfenergy, only: selfec, write_selfec_single_kpoint, & 
    generate_frequency_grid_for_correlation_self_energy
  use precision, only: i32, dp

  implicit none

  private

  public :: execute_task_sigmac

  integer(i32), parameter :: max_length = 30

  character(len=*), parameter :: task_name = "sigmac"

  !> Interface to the parameters defined in the input file
  type task_sigmac_parameters
    private
    type(kpoints_sets) :: k_points
    integer(i32) :: n_omega
    real(dp) :: eigenvalue_cutoff_Coulomb_matrix
    character(len=max_length) :: output_format
  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_kpt )
  class(task_sigmac_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_kpt

  call this%sanity_checks( gw_inp )
  call this%k_points%parse_input( gw_inp%taskGroup%sigmac%kpointsarray, n_kpt )
  this%n_omega = gw_inp%freqgrid%nomeg
  this%output_format = trim( gw_inp%taskGroup%outputFormat )
  this%eigenvalue_cutoff_Coulomb_matrix = gw_inp%barecoul%barcevtol
end subroutine

subroutine sanity_checks( this, gw_inp )
  class(task_sigmac_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%sigmac), &
    'Element sigmac must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%freqgrid), &
    'Element freqgrid must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%barecoul), &
    'Element barecoul must be present when executing '//'"'//task_name//'"' )
end subroutine


subroutine sanity_check_epsilon_and_barc()
  call terminate_if_false( size( barc, 2 ) == size( epsilon, 1 ), &
    'Coulomb matrix and the inverse of epsilon have incompatible sizes' )
end subroutine

subroutine sanity_check_frequencies_of_epsilon()
  call terminate_if_false( freq%nomeg == size( epsilon, 3 ), &
    'Number of frequency points read for the inverse of epsilon incompatible with current calculation' )
end subroutine


!> Obtain the correlation part of the self-energy
subroutine execute_task_sigmac( n_kpoints_max, qpoints, file_format )
  !> Maximum number of irreducible k-points
  integer(i32), intent(in) :: n_kpoints_max
  !> List of q-points in cartesian coordinates
  real(dp), intent(in) :: qpoints(:, :)
  !> Format of input/output files
  character(len=*), intent(in) :: file_format

  integer(i32) :: ik, i, i_start, i_end, omega_i, omega_f
  integer(i32) :: iq, iq_start, iq_end, n_qpoints
  integer(i32), parameter :: maxlen=80
  character(len=maxlen) :: string
  real(dp), parameter :: tolerance_zero_vector = 1.e-6_dp
  real(dp) :: eigenvalue_cutoff
  type(task_sigmac_parameters) :: input_parameters
  type(mpiinfo) :: mpi_environment_kpoints
  
  call assert( size( qpoints, 1 ) == 3, 'qpoints must have size 3 along 1st dimension' )
  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  call input_parameters%parse_input( input%gw, n_kpoints_max )
  call input_parameters%k_points%obtain_list_of_indexes()
  call mpi_environment_kpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_kpoints, size(input_parameters%k_points%list_of_indexes), &
    i_start, i_end )

  omega_i = 1
  omega_f = input_parameters%n_omega
  n_qpoints = size( qpoints, 2 )
  iq_start = 1
  iq_end = n_qpoints
  eigenvalue_cutoff = max( real_zero, input_parameters%eigenvalue_cutoff_Coulomb_matrix )
  call generate_frequency_grid_for_correlation_self_energy( input%gw )
  do i = i_start, i_end
    ik = input_parameters%k_points%list_of_indexes(i)
    if( mpiglobal%rank == 0 ) then
      write( string, * ) ik
      string = '('//task_name//'): rank 0 -> calculating k-point with ik = ' // trim( string )
      call write_to_gwinfo( string )
    end if
    if( allocated( selfec )) deallocate( selfec )
    allocate( selfec(ibgw:nbgw, omega_i:omega_f, ik:ik), source=zzero )
    do iq = iq_start, iq_end
      call read_sgi_from_file( iq, file_format )
      call calcmpwipw( iq )
      call read_barcev_vmat_from_file( iq, file_format )
      Gamma = gammapoint( qpoints(:,iq), tol=tolerance_zero_vector )
      call calculate_sqrt_bare_coulomb( iq, eigenvalue_cutoff, Gamma )
      call read_inverse_epsilon_from_file( iq, Gamma, file_format )
      call sanity_check_frequencies_of_epsilon()
      call sanity_check_epsilon_and_barc()
      call calcselfc(iq, ik, ik)
    end do
    call write_selfec_single_kpoint( ik, file_format )
  end do
  call delete_coulomb_potential

end subroutine

end module