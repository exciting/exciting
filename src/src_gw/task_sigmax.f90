!> Module for the task sigmax
module task_sigmax
  use asserts, only: assert 
  use constants, only: zzero
  use exciting_mpi, only: mpiinfo
  use gw_io, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use modinput, only: input, gw_type
  use modmpi, only: mpiglobal, distribute_loop, terminate_if_false
  use mod_coulomb_potential, only: delete_coulomb_potential, read_barcev_vmat_from_file
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_product_basis, only: read_sgi_from_file
  use mod_selfenergy, only: selfex, write_selfex_single_kpoint
  use precision, only: dp, i32

  implicit none

  private

  public :: execute_task_sigmax

  character(len=*), parameter :: task_name = "sigmax"

  !> Interface to the parameters defined in the input file
  type task_sigmax_parameters
    private
    type(kpoints_sets) :: k_points
  contains
    procedure :: parse_input, sanity_checks
  end type

contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_kpt )
  class(task_sigmax_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_kpt

  call this%sanity_checks( gw_inp )
  call this%k_points%parse_input( gw_inp%taskGroup%sigmax%kpointsarray, n_kpt )

end subroutine


subroutine sanity_checks( this, gw_inp )
  class(task_sigmax_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%sigmax), &
    'Element sigmax must be present when executing '//'"'//task_name//'"' )

end subroutine


!> Execute a task sigmax calculation
subroutine execute_task_sigmax( first_band, last_band, n_kpoints_max, qpoints, file_format )
  !> Index of the first KS band for which sigmax is evaluated
  integer(i32), intent(in) :: first_band
  !> Index of the last KS band for which sigmax is evaluated
  integer(i32), intent(in) :: last_band
  !> Maximum number of irreducible k-points
  integer(i32), intent(in) :: n_kpoints_max
  !> List of q-points
  real(dp), intent(in) :: qpoints(:, :)
  !> Format of input/ouput files
  character(len=*), intent(in) :: file_format

  integer(i32) :: i, i_start, i_end, ik, iq, iq_start, iq_end, n_qpoints
  integer(i32), parameter :: maxlen=80
  character(len=maxlen) :: string
  real(dp), parameter :: tolerance_zero_vector = 1.e-6_dp
  type(task_sigmax_parameters) :: input_parameters
  type(mpiinfo) :: mpi_environment_kpoints

  call assert( size( qpoints, 1 ) == 3, 'qpoints must have size 3 along 1st dimension' )
  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  call input_parameters%parse_input( input%gw, n_kpoints_max )
  call input_parameters%k_points%obtain_list_of_indexes()
  call mpi_environment_kpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_kpoints, size(input_parameters%k_points%list_of_indexes), &
    i_start, i_end )
  n_qpoints = size( qpoints, 2 )
  iq_start = 1
  iq_end = n_qpoints
  
  do i = i_start, i_end
    ik = input_parameters%k_points%list_of_indexes(i)
    if( mpiglobal%rank == 0 ) then
      write( string, * ) ik
      string = '('//task_name//'): rank 0 -> calculating k-point with ik = ' // trim( adjustl( string ) )
      call write_to_gwinfo( string )
    end if
    if( allocated(selfex) ) deallocate( selfex )
    allocate( selfex(first_band:last_band, ik:ik ) , source=zzero )
    do iq = iq_start, iq_end
      call read_barcev_vmat_from_file( iq, file_format )
      call read_sgi_from_file( iq, file_format )
      call calcmpwipw( iq )
      Gamma = gammapoint( qpoints(:, iq), tol=tolerance_zero_vector )
      call calcselfx( iq, ik, ik )
    end do
    call write_selfex_single_kpoint( ik, file_format )
  end do
  call delete_coulomb_potential

end subroutine

end module