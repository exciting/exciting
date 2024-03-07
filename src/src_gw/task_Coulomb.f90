!> This module contains classes and subroutines needed to execute the 
!> task `Coulomb`, used as an element of `taskGroup` in `gw`
module task_Coulomb
  use exciting_mpi, only: mpiinfo
  use gw_io, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use modgw, only: fgw, kqset, Gqset
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, mpiglobal, distribute_loop
  use mod_coulomb_potential, only: calculate_sqrt_bare_coulomb, write_barc_to_file
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_product_basis, only: matsiz, locmatsiz, write_sgi_to_file
  use precision, only: i32, dp
  implicit none

  private

  character(len=*), parameter :: task_name = "Coulomb"

  public :: execute_task_Coulomb

  !> Interface to the parameters defined in the input file
  type task_Coulomb_parameters
    !> indexes of the q-points for which the bare Coulomb matrix must be calculated
    type(kpoints_sets) :: q_points
    !> type of cutoff used to obtain the bare Coulomb matrix
    character(len=20) :: Coulomb_cutoff_type
    !> threshold to include/eliminate eigenvectors when building the bare Coulomb matrix
    real(dp) :: Coulomb_eigenvalue_tol
  contains
    procedure :: parse_input, sanity_checks
  end type

contains 

!> Subroutine to be invoked when task `Coulomb` must be executed
subroutine execute_task_Coulomb( n_qpoints_max, is_output_format_binary )
  !> Maximum number of q-points for the system that is being calculated
  integer(i32), intent(in) :: n_qpoints_max
  !> When true, print output files in binary format
  logical, intent(in) :: is_output_format_binary

  integer(i32) :: iq, i, i_start, i_end
  integer(i32), parameter :: maxlen=60
  character(len=maxlen) :: string
  type(mpiinfo) :: mpi_environment_qpoints
  type(task_Coulomb_parameters) :: input_parameters

  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )

  call input_parameters%parse_input( input%gw, n_qpoints_max )
  call input_parameters%q_points%obtain_list_of_indexes()
  call mpi_environment_qpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_qpoints, size(input_parameters%q_points%list_of_indexes), &
    i_start, i_end )

  do i = i_start, i_end
    iq = input_parameters%q_points%list_of_indexes(i)
    ! Print a message to GW_INFO about the current q-point
    if( mpiglobal%rank == 0 ) then
      write( string, * ) '('//task_name//'): q-point cycle, iq = ', iq
      call write_to_gwinfo( string )
    end if
    Gamma = gammapoint(kqset%vqc(:,iq))
    matsiz = locmatsiz+Gqset%ngk(1, iq )
    ! Obtain an orthonormal set of interstitial plane waves
    call diagsgi( iq )
    call write_sgi_to_file( iq, is_output_format_binary )
    ! Calculates the matrix elements between PW's and orthonormalized IPW's
    call calcmpwipw( iq )
    call calculate_sqrt_bare_coulomb( iq, input_parameters%Coulomb_eigenvalue_tol )
    call write_barc_to_file( iq, is_output_format_binary )
  end do
end subroutine

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_qpt )
  class(task_Coulomb_parameters), intent(inout) :: this
  !> type with the variables given in the input file (only gw element)
  type(gw_type):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_qpt

  integer(i32) :: n_sets

  call this%sanity_checks( gw_inp )
  call this%q_points%parse_input( gw_inp%taskGroup%Coulomb%qpointsarray, n_qpt )
  this%Coulomb_cutoff_type = trim( gw_inp%barecoul%cutofftype )
  this%Coulomb_eigenvalue_tol = gw_inp%barecoul%barcevtol
end subroutine

!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  class(task_Coulomb_parameters), intent(in) :: this
  !> type with the variables given in the input file (inside the `gw` element)
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%Coulomb), &
    'Element Coulomb must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%barecoul), &
    'Element barecoul must be present when executing '//'"'//task_name//'"' )

end subroutine

end module