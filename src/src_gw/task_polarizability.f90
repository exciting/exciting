!> This module contains classes and subroutines needed to compute the polarizability, 
!> used as an element of `taskGroup` in `gw`
module task_polarizability

  use constants, only: real_zero
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use math_utils, only: all_zero
  use modgw, only: kqset, kset, Gqset
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, mpiglobal, distribute_loop
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: gammapoint
  use mod_product_basis, only: locmatsiz, matsiz, read_sgi_from_file
  use mod_polarizability, only: init_polarizability, compute_polarizability_at_q, write_polarizability_to_file, delete_polarizability
  use precision, only: dp, i32
  use to_char_conversion, only: to_char

  implicit none

  private

  public :: execute_task_polarizability

  character(len=*), parameter :: task_name = "polarizability"

  !> Interface to the parameters defined in the input file
  type task_polarizability_parameters
    private
    type(kpoints_sets) :: q_points
    logical :: usingIrreducibleWedge
    integer(i32) :: n_omega
  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_qpt )
  class(task_polarizability_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_qpt
  
  call this%sanity_checks( gw_inp )
  call this%q_points%parse_input( gw_inp%taskGroup%polarizability%qpointsarray, n_qpt )
  this%usingIrreducibleWedge = gw_inp%taskGroup%polarizability%usingIrreducibleWedge
  this%n_omega = gw_inp%freqgrid%nomeg
end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  class(task_polarizability_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%polarizability), &
    'Element polarizability must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%freqgrid), &
    'Element freqgrid must be present when executing '//'"'//task_name//'"' )
end subroutine


!> Subroutine to be invoked when task `polarizability` must be executed
subroutine execute_task_polarizability( n_qpoints_max, file_format )
  integer(i32), intent(in) :: n_qpoints_max
  character(len=*), intent(in) :: file_format

  integer(i32) :: iq, iq_reducible, iq_output
  integer(i32) :: i, i_start, i_end
  integer(i32) :: omega_i, omega_f
  type(task_polarizability_parameters) :: input_parameters
  type(mpiinfo) :: mpi_environment_qpoints
  logical :: is_Gamma

  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  call input_parameters%parse_input( input%gw, n_qpoints_max )
  call input_parameters%q_points%obtain_list_of_indexes()
  call mpi_environment_qpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_qpoints, size(input_parameters%q_points%list_of_indexes), &
    i_start, i_end )
  
  omega_i = 1
  omega_f = input_parameters%n_omega
  
  ! Attention: calcpmatgw makes use of MPI parallelization and calls a mpi_barrier
  if( isGammaInList( kqset%vqc(:,input_parameters%q_points%list_of_indexes) ) ) call calcpmatgw
  do i = i_start, i_end
    iq = input_parameters%q_points%list_of_indexes(i) ! This refers always to the list either full or irreducible
    if (input_parameters%usingIrreducibleWedge) then
      if( mpiglobal%rank == 0) call write_to_gwinfo( '('//task_name//'): q-point cycle, iq (irreducible) = ' // to_char(iq) )
      iq_reducible = kset%ikp2ik(iq) ! iq is the index of the irreducible q-point; iq_reducible is the index in the reducible q-point list
      iq_output = iq ! We are outputing the files with the irreducible wedge numbering
    else
      if( mpiglobal%rank == 0) call write_to_gwinfo( '('//task_name//'): q-point cycle, iq = ' // to_char(iq) )
      iq_reducible = iq ! iq is the index in the full BZ
      iq_output = iq_reducible ! We are outputing files with full BZ numbering
    end if

    ! We need to set matsiz to the appropiate value
    matsiz = locmatsiz + Gqset%ngk(1,iq_reducible)
    call read_sgi_from_file( iq_reducible, file_format )
    call calcmpwipw( iq_reducible )
    is_Gamma = gammapoint( kqset%vqc(:, iq_reducible), tol=1.e-6_dp )
    call init_polarizability( matsiz, omega_i, omega_f, is_Gamma)
    call compute_polarizability_at_q( iq_reducible, omega_i, omega_f, is_Gamma)
    call write_polarizability_to_file( iq_output, is_Gamma, file_format, input_parameters%usingIrreducibleWedge)
    call delete_polarizability()

  end do

end subroutine

!> Check if the gamma point is among the q-points
pure logical function isGammaInList( q_points )
  !> List of q-points (first dimension has the 3 cartesian coordinates)
  real(dp), intent(in) :: q_points(:, :)

  integer(i32) :: i

  isGammaInList = .false.
  do i = 1, size( q_points, 2 )
    isGammaInList = all_zero( q_points(:, i), tol=1.e-6_dp )
    if( isGammaInList ) exit
  end do
end function

end module task_polarizability
