!> This module contains classes and subroutines needed to execute the 
!> task `epsilon`, used as an element of `taskGroup` in `gw`
module task_epsilon
  use asserts, only: assert
  use constants, only: real_zero
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage
  use gw_io, only: write_to_file
  use math_utils, only: all_zero
  use modgw, only: kqset, kset
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, mpiglobal, distribute_loop
  use mod_coulomb_potential, only: delete_coulomb_potential, read_barcev_vmat_from_file, calculate_sqrt_bare_coulomb
  use mod_dielectric_function, only: write_epsilon_to_file, init_dielectric_function, delete_dielectric_function
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_product_basis, only: mbsiz, matsiz, read_sgi_from_file, mpwipw
  use mod_selfenergy, only: singc1, singc2
  use mod_polarizability, only: init_polarizability, from_polarizability_to_epsilon, read_polarizability_from_file, delete_polarizability
  use precision, only: dp, i32
  use to_char_conversion, only: to_char

  implicit none

  private

  public :: execute_task_epsilon

  integer(i32), parameter :: max_length = 30

  character(len=*), parameter :: task_name = "epsilon"

  !> Interface to the parameters defined in the input file
  type task_epsilon_parameters
    private
    type(kpoints_sets) :: q_points
    logical :: usingIrreducibleWedge
    logical :: buildFromPolarizability
    integer(i32) :: n_omega
    character(len=max_length) :: output_format
    real(dp) :: eigenvalue_cutoff_Coulomb_matrix
  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_qpt )
  class(task_epsilon_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_qpt
  
  call this%sanity_checks( gw_inp )
  call this%q_points%parse_input( gw_inp%taskGroup%epsilon%qpointsarray, n_qpt )
  this%usingIrreducibleWedge = gw_inp%taskGroup%epsilon%usingIrreducibleWedge
  this%buildFromPolarizability = gw_inp%taskGroup%epsilon%buildFromPolarizability
  this%n_omega = gw_inp%freqgrid%nomeg
  this%output_format = trim( adjustl( gw_inp%taskGroup%outputFormat ) )
  this%eigenvalue_cutoff_Coulomb_matrix = gw_inp%barecoul%barcevtol
end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  class(task_epsilon_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%epsilon), &
    'Element epsilon must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%barecoul), &
    'Element barecoul must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%scrcoul), &
    'Element scrcoul must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%freqgrid), &
    'Element freqgrid must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( trim(gw_inp%scrcoul%scrtype)=='rpa', &
    'Only the RPA method is currently supported when executing '//'"'//task_name//'"' )
end subroutine


!> Subroutine to be invoked when task `epsilon` must be executed
subroutine execute_task_epsilon( n_qpoints_max, file_format )
  integer(i32), intent(in) :: n_qpoints_max
  character(len=*), intent(in) :: file_format

  integer(i32) :: iq, iq_reducible, iq_io
  integer(i32) :: i, i_start, i_end
  integer(i32) :: omega_i, omega_f
  real(dp) :: eigenvalue_cutoff
  type(task_epsilon_parameters) :: input_parameters
  type(mpiinfo) :: mpi_environment_qpoints
  
  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  call input_parameters%parse_input( input%gw, n_qpoints_max )
  call input_parameters%q_points%obtain_list_of_indexes()
  call mpi_environment_qpoints%init( mpiglobal%comm )
  call distribute_loop( mpi_environment_qpoints, size(input_parameters%q_points%list_of_indexes), &
    i_start, i_end )
  
  omega_i = 1
  omega_f = input_parameters%n_omega
  eigenvalue_cutoff = max( real_zero, input_parameters%eigenvalue_cutoff_Coulomb_matrix )
  
  ! Attention: calcpmatgw makes use of MPI parallelization and calls a mpi_barrier
  ! In the case of building the dielectric matrix from polarizability
  ! the momentum transfer matrix elements should already be present in the folder, 
  ! as they are required by the polarizability task 
  if( isGammaInList( kqset%vqc(:,input_parameters%q_points%list_of_indexes) ) .and. & 
    .not. input_parameters%buildFromPolarizability ) call calcpmatgw
  do i = i_start, i_end
    iq = input_parameters%q_points%list_of_indexes(i) ! This refers always to the list either full or irreducible
    if (input_parameters%usingIrreducibleWedge) then
      if( mpiglobal%rank == 0) call write_to_gwinfo( '('//task_name//'): q-point cycle, iq (irreducible) = ' // to_char(iq) )
      iq_reducible = kset%ikp2ik(iq) ! iq is the index of the irreducible q-point; iq_reducible is the index in the reducible q-point list
      iq_io = iq ! We are outputing the files with the irreducible wedge numbering
    else
      if( mpiglobal%rank == 0) call write_to_gwinfo( '('//task_name//'): q-point cycle, iq = ' // to_char(iq) )
      iq_reducible = iq ! iq is the index in the full BZ
      iq_io = iq_reducible ! We are outputing files with full BZ numbering
    end if

    call read_sgi_from_file( iq_reducible, file_format )
    call calcmpwipw( iq_reducible )
    call read_barcev_vmat_from_file( iq_reducible, file_format )
    Gamma = gammapoint( kqset%vqc(:, iq_reducible), tol=1.e-6_dp )
    call calculate_sqrt_bare_coulomb( iq_reducible, eigenvalue_cutoff, Gamma )
    call init_dielectric_function( mbsiz, omega_i, omega_f, Gamma )
    if (.not. input_parameters%buildFromPolarizability ) then
      call calcepsilon( iq_reducible, omega_i, omega_f )
    else
      call init_polarizability( matsiz, omega_i, omega_f, .false.)
      call read_polarizability_from_file(iq_io, Gamma, file_format, input_parameters%usingIrreducibleWedge)
      call from_polarizability_to_epsilon(iq_reducible, Gamma, omega_i, omega_f)
      call delete_polarizability()
    endif
    call write_epsilon_to_file( iq_io, Gamma, file_format, input_parameters%usingIrreducibleWedge)
  end do

  call deallocate_global_arrays
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


!> Deallocate global arrays needed to obtain the dielectric matrix
subroutine deallocate_global_arrays
  call delete_dielectric_function( Gamma=.true. )
  call delete_coulomb_potential
end subroutine

end module
