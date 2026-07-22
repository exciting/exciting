!> This module contains classes and subroutines needed to compute input
!> for the coupled cluster code cc4s.
module task_cc4sInterface

  use constants, only: real_zero
  use exciting_mpi, only: mpiinfo
  use gw_info, only: write_to_gwinfo_boxmessage
  use modgw, only: kqset, Gqset
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, mpiglobal
  use mod_kqpts, only: kpoints_sets
  use modgw,                          only: mbsiz
  use mod_misc_gw, only: gammapoint
  use mod_product_basis, only: locmatsiz, matsiz, read_sgi_from_file
  use mod_eigenvalue_occupancy, only: efermi
  use mod_cc4sInterface, only: init_coulomb_vertex, delete_coulomb_vertex, &
                               compute_coulomb_vertex, prepare_dft_eigenvalues, &
                               write_coulomb_vertex_info_to_yaml, write_coulomb_vertex_to_file, &
                               write_eigenenergies_to_yaml, write_scf_energies_to_file, &
                               write_orbital_properties_to_yaml
  use mod_cc4sInterface, only: file_name_coulomb_vertex, yaml_name_coulomb_vertex, &
                               file_name_eigen_energies, yaml_name_eigen_energies, &
                               yaml_name_orbital_properties, coulomb_vertex, &
                               write_test_output
  use mod_coulomb_potential, only: delete_coulomb_potential, calculate_bare_coulomb, &
                                   calculate_sqrt_bare_coulomb                            
  use precision, only: dp, i32
  use mod_bands, only: nstdf, nomax


  implicit none

  private

  public :: execute_task_cc4sInterface

  character(len=*), parameter :: task_name = "cc4sInterface"

  !> Interface to the parameters defined in the input file
  type task_cc4sInterface_parameters
    private
    type(kpoints_sets) :: q_points 
    type(kpoints_sets) :: k_points 
    integer(i32) :: n_omega 
    integer(i32) :: n_last_frozen
    real(dp) :: eigenvalue_cutoff_Coulomb_matrix
    logical :: test_run

  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp )
  class(task_cc4sInterface_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in):: gw_inp
  
  call this%sanity_checks( gw_inp )
  call this%q_points%parse_input( gw_inp%taskGroup%cc4sInterface%qpointsarray, 1) !For now only one?
  call this%k_points%parse_input( gw_inp%taskGroup%cc4sInterface%kpointsarray, 1 ) !For now only one?
  this%n_last_frozen = gw_inp%taskGroup%cc4sInterface%nLastFrozen
  this%eigenvalue_cutoff_Coulomb_matrix = gw_inp%barecoul%barcevtol
  this%n_omega = gw_inp%freqgrid%nomeg
  this%test_run = gw_inp%taskGroup%cc4sInterface%testRun
end subroutine


!> Perform sanity checks on the input parameters in the `gw` element
subroutine sanity_checks( this, gw_inp )
  class(task_cc4sInterface_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%cc4sInterface), &
    'Element cc4sInterface must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%freqgrid), &
    'Element freqgrid must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( associated(gw_inp%barecoul), &
    'Element barecoul must be present when executing '//'"'//task_name//'"' )

end subroutine


!> Subroutine to be invoked when task `cc4sInterface` must be executed
subroutine execute_task_cc4sInterface( file_format )
  !> format used to write the files needed by exciting
  character(len=*), intent(in) :: file_format

  integer(i32) :: iq
  integer(i32) :: omega_i, omega_f
  type(task_cc4sInterface_parameters) :: input_parameters
  logical :: is_Gamma
  real(dp) :: eigenvalue_cutoff_Coulomb_matrix
  real(dp), allocatable :: scf_eigenvalues_flattened(:)
  real(dp), allocatable :: scf_eigenvalues(:,:,:)
  integer(i32), allocatable :: sorted_ids_scf_eigvals(:)
  real(dp) :: e_fermi


  if( mpiglobal%rank == 0 ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  call input_parameters%parse_input( input%gw )
  
  omega_i = 1
  omega_f = input_parameters%n_omega

  ! only implemented for a single k- / q-point 
  call input_parameters%q_points%obtain_list_of_indexes()
  call input_parameters%k_points%obtain_list_of_indexes()
  call terminate_if_false( (size(input_parameters%k_points%list_of_indexes(:))==1), &
  'Currently only one k-point can be used when executing '//'"'//task_name//'"' )
   call terminate_if_false( (size(input_parameters%q_points%list_of_indexes(:))==1), &
  'Currently only one q-point can be used when executing '//'"'//task_name//'"' )

  iq = 1
  
  ! We need to set matsiz to the appropiate value
  matsiz = locmatsiz + Gqset%ngk(1,iq)
  call read_sgi_from_file( iq, file_format )
  call calcmpwipw( iq )
  call calculate_bare_coulomb( iq )
  is_Gamma = gammapoint( kqset%vqc(:, iq), tol=1.e-6_dp )
  eigenvalue_cutoff_Coulomb_matrix = max( real_zero, input_parameters%eigenvalue_cutoff_Coulomb_matrix )
  call calculate_sqrt_bare_coulomb( iq, eigenvalue_cutoff_Coulomb_matrix, is_Gamma )

  call init_coulomb_vertex(mbsiz, nstdf) 
  call compute_coulomb_vertex(iq, omega_i, omega_f, mbsiz, input_parameters%n_last_frozen) 
  call prepare_dft_eigenvalues(nomax, nstdf-nomax, scf_eigenvalues, scf_eigenvalues_flattened, sorted_ids_scf_eigvals)

  if (input_parameters%test_run) then
    call write_test_output(iq, coulomb_vertex)
  else
    call write_scf_energies_to_file(scf_eigenvalues, nstdf, 1, file_name_eigen_energies)
    call write_coulomb_vertex_to_file(iq, coulomb_vertex, file_name_coulomb_vertex)
  end if 

  call write_coulomb_vertex_info_to_yaml(nstdf, mbsiz+1, 1.0_dp, yaml_name_coulomb_vertex)
  call write_orbital_properties_to_yaml(sorted_ids_scf_eigvals, nstdf, yaml_name_orbital_properties)
  call write_eigenenergies_to_yaml(1.0_dp, efermi, scf_eigenvalues_flattened, yaml_name_eigen_energies)

  call delete_coulomb_potential 
  call delete_coulomb_vertex()

end subroutine execute_task_cc4sInterface


end module task_cc4sInterface
