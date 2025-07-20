!> This module contains classes and subroutines needed to execute the 
!> task `epsilon`, used as an element of `taskGroup` in `gw`
module task_epsilon
  use asserts, only: assert
  use constants, only: real_zero
  use calculate_dielectric_function, only: calcepsilon, epsilon_indexes
  use exciting_mpi, only: mpiinfo, xmpi_gather
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage, write_to_gwinfo_parallelization_info, &
    write_to_gwinfo_progress, progress_calculate_epsilon, progress_coulomb, progress_ipw, progress_sgi, &
    progress_write_epsilon, write_to_gwinfo_table_with_index_map, &
    q_points_numbering, q_points_numbering_abbr, q_points_indexes, q_points_indexes_abbr, &
    k_points_indexes, k_points_indexes_abbr, empty_bands_indexes, empty_bands_abbr
  use mod_coulomb_potential, only: barc, delete_coulomb_potential, read_barcev_vmat_from_file, calculate_sqrt_bare_coulomb
  use mod_dielectric_function, only: write_epsilon_to_file, init_dielectric_function, delete_dielectric_function
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_mpi_gw, only: mpi_domain, indexes_parallelization, pack_parallelization_indexes, &
    unpack_parallelization_indexes, define_mpi_domains
  use mod_polarizability, only: init_polarizability, from_polarizability_to_epsilon, read_polarizability_from_file, delete_polarizability
  use mod_product_basis, only: mbsiz, matsiz, read_sgi_from_file, mpwipw
  use mod_selfenergy, only: singc1, singc2
  use modinput, only: input, gw_type
  use modmpi, only: terminate_if_false, mpiglobal
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
#include "offload.fpp"

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
    !> Number of MPI Domains to split over q-points
    integer(i32) :: n_MPI_Domains_qpoints
    !> Number of MPI Domains to split over k-points
    integer(i32) :: n_MPI_Domains_kpoints
    !> If true, print the polarizability factor \(F_{nm}(\mathbf{k},\mathbf{q},\omega)\)
    logical :: print_Polarizability_Factor
    real(dp) :: eigenvalue_cutoff_Coulomb_matrix
  contains
    procedure :: parse_input, sanity_checks
  end type
contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_qpt )
  class(task_epsilon_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_qpt

  call this%sanity_checks( gw_inp )
  call this%q_points%parse_input( gw_inp%taskGroup%epsilon%qpointsarray, n_qpt )
  this%usingIrreducibleWedge = gw_inp%taskGroup%epsilon%usingIrreducibleWedge
  this%buildFromPolarizability = gw_inp%taskGroup%epsilon%buildFromPolarizability
  this%n_omega = gw_inp%freqgrid%nomeg
  this%output_format = trim( adjustl( gw_inp%taskGroup%outputFormat ) )
  this%n_MPI_Domains_qpoints = gw_inp%taskGroup%epsilon%MPIDomainsQpoints
  this%n_MPI_Domains_kpoints = gw_inp%taskGroup%epsilon%MPIDomainsKpoints
  this%print_Polarizability_Factor = gw_inp%taskGroup%epsilon%printPolarizabilityFactor
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
  call terminate_if_false( gw_inp%taskGroup%epsilon%MPIDomainsQpoints > 0, &
    'MPIDomainsQpoints must be positive' )
  call terminate_if_false( gw_inp%taskGroup%epsilon%MPIDomainsKpoints > 0, &
    'MPIDomainsKpoints must be positive' )
end subroutine


!> Subroutine to be invoked when task `epsilon` must be executed
subroutine execute_task_epsilon( all_q_points, idx_reduced_qpt, n_kpoints, first_empty_state, last_empty_state, file_format, dry_run )
  !> Array containing all q-points
  real(dp), intent(in) :: all_q_points(:, :)
  !> Array containing all reduced q-points needed in this calculation
  integer(i32), intent(in) :: idx_reduced_qpt(:)
  !> Number of k-points
  integer(i32), intent(in) :: n_kpoints
  !> Index of the first unoccupied state
  integer(i32), intent(in) :: first_empty_state
  !> Index of the last unoccupied state
  integer(i32), intent(in) :: last_empty_state
  !> Format of input/output files
  character(len=*), intent(in) :: file_format
  !> If `.true.`, a dry-run only is required
  logical, intent(in) :: dry_run

  integer(i32) :: iq, iq_reducible, iq_io, i, omega_i, omega_f
  integer(i32) :: rank_to_write
  real(dp), parameter :: tol = 1.0e-6_dp
  logical :: myrank_writes_GWINFO, myrank_writes_EPSILON
  logical :: write_progress
  real(dp) :: eigenvalue_cutoff
  type(task_epsilon_parameters) :: input_parameters
  type(mpi_domain) :: mpi_qpoints, mpi_kpoints
  type(indexes_parallelization) :: empty_states
  type(indexes_parallelization) :: frequencies
  
  call assert( size( all_q_points, 1 ) == 3, 'qpoints must have size 3 along 1st dimension' )
  rank_to_write = mpiglobal%root
  myrank_writes_GWINFO = ( mpiglobal%rank == rank_to_write )
  write_progress = myrank_writes_GWINFO
  if( myrank_writes_GWINFO ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )

  call input_parameters%parse_input( input%gw, size( idx_reduced_qpt ) )
  call input_parameters%q_points%obtain_list_of_indexes()
  if( myrank_writes_GWINFO ) call write_to_gwinfo_table_with_index_map( &
    input_parameters%q_points%list_of_indexes, &
    [character(len=max(len(q_points_numbering_abbr), len(q_points_indexes_abbr)))::q_points_numbering_abbr, q_points_indexes_abbr], &
    [character(len=max(len(q_points_numbering), len(q_points_indexes)))::q_points_numbering, q_points_indexes] )
  
  call mpi_qpoints%index%set_global_first_last( 1, size(input_parameters%q_points%list_of_indexes) )
  call mpi_kpoints%index%set_global_first_last( 1, n_kpoints )
  call empty_states%set_global_first_last( first_empty_state, last_empty_state )
  call define_mpi_domains( mpiglobal, input_parameters%n_MPI_Domains_qpoints, input_parameters%n_MPI_Domains_kpoints,&
    mpi_qpoints, mpi_kpoints, empty_states )
  call write_parallelization_info( mpiglobal, mpi_qpoints%index, mpi_kpoints%index, &
    empty_states, rank_to_write )

  myrank_writes_EPSILON = ( mpi_qpoints%mpi_environment%rank == rank_to_write )
  omega_i = 1
  omega_f = input_parameters%n_omega
  frequencies = indexes_parallelization( omega_i, omega_f, 1, input_parameters%n_omega )
  eigenvalue_cutoff = max( real_zero, input_parameters%eigenvalue_cutoff_Coulomb_matrix )

  if( dry_run ) then
    if( myrank_writes_GWINFO ) call write_dry_run_info( last_empty_state-first_empty_state+1, input_parameters%n_omega )
  else
    ! Attention: calcpmatgw makes use of MPI parallelization and calls a mpi_barrier
    ! In the case of building the dielectric matrix from polarizability
    ! the momentum transfer matrix elements should already be present in the folder, 
    ! as they are required by the polarizability task 
    if( isGammaInList( all_q_points(:,input_parameters%q_points%list_of_indexes) ) .and. & 
      .not. input_parameters%buildFromPolarizability ) call calcpmatgw
    do i = mpi_qpoints%index%my_first, mpi_qpoints%index%my_last
      iq = input_parameters%q_points%list_of_indexes(i)
      if (input_parameters%usingIrreducibleWedge) then
        iq_reducible = idx_reduced_qpt(iq) ! iq is the index of the irreducible q-point; iq_reducible is the index in the reducible q-point list
        iq_io = iq ! We are outputing the files with the irreducible wedge numbering
        if( myrank_writes_GWINFO ) then
          call write_to_gwinfo( '('//task_name//'): q-point cycle, '// &
            trim(q_points_numbering_abbr) // ' = ' // to_char(i) // ', ' // &
            trim(q_points_indexes_abbr) // ' (irreducible) = ' // to_char(iq) )
          call write_to_gwinfo_progress( progress_sgi )
        end if
      else
        iq_reducible = idx_reduced_qpt(iq) ! iq is the index in the full BZ
        iq_io = iq_reducible ! We are outputing files with full BZ 
        if( myrank_writes_GWINFO ) then
          call write_to_gwinfo( '('//task_name//'): q-point cycle, '// &
            trim(q_points_numbering_abbr) // ' = ' // to_char(i) // ', ' // &
            trim(q_points_indexes_abbr) // ' = ' // to_char(iq) )
          call write_to_gwinfo_progress( progress_sgi )
        end if
      end if
      call read_sgi_from_file( iq_reducible, file_format )
      if( myrank_writes_GWINFO ) call write_to_gwinfo_progress( progress_ipw )
      call calcmpwipw( iq_reducible )
      if( myrank_writes_GWINFO ) call write_to_gwinfo_progress( progress_coulomb )
      call read_barcev_vmat_from_file( iq_reducible, file_format )
      Gamma = gammapoint( all_q_points(:, iq_reducible), tol=tol )
      call calculate_sqrt_bare_coulomb( iq_reducible, eigenvalue_cutoff, Gamma )
      call init_dielectric_function( mbsiz, omega_i, omega_f, Gamma )
      if( myrank_writes_GWINFO ) call write_to_gwinfo_progress( progress_calculate_epsilon )
      if (.not. input_parameters%buildFromPolarizability ) then
        call calcepsilon( iq_reducible, epsilon_indexes( mpi_kpoints%index, empty_states, frequencies ), &
          write_progress, mpi_qpoints%mpi_environment, input_parameters%print_Polarizability_Factor, file_format )
      else
        call init_polarizability( matsiz, omega_i, omega_f, .false.)
        call read_polarizability_from_file(iq_io, Gamma, file_format, input_parameters%usingIrreducibleWedge)
        call from_polarizability_to_epsilon(iq_reducible, Gamma, omega_i, omega_f)
        call delete_polarizability()
      end if
      if( myrank_writes_GWINFO ) call write_to_gwinfo_progress( progress_write_epsilon )
      if( myrank_writes_EPSILON ) call write_epsilon_to_file( iq_io, Gamma, file_format, input_parameters%usingIrreducibleWedge )
      if( myrank_writes_GWINFO ) call write_to_gwinfo('')
    end do
  end if

  call deallocate_global_arrays
end subroutine


!> (private) Write parallelization info to `GW_INFO.OUT`
subroutine write_parallelization_info( mpi_global, my_qpoints, my_kpoints, my_bands, rank_to_write_GWINFO )
  !> Type with the global MPI environment
  type(mpiinfo), intent(in) :: mpi_global
  !> Type with the first and last q-point indexes of this rank
  type(indexes_parallelization), intent(in) :: my_qpoints
  !> Type with the first and last k-point indexes of this rank
  type(indexes_parallelization), intent(in) :: my_kpoints
  !> Type with the first and last band indexes of this rank
  type(indexes_parallelization), intent(in) :: my_bands
  !> The rank that writes into `GW_INFO.OUT`
  integer(i32), intent(in) :: rank_to_write_GWINFO
  
  integer(i32) :: i, n_procs
  integer(i32), parameter :: n_variables_represented = 3 !q-points, k-points, and empty bands
  integer(i32), allocatable :: indexes_from_all_procs(:)
  integer(i32), allocatable :: ranks(:)

  n_procs = mpi_global%procs
  call xmpi_gather( mpi_global, pack_parallelization_indexes([my_qpoints, my_kpoints, my_bands]), &
    indexes_from_all_procs )
  if( mpi_global%rank == rank_to_write_GWINFO ) then 
    ranks = [ (i-1, i = 1, n_procs) ]
    call write_to_gwinfo_parallelization_info( ranks, &
      reshape(unpack_parallelization_indexes(indexes_from_all_procs), shape=[n_variables_represented,n_procs]), &
      [character(max(len(q_points_numbering_abbr),len(k_points_indexes_abbr),len(empty_bands_abbr)))::q_points_numbering_abbr,k_points_indexes_abbr,empty_bands_abbr], &
      [character(max(len(q_points_indexes),len(k_points_indexes),len(empty_bands_indexes)))::q_points_indexes,k_points_indexes,empty_bands_indexes] )
  end if
end subroutine


!> Check if the gamma point is among the q-points
pure logical function isGammaInList( q_points )
  !> List of q-points (first dimension has the 3 cartesian coordinates)
  real(dp), intent(in) :: q_points(:, :)

  integer(i32) :: i

  isGammaInList = .false.
  do i = 1, size( q_points, 2 )
    isGammaInList = gammapoint( q_points(:, i), tol=1.e-6_dp )
    if( isGammaInList ) exit
  end do
end function

!> Deallocate global arrays needed to obtain the dielectric matrix
subroutine deallocate_global_arrays
  call delete_dielectric_function( Gamma=.true. )
  call delete_coulomb_potential
  OMP_OFFLOAD target exit data map(delete: barc) if(allocated(barc))
  if (allocated(barc)) deallocate(barc)
  OMP_OFFLOAD target exit data map(delete: mpwipw) if(allocated(mpwipw))
  if (allocated(mpwipw)) deallocate(mpwipw)
end subroutine

!> (private) Write general information about memory usage in the context of a dry run
subroutine write_dry_run_info( n_empty_states, n_omega )
  use gw_memory, only: close_file_memory_usage, open_file_memory_usage, write_memory_usage, &
    field_and_values, field_sgi, field_mpwipw, field_vmat, field_epsilon, field_fnm, field_minmmat, &
    field_minm, field_temp_calcminm2
  use mod_bands, only: n_occupied_bands => nomax
  use mod_core_states, only: n_core_states => ncg
  use mod_product_basis, only: matsizmax
  use modgw, only: Gkqset, Gqbarc, Gqset, kqset, mblksiz
  use modinput, only: input

  integer(i32), intent(in) :: n_empty_states
  integer(i32), intent(in) :: n_omega

  integer(i32), parameter :: mega_byte = 1024**2
  integer(i32) :: n_dim, m_dim, ngq_block_size
  integer(i32), parameter :: n_descriptors = 8
  type(field_and_values) :: descriptors(n_descriptors)

  call open_file_memory_usage( )
  call descriptors(1)%init( field_sgi, [Gqset%ngkmax, Gqset%ngkmax], mega_byte )
  call descriptors(2)%init( field_mpwipw, [Gqset%ngkmax, Gqbarc%ngkmax], mega_byte )
  call descriptors(3)%init( field_vmat, [matsizmax, matsizmax], mega_byte )
  call descriptors(4)%init( field_epsilon, [matsizmax, matsizmax, n_omega], mega_byte )
  n_dim = n_occupied_bands
  if (input%gw%coreflag=='all') n_dim = n_dim + n_core_states
  call descriptors(5)%init( field_fnm, [n_dim, n_empty_states, n_omega, kqset%nkpt], mega_byte ) 
  m_dim = min(n_empty_states, mblksiz)
  call descriptors(6)%init( field_minmmat, [matsizmax, n_dim, m_dim], mega_byte ) 
  call descriptors(7)%init( field_minm, [matsizmax, n_dim, m_dim], mega_byte ) 
  ngq_block_size = merge(Gqset%ngkmax, min(input%gw%GBatchCount, Gqset%ngkmax), input%gw%GBatchCount <= 0)
  call descriptors(8)%init( field_temp_calcminm2, [Gkqset%ngkmax, Gkqset%ngkmax, ngq_block_size], mega_byte ) 
  call write_memory_usage( "task = " // task_name // ", Memory estimates (in MB)", descriptors )
  call close_file_memory_usage( )
end subroutine

end module
