module task_sigmac
#include "asserts.fpp"
  use calculate_correlation_self_energy, only: calcselfc, sigmac_indexes
  use constants, only: zzero, real_zero
  use exciting_mpi, only: mpiinfo, xmpi_gather
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage, write_to_gwinfo_parallelization_info, &
    write_to_gwinfo_progress_bar, write_to_gwinfo_table_with_index_map, &
    k_points_numbering, k_points_numbering_abbr, k_points_indexes, k_points_indexes_abbr, &
    q_points_indexes, q_points_indexes_abbr, bands_indexes, bands_abbr
  use math_utils, only: all_zero
  use modinput, only: input, gw_type
  use modgw, only: freq, ibgw, nbgw
  use modmpi, only: mpiglobal, terminate_if_false, distribute_loop
  use mod_coulomb_potential, only: read_barcev_vmat_from_file, calculate_sqrt_bare_coulomb, delete_coulomb_potential, barc
  use mod_dielectric_function, only: read_inverse_epsilon_from_file, epsilon
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_mpi_gw, only: define_mpi_domains, indexes_parallelization, mpi_domain, &
    mpi_sum_array, pack_parallelization_indexes, unpack_parallelization_indexes
  use mod_product_basis, only: read_sgi_from_file, mpwipw
  use mod_selfenergy, only: selfec, write_selfec_single_kpoint, & 
    generate_frequency_grid_for_correlation_self_energy
  use precision, only: i32, dp
  use to_char_conversion, only: to_char
  use mod_offdiagonal_selfenergy, only: init_offdiagonal_selfenergy_correlation, mpi_reduce_offdiagonal_selfenergy_correlation, &
                                        write_offdiagonal_selfenergy_correlation, delete_offdiagonal_selfenergy                                   
#include "offload.fpp"

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
    !> Number of MPI Domains to split over q-points
    integer(i32) :: n_MPI_Domains_qpoints
    !> Number of MPI Domains to split over k-points
    integer(i32) :: n_MPI_Domains_kpoints
    !> Flag determining if the offdiagonal terms of the correlation
    !> self-energy are computed
    logical :: offdiagonal
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
  this%n_MPI_Domains_kpoints = gw_inp%taskGroup%sigmac%MPIDomainsKpoints
  this%n_MPI_Domains_qpoints = gw_inp%taskGroup%sigmac%MPIDomainsQpoints
  this%eigenvalue_cutoff_Coulomb_matrix = gw_inp%barecoul%barcevtol
  this%offdiagonal = gw_inp%taskGroup%sigmac%offdiagonal 
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
  call terminate_if_false( associated(gw_inp%selfenergy), &
    'Element selfenergy must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false( gw_inp%selfenergy%method == 'ac', &
    'Task "' // task_name // '" only implemented for ac as selfenergy method' )
  call terminate_if_false(gw_inp%taskGroup%sigmac%MPIDomainsKpoints > 0, &
    'MPIDomainsKpoints must be positive when executing "' // task_name // '"')
  call terminate_if_false(gw_inp%taskGroup%sigmac%MPIDomainsQpoints > 0, &
    'MPIDomainsQpoints must be positive when executing "' // task_name // '"')
end subroutine


subroutine sanity_check_epsilon_and_barc()
  call terminate_if_false( size( barc, 2 ) == size( epsilon, 1 ), &
    'Coulomb matrix and the inverse of epsilon have incompatible sizes' )
end subroutine

subroutine sanity_check_frequencies_of_epsilon()
  call terminate_if_false( freq%nomeg == size( epsilon, 3 ), &
    'Number of frequency points read for the inverse of epsilon incompatible with current calculation' )
end subroutine


subroutine execute_task_sigmac( n_kpoints_max, qpoints, first_state, last_state, file_format )
  integer(i32), intent(in) :: n_kpoints_max
  !> List of q-points in cartesian coordinates
  real(dp), intent(in) :: qpoints(:, :)
  !> Index of the first state. It should usually be equal to 1
  integer(i32), intent(in) :: first_state
  !> Index of the last unoccupied state. If larger than `nempty`, then core states are included in the calculation
  integer(i32), intent(in) :: last_state
  character(len=*), intent(in) :: file_format

  integer(i32) :: ik, i, omega_i, omega_f
  integer(i32) :: iq, iq_start, iq_end, n_qpoints, rank_to_write 
  logical :: myrank_writes_GWINFO, myrank_writes_SIGMAC
  real(dp), parameter :: tolerance_zero_vector = 1.e-6_dp
  real(dp) :: ti, tf, t_acc, fraction
  real(dp) :: eigenvalue_cutoff
  type(task_sigmac_parameters) :: input_parameters
  type(mpi_domain) :: mpi_qpoints, mpi_kpoints
  type(indexes_parallelization) :: bands
  
  CALL_ASSERT( size( qpoints, 1 ) == 3, 'qpoints must have size 3 along 1st dimension' )

  rank_to_write = mpiglobal%root
  myrank_writes_GWINFO = ( mpiglobal%rank == rank_to_write )
  if( myrank_writes_GWINFO ) call write_to_gwinfo_boxmessage( '=', 'task: '//task_name )
  
  call input_parameters%parse_input( input%gw, n_kpoints_max )
  call input_parameters%k_points%obtain_list_of_indexes()
  if( myrank_writes_GWINFO ) call write_to_gwinfo_table_with_index_map( &
    input_parameters%k_points%list_of_indexes, &
    [character(len=max(len(k_points_numbering_abbr), len(k_points_indexes_abbr)))::k_points_numbering_abbr, k_points_indexes_abbr], &
    [character(len=max(len(k_points_numbering), len(k_points_indexes)))::k_points_numbering, k_points_indexes] )

  call mpi_kpoints%index%set_global_first_last( 1, size(input_parameters%k_points%list_of_indexes) )
  n_qpoints = size( qpoints, 2 )
  call mpi_qpoints%index%set_global_first_last( 1, n_qpoints )
  call bands%set_global_first_last( first_state, last_state )
  call define_mpi_domains( mpiglobal, input_parameters%n_MPI_Domains_kpoints, input_parameters%n_MPI_Domains_qpoints,&
    mpi_kpoints, mpi_qpoints, bands )
  call write_parallelization_info( mpiglobal, mpi_kpoints%index, mpi_qpoints%index, &
    bands, rank_to_write )

  myrank_writes_SIGMAC = ( mpi_kpoints%mpi_environment%rank == rank_to_write )
  omega_i = 1
  omega_f = input_parameters%n_omega
  
  iq_start = mpi_qpoints%index%my_first
  iq_end = mpi_qpoints%index%my_last

  eigenvalue_cutoff = max( real_zero, input_parameters%eigenvalue_cutoff_Coulomb_matrix )
  
  call generate_frequency_grid_for_correlation_self_energy( input%gw )
  do i = mpi_kpoints%index%my_first, mpi_kpoints%index%my_last
    ik = input_parameters%k_points%list_of_indexes(i)
    if( myrank_writes_GWINFO ) then
      call write_to_gwinfo( '('//task_name//'): k-point cycle, '// &
        trim(k_points_numbering_abbr) // ' = ' // to_char(i) // ', ' // &
        trim(k_points_indexes_abbr) // ' = ' // to_char(ik) )
    end if
    if( allocated( selfec )) deallocate( selfec )
    allocate( selfec(ibgw:nbgw, omega_i:omega_f, ik:ik), source=zzero )
    if (input_parameters%offdiagonal) call init_offdiagonal_selfenergy_correlation(ik, ik, omega_i, omega_f)
    call timesec(tf)
    t_acc = 0._dp
    do iq = iq_start, iq_end
      ti = tf
      call read_sgi_from_file( iq, file_format )
      call calcmpwipw( iq )
      call read_barcev_vmat_from_file( iq, file_format )
      Gamma = gammapoint( qpoints(:,iq), tol=tolerance_zero_vector )
      call calculate_sqrt_bare_coulomb( iq, eigenvalue_cutoff, Gamma )
      call read_inverse_epsilon_from_file( iq, Gamma, file_format )
      call sanity_check_frequencies_of_epsilon()
      call sanity_check_epsilon_and_barc()
      OMP_OFFLOAD target data map(alloc: epsilon)
      call calcselfc(iq, sigmac_indexes( indexes_parallelization(ik, ik, ik, ik), bands ), input_parameters%offdiagonal)
      OMP_OFFLOAD end target data
      call timesec(tf)
      if( myrank_writes_GWINFO ) then
        t_acc = t_acc + (tf-ti)
        fraction = real(iq - iq_start + 1, dp)/(iq_end - iq_start + 1)
        call write_to_gwinfo_progress_bar( t_acc, fraction, identation_level=1 )
      end if
    end do
    call mpi_sum_array( selfec, mpi_kpoints%mpi_environment, all_reduce=.false.)
    if (input_parameters%offdiagonal) call mpi_reduce_offdiagonal_selfenergy_correlation(mpi_kpoints%mpi_environment)
    if( myrank_writes_SIGMAC ) call write_selfec_single_kpoint( ik, file_format )
    if( myrank_writes_SIGMAC .and. input_parameters%offdiagonal) call write_offdiagonal_selfenergy_correlation( ik, file_format )
    if (input_parameters%offdiagonal) call delete_offdiagonal_selfenergy()
  end do

  ! Delete global arrays
  call delete_coulomb_potential
  OMP_OFFLOAD target exit data map(delete: barc) if(allocated(barc))
  if (allocated(barc)) deallocate(barc)
  OMP_OFFLOAD target exit data map(delete: mpwipw) if(allocated(mpwipw))
  if (allocated(mpwipw)) deallocate(mpwipw)
  OMP_OFFLOAD target exit data map(delete: epsilon) if(allocated(epsilon))
  if (allocated(epsilon)) deallocate(epsilon)

end subroutine


!> (private) Write parallelization info to `GW_INFO.OUT`
subroutine write_parallelization_info( mpi_global, my_kpoints, my_qpoints, my_bands, rank_to_write_GWINFO )
  !> Type with the global MPI environment
  type(mpiinfo), intent(in) :: mpi_global
  !> Type with the first and last k-point indexes of this rank
  type(indexes_parallelization), intent(in) :: my_kpoints
  !> Type with the first and last q-point indexes of this rank
  type(indexes_parallelization), intent(in) :: my_qpoints
  !> Type with the first and last band indexes of this rank
  type(indexes_parallelization), intent(in) :: my_bands
  !> The rank that writes into `GW_INFO.OUT`
  integer(i32), intent(in) :: rank_to_write_GWINFO
  
  integer(i32) :: i, n_procs
  integer(i32), parameter :: n_variables_represented = 3 !q-points, k-points, and bands
  integer(i32), allocatable :: indexes_from_all_procs(:)
  integer(i32), allocatable :: ranks(:)

  n_procs = mpi_global%procs
  call xmpi_gather( mpi_global, pack_parallelization_indexes([my_kpoints, my_qpoints, my_bands]), &
    indexes_from_all_procs )
  if( mpi_global%rank == rank_to_write_GWINFO ) then 
    ranks = [ (i-1, i = 1, n_procs) ]
    call write_to_gwinfo_parallelization_info( ranks, &
      reshape(unpack_parallelization_indexes(indexes_from_all_procs), shape=[n_variables_represented,n_procs]), &
      [character(max(len(k_points_numbering_abbr),len(q_points_indexes_abbr),len(bands_abbr)))::k_points_numbering_abbr,q_points_indexes_abbr,bands_abbr], &
      [character(max(len(k_points_numbering),len(q_points_indexes),len(bands_indexes)))::k_points_numbering,q_points_indexes,bands_indexes] )
  end if
end subroutine
end module
