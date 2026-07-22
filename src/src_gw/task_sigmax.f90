!> Module for the task sigmax
module task_sigmax
#include "asserts.fpp"
  use constants, only: zzero
  use calculate_exchange_self_energy, only: calcselfx
  use exciting_mpi, only: mpiinfo, xmpi_gather
  use gw_info, only: write_to_gwinfo, write_to_gwinfo_boxmessage, write_to_gwinfo_parallelization_info, &
    write_to_gwinfo_progress_bar, write_to_gwinfo_table_with_index_map, &
    k_points_numbering, k_points_numbering_abbr, k_points_indexes, k_points_indexes_abbr, &
    q_points_indexes, q_points_indexes_abbr
  use math_utils, only: all_zero
  use modinput, only: input, gw_type
  use modmpi, only: mpiglobal, terminate_if_false
  use mod_coulomb_potential, only: barc, delete_coulomb_potential, read_barcev_vmat_from_file
  use mod_kqpts, only: kpoints_sets
  use mod_misc_gw, only: Gamma, gammapoint
  use mod_mpi_gw, only: mpi_domain, indexes_parallelization, define_mpi_domains, &
    mpi_sum_array, pack_parallelization_indexes, unpack_parallelization_indexes
  use mod_product_basis, only: read_sgi_from_file, mpwipw
  use mod_selfenergy, only: selfex, write_selfex_single_kpoint
  use mod_offdiagonal_selfenergy, only: init_offdiagonal_selfenergy_exchange, delete_offdiagonal_selfenergy, &
                                        mpi_reduce_offdiagonal_selfenergy_exchange, write_offdiagonal_selfenergy_exchange
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
#include "offload.fpp"
  implicit none

  private

  public :: execute_task_sigmax

  character(len=*), parameter :: task_name = "sigmax"

  !> Interface to the parameters defined in the input file
  type task_sigmax_parameters
    private
    type(kpoints_sets) :: k_points
    !> Number of MPI Domains to split over k-points
    integer(i32) :: n_MPI_Domains_kpoints
    !> Flag determining if the offdiagonal terms of the exchange
    !> self-energy are computed
    logical :: offdiagonal
  contains
    procedure :: parse_input, sanity_checks
  end type

contains

!> Interface to the parameters defined in the input file
subroutine parse_input( this, gw_inp, n_kpt )
  class(task_sigmax_parameters), intent(inout) :: this
  !> type with the variables given in the input file
  type(gw_type), intent(in) :: gw_inp
  !> maximum number of q-points
  integer(i32), intent(in) :: n_kpt

  call this%sanity_checks( gw_inp )
  call this%k_points%parse_input( gw_inp%taskGroup%sigmax%kpointsarray, n_kpt )
  this%n_MPI_Domains_kpoints = gw_inp%taskGroup%sigmax%MPIDomainsKpoints
  this%offdiagonal = gw_inp%taskGroup%sigmax%offdiagonal 

end subroutine


subroutine sanity_checks( this, gw_inp )
  class(task_sigmax_parameters), intent(in) :: this
  !> type with the variables given in the input file
  type(gw_type):: gw_inp

  call terminate_if_false( associated(gw_inp%taskGroup%sigmax), &
    'Element sigmax must be present when executing '//'"'//task_name//'"' )
  call terminate_if_false(gw_inp%taskGroup%sigmax%MPIDomainsKpoints > 0, &
    'MPIDomainsKpoints must be positive when executing "' // task_name // '"')

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

  integer(i32) :: i, ik, iq, n_qpoints, rank_to_write
  real(dp), parameter :: tolerance_zero_vector = 1.e-6_dp
  real(dp) :: ti, tf, t_acc
  logical :: myrank_writes_GWINFO, myrank_writes_SIGMAX
  type(task_sigmax_parameters) :: input_parameters
  type(mpi_domain) :: mpi_kpoints
  type(indexes_parallelization) :: q_points

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
  call q_points%set_global_first_last( 1, n_qpoints )
  call define_mpi_domains( mpiglobal, input_parameters%n_MPI_Domains_kpoints, mpi_kpoints, q_points )
  call write_parallelization_info( mpiglobal, mpi_kpoints%index, q_points, rank_to_write )
  myrank_writes_SIGMAX = ( mpi_kpoints%mpi_environment%rank == rank_to_write )
  
  do i = mpi_kpoints%index%my_first, mpi_kpoints%index%my_last
    ik = input_parameters%k_points%list_of_indexes(i)
    if( myrank_writes_GWINFO ) then
      call write_to_gwinfo( '('//task_name//'): k-point cycle, '// &
      trim(k_points_numbering_abbr) // ' = ' // to_char(i) // ', ' // &
      trim(k_points_indexes_abbr) // ' = ' // to_char(ik) )
    end if
    if( allocated(selfex) ) deallocate( selfex )
    allocate( selfex(first_band:last_band, ik:ik ), source=zzero )
    if (input_parameters%offdiagonal) call init_offdiagonal_selfenergy_exchange(ik, ik)

    t_acc = 0.0_dp
    call timesec(tf)
    do iq = q_points%my_first, q_points%my_last
      ti = tf
      call read_barcev_vmat_from_file( iq, file_format )
      call read_sgi_from_file( iq, file_format )
      call calcmpwipw( iq )
      Gamma = gammapoint( qpoints(:, iq), tol=tolerance_zero_vector )
      call calcselfx( iq, ik, ik, input_parameters%offdiagonal )
      call timesec(tf)
      if( myrank_writes_GWINFO ) then
        t_acc = t_acc + (tf-ti)
        call write_to_gwinfo_progress_bar( t_acc, real(iq-q_points%my_first+1,dp)/(q_points%my_last-q_points%my_first+1), identation_level=2 )
      end if
    end do
    call mpi_sum_array( selfex, mpi_kpoints%mpi_environment, all_reduce=.false. )
    if (input_parameters%offdiagonal) call mpi_reduce_offdiagonal_selfenergy_exchange(mpi_kpoints%mpi_environment)

    if( myrank_writes_SIGMAX ) call write_selfex_single_kpoint( ik, file_format )
    if( myrank_writes_SIGMAX .and. input_parameters%offdiagonal) call write_offdiagonal_selfenergy_exchange( ik, file_format )
    if( myrank_writes_GWINFO ) call write_to_gwinfo( '' )
    if (input_parameters%offdiagonal) call delete_offdiagonal_selfenergy()
  end do
  
  ! Clean global variables
  call delete_coulomb_potential
  OMP_OFFLOAD target exit data map(delete: barc) if(allocated(barc))
  if (allocated(barc)) deallocate(barc)
  OMP_OFFLOAD target exit data map(delete: mpwipw) if(allocated(mpwipw))
  if (allocated(mpwipw)) deallocate(mpwipw)


end subroutine


!> (private) Write parallelization info to `GW_INFO.OUT`
subroutine write_parallelization_info( mpi_global, my_kpoints, my_qpoints, rank_to_write_GWINFO )
  !> Type with the global MPI environment
  type(mpiinfo), intent(in) :: mpi_global
  !> Type with the first and last k-point indexes of this rank
  type(indexes_parallelization), intent(in) :: my_kpoints
    !> Type with the first and last q-point indexes of this rank
  type(indexes_parallelization), intent(in) :: my_qpoints
  !> The rank that writes into `GW_INFO.OUT`
  integer(i32), intent(in) :: rank_to_write_GWINFO
  
  integer(i32) :: i, n_procs
  integer(i32), parameter :: n_variables_represented = 2 !k-points and q-points
  integer(i32), allocatable :: indexes_from_all_procs(:)
  integer(i32), allocatable :: ranks(:)
  type(indexes_parallelization), allocatable :: set_with_indexes(:, :)
  character(len=:), allocatable :: symbols(:)
  character(len=:), allocatable :: names(:)

  n_procs = mpi_global%procs
  ! All MPI procs send their [ki, kf], [qi, qf] ranges to rank 0
  call xmpi_gather( mpi_global, pack_parallelization_indexes([my_kpoints, my_qpoints]), &
    indexes_from_all_procs )
  if( mpi_global%rank == rank_to_write_GWINFO ) then 
    ranks = [ (i-1, i = 1, n_procs) ]
    ! Convert [ki, kf] and [qi, qf] ranges from all MPI procs. (`indexes_from_all_procs`)
    ! into a matrix of type `indexes_parallelization`.
    set_with_indexes = reshape( unpack_parallelization_indexes(indexes_from_all_procs), shape=[n_variables_represented,n_procs] )
    ! Create an array containing `k_points_numbering_abbr` and `q_points_indexes_abbr`,
    ! using the larger of the two string lengths to define each string size in the array
    symbols = [character(max(len(k_points_numbering_abbr),len(q_points_indexes_abbr)))::k_points_numbering_abbr, q_points_indexes_abbr]
    ! Perform the same operation as for `symbols`, but using `k_points_numbering` and `q_points_indexes` instead
    names = [character(max(len(k_points_numbering),len(q_points_indexes))):: k_points_numbering, q_points_indexes]
    call write_to_gwinfo_parallelization_info( ranks, set_with_indexes, symbols, names )
  end if
end subroutine

end module
