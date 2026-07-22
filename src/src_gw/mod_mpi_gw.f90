!> Module for MPI features used by GW
module mod_mpi_gw
#include "asserts.fpp"
    use exciting_mpi, only: mpiinfo, xmpi_allreduce, xmpi_reduce
    use modmpi, only: ierr, firstofset, lastofset, distribute_loop, terminate_if_false
#ifdef MPI
    ! Remark(Ronaldo): using mpi instead of mpi_f08 leads to a seg. fault when using
    ! assumed rank arrays in calls to openmpi-MPI subroutines
    ! A similar remark is in `src/src_xs/src_rttddft/rttddft_io_parallel.f90`
    use mpi_f08, only: MPI_COMM, MPI_MAX_PROCESSOR_NAME, MPI_INTEGER, MPI_SUM, &
      MPI_COMM_RANK, MPI_BARRIER, MPI_COMM_WORLD, MPI_DOUBLE_COMPLEX, MPI_IN_PLACE, &
      MPI_Comm_split, MPI_Allreduce, MPI_Reduce, MPI_Datatype
#endif
    use precision, only: dp, i32, long_int
    use to_char_conversion, only: to_char

    implicit none

    private 

    public :: define_mpi_domains, &
              mpi_sum_array, &
              pack_parallelization_indexes, &
              unpack_parallelization_indexes

    integer(i32), public :: iqstart, iqend
    integer(i32), public :: iomstart, iomend

    !> Type to encapsulate the first and the last parallelization indexes
    !> that are treated by an MPI process
    type, public :: indexes_parallelization
      !> first index of this MPI process
      integer(i32) :: my_first
      !> last index of this MPI process
      integer(i32) :: my_last
      !> global first index
      integer(i32) :: global_first
      !> global last index
      integer(i32) :: global_last
    contains
      procedure :: set_my_first_last
      procedure :: set_global_first_last
      generic :: shift_first_last => shift_first_last_by_an_integer, shift_first_last_with_global_first
      procedure, private :: shift_first_last_by_an_integer
      procedure, private :: shift_first_last_with_global_first
    end type
    
    !> Type to encapsulate an MPI domain
    type, public :: mpi_domain
      !> Color: is a number to label this MPI Domain
      integer(i32) :: color
      !> First and last indexes belonging to this MPI Domain
      type(indexes_parallelization) :: index
      !> The MPI environment
      type(mpiinfo) :: mpi_environment
    contains
      procedure, private :: split_from => split_mpi_domain
    end type

    interface define_mpi_domains
      module procedure :: define_mpi_domains_AB
      module procedure :: define_mpi_domains_ABC
    end interface

contains
  pure subroutine set_my_first_last( this, first, last )
    class(indexes_parallelization), intent(inout) :: this
    integer(i32), intent(in) :: first
    integer(i32), intent(in) :: last

    this%my_first = first
    this%my_last = last
  end subroutine

  pure subroutine set_global_first_last( this, first, last )
    class(indexes_parallelization), intent(inout) :: this
    integer(i32), intent(in) :: first
    integer(i32), intent(in) :: last

    this%global_first = first
    this%global_last = last
  end subroutine

  !> Apply a shift to `my_first` and `my_last`
  pure subroutine shift_first_last_by_an_integer( this, shift )
    class(indexes_parallelization), intent(inout) :: this
    integer(i32), intent(in) :: shift

    this%my_first = this%my_first + shift
    this%my_last = this%my_last + shift
  end subroutine

  !> Use `global_first` to shift `my_first` and `my_last`
  pure subroutine shift_first_last_with_global_first( this )
    class(indexes_parallelization), intent(inout) :: this

    call this%shift_first_last_by_an_integer( this%global_first - 1 )
  end subroutine

  !> Split an MPI domain in subdomains
  subroutine split_mpi_domain( this, mpi_type_to_split, n_domains )
    class(mpi_domain), intent(inout) :: this
    !> Type with the MPI Communicator to be split
    type(mpiinfo), intent(in)  :: mpi_type_to_split
    !> Number of MPI (sub)domains to split the current MPI domain
    integer(i32), intent(in) :: n_domains
    
    integer(i32) :: n_tasks, n_effective_domains

    CALL_ASSERT( this%index%global_first <= this%index%global_last, 'First index must be <= than last one' )
    call terminate_if_false(n_domains > 0, 'Number of MPI domains must be positive')
    n_effective_domains = min(n_domains, mpi_type_to_split%procs)
    if (mpi_type_to_split%rank == 0 .and. n_effective_domains < n_domains) then
      call warning('Warning(mod_mpi_gw): Requested ' // trim(to_char(n_domains)) // &
        ' MPI domains, reducing to available rank count: ' // trim(to_char(n_effective_domains)))
    end if
    n_tasks = this%index%global_last - this%index%global_first + 1
    this%color = get_color( mpi_type_to_split%procs, n_effective_domains, mpi_type_to_split%rank )
    call mpi_split( mpi_type_to_split, this%color, this%mpi_environment, .true. )
    call this%index%set_my_first_last( firstofset( this%color, n_tasks, n_effective_domains ), &
                                  lastofset( this%color, n_tasks, n_effective_domains ) )
    call this%index%shift_first_last()
  end subroutine

  !> Get the color of a rank when splitting `n` MPI processes into `n_domains`.   
  !> In MPI, a color is an integer used in functions like `MPI_Comm_split` 
  !> to divide processes into subgroups. Processes with the same color 
  !> are grouped into the same communicator, while those with different 
  !> colors are assigned to separate communicators.
  pure function get_color( n_procs_to_split, n_domains, my_rank ) result(color)
    !> Number of processes to split
    integer(i32), intent(in) :: n_procs_to_split
    !> Number of MPI domains
    integer(i32), intent(in) :: n_domains
    !> Rank of this MPI process
    integer(i32), intent(in) :: my_rank
    !> Color that this rank will have after splitting.
    !> This is a label for the MPI Domain that my_rank will belong to
    integer(i32) :: color

    integer(i32) :: procs_per_group
    
    procs_per_group = n_procs_to_split/n_domains
    if( modulo(n_procs_to_split, n_domains) /= 0 ) procs_per_group = procs_per_group + 1
    color = my_rank/procs_per_group
  end function
  
  !> Split an MPI communicator into groups.
  !> This subroutine is basically a wrapper to `MPI_Comm_split`
  subroutine mpi_split( mpi_type_to_split, color, new_mpi_type, use_rank_as_key )
    !> Type with the MPI Communicator to be split
    type(mpiinfo), intent(in)  :: mpi_type_to_split
    !> Color used to define the MPI domains
    integer(i32), intent(in)   :: color
    !> Type that will contain the new MPI communicator after the split
    type(mpiinfo), intent(out) :: new_mpi_type
    !> If true, use the rank as key to define the order of ranks in the new MPI domains
    logical, intent(in) :: use_rank_as_key

#ifdef MPI
    integer(i32) :: ierror, key
    type(MPI_COMM) :: new_handle
    
    key = merge( mpi_type_to_split%rank, 0, use_rank_as_key )
    call MPI_Comm_split( MPI_COMM(mpi_type_to_split%comm), color, key, new_handle, ierror )
    call new_mpi_type%init( new_handle%mpi_val )
#else
    call new_mpi_type%init( mpi_type_to_split%comm )
#endif
  end subroutine

  !> Define a set of MPI Domains
  subroutine define_mpi_domains_ABC( mpi_global, n_Domains_A, n_Domains_B, mpi_A, mpi_B, C )
    !> Type with the global MPI environment
    type(mpiinfo), intent(in) :: mpi_global
    !> Number of MPI domains to define `mpi_A`
    integer(i32), intent(in) :: n_Domains_A 
    !> Number of MPI domains to define `mpi_B`
    integer(i32), intent(in) :: n_Domains_B
    !> The 1st MPI domain to be defined by splitting mpi_global
    !> On entry, `mpi_A%index%global_first` and `mpi_A%index%global_last` must be already initialized
    type(mpi_domain), intent(inout) :: mpi_A
    !> The 2nd MPI domain to be defined by splitting `mpi_A`
    !> On entry, `mpi_B%index%global_first` and `mpi_B%index%global_last` must be already initialized
    type(mpi_domain), intent(inout) :: mpi_B
    !> Type that encapsulates first and last indexes (of a variable C) to be treated in my rank
    !> On entry, must have `global_first` and `global_last` already initialized
    type(indexes_parallelization), intent(inout) :: C

    integer(i32) :: n_C, i_start, i_end

    call mpi_A%split_from( mpi_global, n_Domains_A )
    call mpi_B%split_from( mpi_A%mpi_environment, n_Domains_B )
    n_C = C%global_last - C%global_first + 1
    call distribute_loop( mpi_B%mpi_environment, n_C, i_start, i_end )
    call C%set_my_first_last( i_start, i_end )
    call C%shift_first_last()

  end subroutine

  !> Define a set of MPI Domains
  subroutine define_mpi_domains_AB( mpi_global, n_Domains_A, mpi_A, B )
    !> Type with the global MPI environment
    type(mpiinfo), intent(in) :: mpi_global
    !> Number of MPI domains to define `mpi_A`
    integer(i32), intent(in) :: n_Domains_A 
    !> The 1st MPI domain to be defined by splitting mpi_global
    !> On entry, `mpi_A%index%global_first` and `mpi_A%index%global_last` must be already initialized
    type(mpi_domain), intent(inout) :: mpi_A
    !> Type that encapsulates first and last indexes (of a variable B) to be treated in my rank
    !> On entry, must have `global_first` and `global_last` already initialized
    type(indexes_parallelization), intent(inout) :: B

    integer(i32) :: n_B, i_start, i_end

    call mpi_A%split_from( mpi_global, n_Domains_A )
    n_B = B%global_last - B%global_first + 1
    call distribute_loop( mpi_A%mpi_environment, n_B, i_start, i_end )
    call B%set_my_first_last( i_start, i_end )
    call B%shift_first_last()

  end subroutine

  !> Pack the first and last indexes into an array of integers
  pure function pack_parallelization_indexes( a ) result(indexes_packed)
    !> Array with the `a`-first and `a`-last indexes
    type(indexes_parallelization), contiguous, intent(in) :: a(:)

    integer(i32) :: i, n
    integer(i32), allocatable :: indexes_packed(:)

    n = size( a )
    allocate( indexes_packed(2*n) )
    do i  = 1, n
      indexes_packed(2*i-1) = a(i)%my_first 
      indexes_packed(2*i) = a(i)%my_last
    end do
  end function

  !> Unpack first and last indexes from an array of integers to an array of 
  !> type(indexes_parallelization). It is the opposite of [[pack_parallelization_indexes]]
  pure function unpack_parallelization_indexes( indexes_packed ) result(a)
    !> Array with all indexes packed
    integer(i32), intent(in) :: indexes_packed(:)
    !> Type with the `a`-first and `a`-last indexes
    type(indexes_parallelization), allocatable :: a(:)

    integer(i32) :: i, n

    n = size( indexes_packed )/2
    allocate( a(n) )
    do i  = 1, n
      call a(i)%set_my_first_last( indexes_packed(2*i-1), indexes_packed(2*i) )
    end do
  end function

  !> Wrapper for MPI_Reduce and MPI_Allreduce
  subroutine mpi_sum_array( array, mpi_env, all_reduce )
    !> The array to which we want to perform a reduction (sum over all MPI ranks)
    complex(dp), contiguous, target, intent(inout) :: array(..)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> Perform MPI_Allreduce, if `all_reduce` is true, or MPI_Reduce, if false.
    logical, intent(in) :: all_reduce

    if ( all_reduce ) then
      call xmpi_allreduce( array, mpi_env )
    else
      call xmpi_reduce( array, mpi_env )
    end if
  end subroutine

end module
