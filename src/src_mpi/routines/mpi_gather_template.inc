! TEMPLATE
#if !defined(TYPE1) || !defined(PRECISION1) || !defined(MPIDATATYPE1) || !defined(RANK_DECL) || !defined(FULL_EXT) || !defined(ALLOC_GATHER) || !defined(ALLOC_GATHERV)
#error "Template macros must be defined before inclusion"
#endif

#define PASTE(X) X
#define PASTE2(X) PASTE(X)_
#define CONCATENATE(X,Y) PASTE2(X)Y
#define GENERATENAME(name) CONCATENATE(name,FULL_EXT)

!> Wrapper for `mpi_gather`. Allocates receive_buffer internally.
subroutine GENERATENAME(mpi_gather)( mpi_env, send_buffer, receive_buffer )
  type(mpiinfo), intent(in) :: mpi_env
  TYPE1(PRECISION1), contiguous, intent(in) :: send_buffer RANK_DECL
  TYPE1(PRECISION1), allocatable, intent(out) :: receive_buffer RANK_DECL
#ifdef MPI
  integer(i32) :: ierr
  integer(long_int) :: send_count
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1

  send_count = size( send_buffer, kind = long_int )

#ifdef USE_MPI4_LARGE_COUNTS
  if ( mpi_env%is_root ) then
    ALLOC_GATHER( send_buffer, receive_buffer, mpi_env%procs )
    call mpi_gather( send_buffer, int( send_count, kind = MPI_COUNT_KIND ), mpi_type, &
      receive_buffer, int( send_count, kind = MPI_COUNT_KIND ), mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  else
    call mpi_gather( send_buffer, int( send_count, kind = MPI_COUNT_KIND ), mpi_type, &
      send_buffer, int( send_count, kind = MPI_COUNT_KIND ), mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  end if
#else
  if ( send_count > int( huge(0_i32), kind = long_int) ) call terminate_mpi_env( mpi_env, &
      "mpi_gather buffer exceeds 32-bit limits. Native large-count support (MPI 4.0+) is required." )

  if ( mpi_env%is_root ) then
    ALLOC_GATHER( send_buffer, receive_buffer, mpi_env%procs )
    call mpi_gather( send_buffer, int( send_count, kind = i32 ), mpi_type, &
      receive_buffer, int( send_count, kind = i32 ), mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  else
    call mpi_gather( send_buffer, int( send_count, kind = i32 ), mpi_type, &
      send_buffer, int( send_count, kind = i32 ), mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  end if
#endif

#else
  receive_buffer = send_buffer
#endif
end subroutine


!> Wrapper for `mpi_gatherv`. Allocates receive_buffer internally.
subroutine GENERATENAME(mpi_gatherv)( mpi_env, send_buffer, receive_buffer )
  type(mpiinfo), intent(in) :: mpi_env
  TYPE1(PRECISION1), contiguous, intent(in) :: send_buffer RANK_DECL
  TYPE1(PRECISION1), allocatable, intent(out) :: receive_buffer RANK_DECL
#ifdef MPI
  integer(long_int) :: send_count, total_elements
  integer(long_int), allocatable :: receive_counts_long(:)
  integer(i32) :: i, ierr
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1

#ifdef USE_MPI4_LARGE_COUNTS
  integer(MPI_COUNT_KIND), allocatable :: receive_counts(:)
  integer(MPI_ADDRESS_KIND), allocatable :: displacements(:)
#else
  integer(i32), allocatable :: receive_counts(:), displacements(:)
  logical :: local_err, global_err
#endif

  send_count = size( send_buffer, kind = long_int )
  call xmpi_get_gatherv_counts( mpi_env, send_count, receive_counts_long )
  total_elements = sum( receive_counts_long )

#ifdef USE_MPI4_LARGE_COUNTS
  if ( mpi_env%is_root ) then
    allocate( receive_counts(mpi_env%procs), displacements(mpi_env%procs) )
    receive_counts = int(receive_counts_long, kind=MPI_COUNT_KIND)
    
    displacements(1) = 0_MPI_ADDRESS_KIND
    do i = 1, mpi_env%procs - 1
      displacements(i + 1) = displacements(i) + int(receive_counts(i), kind=MPI_ADDRESS_KIND)
    end do
    
    ALLOC_GATHERV( send_buffer, receive_buffer, total_elements )
    
    call mpi_gatherv( send_buffer, int( send_count, kind = MPI_COUNT_KIND ), mpi_type, &
      receive_buffer, receive_counts, displacements, mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  else
    allocate( receive_counts(1), displacements(1) ) 
    call mpi_gatherv( send_buffer, int( send_count, kind = MPI_COUNT_KIND ), mpi_type, &
      send_buffer, receive_counts, displacements, mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  end if

#else
  local_err = .false.
  if ( send_count > int( huge(0_i32), kind = long_int) ) local_err = .true.
  if ( mpi_env%is_root ) then
    if ( total_elements > int( huge(0_i32), kind = long_int) ) local_err = .true.
  end if

  call mpi_allreduce( local_err, global_err, 1, MPI_LOGICAL, MPI_LOR, mpi_comm(mpi_env%comm), ierr )
  if ( global_err ) call terminate_mpi_env( mpi_env, &
    "mpi_gatherv buffer exceeds 32-bit limits. Native large-count support (MPI 4.0+) is required." )

  if ( mpi_env%is_root ) then

    allocate( receive_counts(mpi_env%procs), displacements(mpi_env%procs) )
    receive_counts = int(receive_counts_long, kind=i32)
    
    displacements(1) = 0
    do i = 1, mpi_env%procs - 1
      displacements(i + 1) = displacements(i) + receive_counts(i)
    end do
    
    ALLOC_GATHERV(send_buffer, receive_buffer, total_elements)
    
    call mpi_gatherv( send_buffer, int( send_count, kind = i32 ), mpi_type, &
      receive_buffer, receive_counts, displacements, mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  else
    allocate( receive_counts(1), displacements(1) ) 
    call mpi_gatherv( send_buffer, int( send_count, kind = i32 ), mpi_type, &
      send_buffer, receive_counts, displacements, mpi_type, mpi_env%root, mpi_comm(mpi_env%comm), ierr )
  end if
#endif

#else
  receive_buffer = send_buffer
#endif
end subroutine

! Clean macros
#undef PASTE
#undef PASTE2 
#undef CONCATENATE
#undef GENERATENAME
#undef TYPE1 
#undef PRECISION1
#undef MPIDATATYPE1
#undef RANK_DECL
#undef FULL_EXT
#undef ALLOC_GATHER
#undef ALLOC_GATHERV