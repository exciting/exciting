! This "template" implements a generic wrapper to mpi_allgatherv
! TEMPLATE
! TYPE1 : type of the input_vector
! PRECISION1 : precision of TYPE1
! MPIDATATYPE1 : MPI_datatype
#if !defined(TYPE1) || !defined(PRECISION1) || !defined(MPIDATATYPE1)
#error "TYPE1, PRECISION1, and MPIDATATYPE1 must be defined before using this template (mpi_allgatherv_template.inc)"
#endif

#define PROCEDURENAME1 mpi_allgatherv_in_place
#define PASTE(X) X
#define PASTE2(X) PASTE(X)_
#define CONCATENATE(X,Y) PASTE2(X)Y
#define EXTENSION CONCATENATE(TYPE1,PRECISION1)
#define GENERATENAME(name) CONCATENATE(name,EXTENSION)


!> Wrapper for calling `mpi_allgatherv`, assuming an in-place gather. 
subroutine GENERATENAME(PROCEDURENAME1)( mpi_env, buffer, chunk_size )
  type(mpiinfo), intent(in) :: mpi_env
  TYPE1(PRECISION1), contiguous, target, intent(inout) :: buffer(..)
  integer(i32), intent(in) :: chunk_size
#ifdef MPI
  integer(i32) :: i, ierr
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1

#ifdef USE_MPI4_LARGE_COUNTS
  ! ====================================================================
  ! PATH 1: Native MPI-4 Large Counts (Stable MPI Libraries)
  ! ====================================================================
  integer(i32), allocatable :: local_receive_counts(:)
  integer(MPI_COUNT_KIND), allocatable :: receive_counts(:)
  integer(MPI_ADDRESS_KIND), allocatable :: displacements(:)

  allocate( receive_counts(mpi_env%procs) )
  allocate( displacements(mpi_env%procs) )

  call xmpi_allgather( mpi_env, chunk_size, local_receive_counts )
  receive_counts = int( local_receive_counts, kind = MPI_COUNT_KIND )

  displacements(1) = 0_MPI_ADDRESS_KIND
  do i = 1, mpi_env%procs - 1
    displacements(i + 1) = displacements(i) + receive_counts(i)
  end do
  
  call mpi_allgatherv( MPI_IN_PLACE, 0_MPI_COUNT_KIND, MPI_DATATYPE_NULL, buffer, receive_counts, &
    displacements, mpi_type, mpi_comm( mpi_env%comm ), ierr )

#else
  ! ====================================================================
  ! PATH 2: Standard 32-bit Integer (Older MPI / Fallback)
  ! ====================================================================
  integer(i32), allocatable :: receive_counts(:)
  integer(i32), allocatable :: displacements(:)
  integer(long_int) :: total_elements

  allocate( displacements(mpi_env%procs) )

  call xmpi_allgather( mpi_env, chunk_size, receive_counts )

  total_elements = sum( int( receive_counts, kind = long_int ) )
  if ( total_elements > int( huge( 0_i32 ), kind = long_int) ) &
    call terminate_mpi_env( mpi_env, "mpi_allgatherv buffer exceeds 32-bit limits. &
    Native large-count support (MPI 4.0+) is required.")

  displacements(1) = 0
  do i = 1, mpi_env%procs - 1
    displacements(i + 1) = displacements(i) + receive_counts(i)
  end do
  
  call mpi_allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buffer, receive_counts, &
    displacements, mpi_type, mpi_comm( mpi_env%comm ), ierr )
#endif

#endif
end subroutine


!> Wrapper for calling `mpi_allgatherv`, assuming an in-place gather and large buffer size.
subroutine GENERATENAME(PROCEDURENAME1)_large_buffer( mpi_env, buffer, chunk_size )
  type(mpiinfo), intent(in) :: mpi_env
  TYPE1(PRECISION1), contiguous, target, intent(inout) :: buffer(..)
  integer(long_int), intent(in) :: chunk_size
#ifdef MPI
  integer(i32) :: i, ierr
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1

#ifdef USE_MPI4_LARGE_COUNTS
  ! ====================================================================
  ! PATH 1: Native MPI-4 Large Counts (Stable MPI Libraries)
  ! ====================================================================
  integer(long_int), allocatable :: local_receive_counts(:)
  integer(MPI_COUNT_KIND), allocatable :: receive_counts(:)
  integer(MPI_ADDRESS_KIND), allocatable :: displacements(:)

  allocate( receive_counts(mpi_env%procs) )
  allocate( displacements(mpi_env%procs) )

  call xmpi_allgather( mpi_env, chunk_size, local_receive_counts )
  receive_counts = int( local_receive_counts, kind = MPI_COUNT_KIND )

  displacements(1) = 0_MPI_ADDRESS_KIND
  do i = 1, mpi_env%procs - 1
    displacements(i + 1) = displacements(i) + receive_counts(i)
  end do
  
  call mpi_allgatherv( MPI_IN_PLACE, 0_MPI_COUNT_KIND, MPI_DATATYPE_NULL, buffer, receive_counts, &
    displacements, mpi_type, mpi_comm( mpi_env%comm ), ierr )

#else
  ! ====================================================================
  ! PATH 2: Standard 32-bit Integer (Older MPI / Fallback)
  ! ====================================================================
  integer(long_int), allocatable :: local_receive_counts(:)
  integer(i32), allocatable :: receive_counts(:)
  integer(i32), allocatable :: displacements(:)
  integer(long_int) :: total_elements

  allocate( receive_counts(mpi_env%procs) )
  allocate( displacements(mpi_env%procs) )

  call xmpi_allgather( mpi_env, chunk_size, local_receive_counts )
  
  total_elements = sum( local_receive_counts )
  if ( total_elements > int( huge( 0_i32 ), kind = long_int) ) &
    call terminate_mpi_env( mpi_env, "mpi_allgatherv buffer exceeds 32-bit limits. &
    Native large-count support (MPI 4.0+) is required.")

  receive_counts = int( local_receive_counts, kind = i32 )

  displacements(1) = 0
  do i = 1, mpi_env%procs - 1
    displacements(i + 1) = displacements(i) + receive_counts(i)
  end do
  
  call mpi_allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, buffer, receive_counts, &
    displacements, mpi_type, mpi_comm( mpi_env%comm ), ierr )
#endif

#endif
end subroutine

! Clean functions
#undef PASTE
#undef PASTE2 
#undef CONCATENATE
#undef EXTENSION
#undef GENERATENAME

! Clean the name of the procedure
#undef PROCEDURENAME1

! Clean templated type and precision
#undef TYPE1 
#undef PRECISION1
#undef MPIDATATYPE1
