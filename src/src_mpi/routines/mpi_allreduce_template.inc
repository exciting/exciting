! This "template" implement a generic wrapper to mpi_allreduce
! TEMPLATE
! TYPE1 : type of the input_vector
! PRECISION1 : precision of TYPE1
! MPIDATATYPE1 : MPI_datatype
#if !defined(TYPE1) || !defined(PRECISION1) || !defined(MPIDATATYPE1)
#error "TYPE1, PRECISION1, and MPIDATATYPE1 must be defined before using this template (mpi_allreduce_template.inc)"
#endif

#define PROCEDURENAME1 mpi_allreduce
#define PASTE(X) X
#define PASTE2(X) PASTE(X)_
#define CONCATENATE(X,Y) PASTE2(X)Y
#define EXTENSION CONCATENATE(TYPE1,PRECISION1)
#define GENERATENAME(name) CONCATENATE(name,EXTENSION) 
#define SMALL small_buffer
#define GENERATENAMESMALL(name) CONCATENATE(GENERATENAME(name),SMALL) 
#define LARGE large_buffer
#define GENERATENAMELARGE(name) CONCATENATE(GENERATENAME(name),LARGE)

!> Wrapper for calling `mpi_allreduce`, assuming an in-place reduce. 
subroutine GENERATENAME(PROCEDURENAME1)( buffer, mpi_env, large_buffer_interface )
  !> Array where the MPI reduction is performed
  TYPE1(PRECISION1), contiguous, target, intent(inout) :: buffer(..)
  !> The MPI environment
  type(mpiinfo), intent(in) :: mpi_env
  !> When present, use the large-buffer interface if the buffer size is larger than a threshold
  logical, optional :: large_buffer_interface
#ifdef MPI
  logical :: big_int_may_be_required
  integer(long_int) :: buffer_size
  TYPE1(PRECISION1), contiguous, pointer :: ptr(:)
  
  big_int_may_be_required = .false.
  if( present(large_buffer_interface) ) big_int_may_be_required = large_buffer_interface
  if ( big_int_may_be_required ) then
    buffer_size = product( [1_long_int*shape(buffer)] )
    if( buffer_size < big_int ) then
      call GENERATENAMESMALL(PROCEDURENAME1)( buffer, mpi_env )
    else
      call c_f_pointer(c_loc(buffer), ptr, [buffer_size])
      call GENERATENAMELARGE(PROCEDURENAME1)( ptr, buffer_size, mpi_env )
    end if
  else
    call GENERATENAMESMALL(PROCEDURENAME1)( buffer, mpi_env )
  end if
#endif
end subroutine

!> (private) Wrapper for calling `mpi_allreduce`, assuming an in-place reduce. 
!> Case: small-sized buffers
subroutine GENERATENAMESMALL(PROCEDURENAME1)( buffer, mpi_env )
  !> Array where the MPI reduction is performed
  TYPE1(PRECISION1), contiguous, intent(inout) :: buffer(..)
  type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
  integer(i32) :: i_err
  integer(i32), parameter :: receive_buffer_ignore = 0
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1
  
  call mpi_allreduce( MPI_IN_PLACE, buffer, size( buffer ), mpi_type, MPI_SUM, MPI_COMM(mpi_env%comm), i_err )
#endif
end subroutine

!> (private) Wrapper for calling `mpi_allreduce`, assuming an in-place reduce. 
!> Case: large-sized buffers
subroutine GENERATENAMELARGE(PROCEDURENAME1)( buffer, buffer_size, mpi_env )
  !> The buffer to which we want to perform a reduction (sum over all MPI ranks)
  TYPE1(PRECISION1), intent(inout) :: buffer(*)
  !> Size of `buffer` (expected not to fit into `i32` integer)
  integer(long_int), intent(in) :: buffer_size
  type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
  integer(i32), parameter :: chunk_size = big_int
  integer(long_int) :: first, last, i_chunck, n_chunks
  n_chunks = buffer_size/chunk_size
  first = 1
  do i_chunck = 1, n_chunks-1
    last = max(first + chunk_size - 1, buffer_size)
    call GENERATENAMESMALL(PROCEDURENAME1)( buffer(first:last), mpi_env )
    first = last + 1
  end do
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