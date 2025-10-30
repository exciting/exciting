! This "template" implement a generic wrapper to mpi_allgatherv
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
  !> The MPI environment
  type(mpiinfo), intent(inout) :: mpi_env
  !> Buffer. On input, contains the data of the current rank
  !> On output, contains the data of all ranks.
  TYPE1(PRECISION1), contiguous, target, intent(inout) :: buffer(..)
  !> Number of elements handled by the current rank.
  integer(i32), intent(in) :: chunk_size
#ifdef MPI
  integer(i32), allocatable :: receive_counts(:), displacements(:)
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1

  call xmpi_allgather( mpi_env, chunk_size, receive_counts )
  call calculate_displacements( mpi_env, receive_counts, displacements )
  
  call mpi_allgatherv( mpi_in_place, 0, MPI_DATATYPE_NULL, buffer, receive_counts, &
    displacements, mpi_type, mpi_comm( mpi_env%comm ), mpi_env%ierr )
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