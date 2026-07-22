! This "template" implements a generic wrapper to mpi_bcast
! TEMPLATE
! TYPE1 : type of the input_vector
! PRECISION1 : precision of TYPE1
! MPIDATATYPE1 : MPI_datatype
#if !defined(TYPE1) || !defined(PRECISION1) || !defined(MPIDATATYPE1)
#error "TYPE1, PRECISION1, and MPIDATATYPE1 must be defined before using this template"
#endif

#define PROCEDURENAME1 mpi_bcast
#define PASTE(X) X
#define PASTE2(X) PASTE(X)_
#define CONCATENATE(X,Y) PASTE2(X)Y
#define EXTENSION CONCATENATE(TYPE1,PRECISION1)
#define GENERATENAME(name) CONCATENATE(name,EXTENSION) 

!> Wrapper for calling `mpi_bcast` automatically handling arbitrary ranks and huge buffers
subroutine GENERATENAME(PROCEDURENAME1)( mpi_env, buffer, bcasting_rank )
  type(mpiinfo), intent(in) :: mpi_env
  TYPE1(PRECISION1), contiguous, target, intent(inout) :: buffer(..)
  integer(i32), optional, intent(in) :: bcasting_rank
#ifdef MPI
  integer(long_int) :: buffer_size
  integer(i32) :: i_err, root
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1
  
  ! Chunking variables for fallback
  integer(long_int) :: first, last, chunk_elements
  TYPE1(PRECISION1), pointer :: buf_1d(:)
  
  root = 0
  if ( present(bcasting_rank) ) root = bcasting_rank
  
  buffer_size = size( buffer, kind = long_int )
  
  if ( buffer_size > big_int ) then
#ifdef USE_MPI4_LARGE_COUNTS
    call mpi_bcast( buffer, int( buffer_size, kind = MPI_COUNT_KIND), &
                    mpi_type, root, mpi_comm(mpi_env%comm), i_err )
#else
    ! Fallback: Map assumed-rank to 1D array and chunk
    call c_f_pointer( c_loc(buffer), buf_1d, [buffer_size] )
    
    first = 1_long_int
    do while (first <= buffer_size)
      last = min( first + big_int - 1_long_int, buffer_size )
      chunk_elements = last - first + 1_long_int
      
      call mpi_bcast( buf_1d(first:last), int( chunk_elements, kind = i32 ), &
        mpi_type, root, mpi_comm(mpi_env%comm), i_err )
      first = last + 1_long_int
    end do
#endif
  else
    call mpi_bcast( buffer, int( buffer_size, kind = i32 ), &
      mpi_type, root, mpi_comm(mpi_env%comm), i_err )
  end if
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