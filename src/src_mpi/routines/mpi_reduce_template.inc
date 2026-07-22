! This "template" implements a generic wrapper to mpi_reduce
! TEMPLATE
! TYPE1 : type of the input_vector
! PRECISION1 : precision of TYPE1
! MPIDATATYPE1 : MPI_datatype
#if !defined(TYPE1) || !defined(PRECISION1) || !defined(MPIDATATYPE1)
#error "TYPE1, PRECISION1, and MPIDATATYPE1 must be defined before using this template (mpi_reduce_template.inc)"
#endif

#define PROCEDURENAME1 mpi_reduce
#define PASTE(X) X
#define PASTE2(X) PASTE(X)_
#define CONCATENATE(X,Y) PASTE2(X)Y
#define EXTENSION CONCATENATE(TYPE1,PRECISION1)
#define GENERATENAME(name) CONCATENATE(name,EXTENSION) 

!> Wrapper for calling `mpi_reduce`, assuming an in-place reduce on the root. 
subroutine GENERATENAME(PROCEDURENAME1)( buffer, mpi_env )
  !> Array where the MPI reduction is performed
  TYPE1(PRECISION1), contiguous, target, intent(inout) :: buffer(..)
  !> The MPI environment
  type(mpiinfo), intent(in) :: mpi_env
#ifdef MPI
  integer(long_int) :: buffer_size
  integer(i32) :: i_err
  type(MPI_Datatype), parameter :: mpi_type = MPIDATATYPE1
  
  ! Use a local array instead of a parameter to avoid strict intent(out) compiler errors
  integer(i32) :: receive_buffer_ignore(1) 
  
  ! Chunking variables for fallback
  integer(long_int) :: first, last, chunk_elements
  TYPE1(PRECISION1), pointer :: buf_1d(:)
  
  buffer_size = size( buffer, kind=long_int )
  
  if ( buffer_size > big_int ) then
#ifdef USE_MPI4_LARGE_COUNTS
    if ( mpi_env%is_root ) then
      call MPI_Reduce( MPI_IN_PLACE, buffer, int( buffer_size, kind = MPI_COUNT_KIND ), &
        mpi_type, MPI_SUM, mpi_env%root, mpi_comm(mpi_env%comm), i_err )
    else
      call MPI_Reduce( buffer, receive_buffer_ignore, int( buffer_size, kind = MPI_COUNT_KIND ), &
        mpi_type, MPI_SUM, mpi_env%root, mpi_comm(mpi_env%comm), i_err )
    end if
#else
    ! Fallback: Map assumed-rank to 1D array and chunk
    call c_f_pointer( c_loc(buffer), buf_1d, [buffer_size] )
    
    first = 1_long_int
    do while ( first <= buffer_size )
      last = min( first + big_int - 1_long_int, buffer_size )
      chunk_elements = last - first + 1_long_int
      
      if ( mpi_env%is_root ) then
        call MPI_Reduce( MPI_IN_PLACE, buf_1d(first:last), int(chunk_elements, kind=i32), &
          mpi_type, MPI_SUM, mpi_env%root, mpi_comm(mpi_env%comm), i_err )
      else
        call MPI_Reduce( buf_1d(first:last), receive_buffer_ignore, int(chunk_elements, kind=i32), &
          mpi_type, MPI_SUM, mpi_env%root, mpi_comm(mpi_env%comm), i_err )
      end if
      
      first = last + 1_long_int
    end do
#endif
  else
    if ( mpi_env%is_root ) then
      call MPI_Reduce( MPI_IN_PLACE, buffer, int( buffer_size, kind = i32 ), &
        mpi_type, MPI_SUM, mpi_env%root, mpi_comm(mpi_env%comm), i_err )
    else
      call MPI_Reduce( buffer, receive_buffer_ignore, int( buffer_size, kind = i32 ), &
        mpi_type, MPI_SUM, mpi_env%root, mpi_comm(mpi_env%comm), i_err )
    end if
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