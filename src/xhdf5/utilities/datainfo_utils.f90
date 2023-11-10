!> Module for utilities to handle information about datasets.
module datainfo_utils
  use iso_fortran_env, only: int32, real32, real64
  use xhdf5_globals
  use xhdf5_error_handling
  use mpi_utils
  use hdf5_datatype


  implicit none 


  private
  public :: datainfo_type


  !> Uninitialized rank
  integer, private, parameter :: uninitialized_rank = -1

  type datainfo_type
    !> Rank of the dataset
    integer(int32) :: rank = uninitialized_rank
    !> Shape of the local data chunk handled by the current MPI rank
    integer(hdf5_size), allocatable :: local_shape(:)
    !> Shape of the full dataset (over all MPI ranks)
    integer(hdf5_size), allocatable :: global_shape(:)
    !> Offset of the local data chunk in context of the whole dataset
    integer(hdf5_ssize), allocatable :: offset(:)

    contains

    procedure :: initialize, sanity_checks, delete
  end type

  contains 

  !> Initialize the dataset properties.
  subroutine initialize(this, local_shape, global_shape, offset, complex_dataset)
    !> Dataset propertis instance
    class(datainfo_type), intent(inout) :: this
    !> Array chunk handled by the current MPI rank.
    integer(int32), intent(in) :: local_shape(:)
    !> Shape of the whole array to be written. If not given, set to [-1].
    integer(int32), intent(in) :: global_shape(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. If not given, set to [-1].
    integer(int32), intent(in) :: offset(:)
    !> Dataset is complex. In this case it is treated as real dataset with an extra rank.
    logical, intent(in) :: complex_dataset

    this%local_shape = local_shape
    this%rank = size(this%local_shape)
    
    ! set rank and shape for scalar input
    if (this%rank == 0) then
      this%local_shape = int([1], hdf5_size)
      this%rank = 1
    end if

    if (sum(global_shape) < 0) then
      this%global_shape = this%local_shape
    else 
      this%global_shape = int(global_shape, hdf5_size)
    end if
    
    if (sum(offset) < 0) then
      this%offset = spread(int(0, hdf5_ssize), 1, this%rank)
    else 
      this%offset = int(offset - 1, hdf5_ssize)
    end if 

    ! complex is treated as real with an extra rank
    if(complex_dataset) then
      this%rank = this%rank + 1
      this%local_shape = [[int(2, hdf5_size)], this%local_shape]
      this%global_shape = [[int(2, hdf5_size)], this%global_shape]
      this%offset = [[int(0, hdf5_ssize)], this%offset]
    end if

  end subroutine initialize


  !> Run a sanity check on the dataset properties instance
  subroutine sanity_checks(this, mpi_comm, calling_routine)
    !> Dataset propertis instance
    class(datainfo_type), intent(inout) :: this
    !> MPI communitcator
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Name of the calling_routine for meaningfull error messages.
    character(*), intent(in) :: calling_routine 

    call xhdf5_assert(mpi_comm, this%rank /= uninitialized_rank, &
                'Error(xhdf5: '// calling_routine //'): datainfo_type instance is not initialized.')

    call xhdf5_assert(mpi_comm, size(this%offset) == this%rank, &
                'Error(xhdf5: '// calling_routine //'): size(this%offset) /= rank.')

    call xhdf5_assert(mpi_comm, size(this%global_shape) == this%rank, &
                'Error(xhdf5: '// calling_routine //'): size(this%global_shape) /= rank.')

    call xhdf5_assert(mpi_comm, all(this%offset >= 0), &
                'Error(xhdf5: '// calling_routine //'): Some elements of this%offset < 0.')

    call xhdf5_assert(mpi_comm, all(this%offset + this%local_shape <= this%global_shape), &
                'Error(xhdf5: '// calling_routine //'): Some elements of this%offset + this%offset > this%global_shape.')
  
  end subroutine sanity_checks

  !> Delete dataset properties
  subroutine delete(this)
    !> Dataset propertis instance
    class(datainfo_type), intent(inout) :: this

    this%rank = uninitialized_rank
    if(allocated(this%local_shape)) deallocate(this%local_shape)
    if(allocated(this%global_shape)) deallocate(this%global_shape)
    if(allocated(this%offset)) deallocate(this%offset)

  end subroutine

end module datainfo_utils