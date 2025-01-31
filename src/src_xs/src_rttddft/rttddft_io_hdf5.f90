!> Module for reading/writing files in HDF5 format
module rttddft_io_hdf5
  use asserts, only: assert
  use math_utils, only: all_close
  use mod_mpi_env, only: mpiinfo
  use modmpi, only: terminate_if_false
  use precision, only: dp, i32
  use rttddft_arrays_utils, only: map_array_to_pointer
  use xhdf5, only: xhdf5_type

  implicit none

  private

  public :: read_array_hdf5, write_array_hdf5

  interface read_array_hdf5
    module procedure :: read_array_rank4
    module procedure :: read_array_rank5
  end interface

  interface write_array_hdf5
    module procedure :: write_array_rank4
    module procedure :: write_array_rank5
  end interface

  character(len=*), parameter :: dims_file_name = 'dimensions'

contains
!> Read an array of rank=4 stored in HDF5 format
subroutine read_array_rank4( h5file, h5path, dataset_name, first, array, descriptors, descriptors_name, mpi_env )
  !> Name of the HDF5 file to read from.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to read from.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> first index along 4th dim (needed to determine offsets)
  integer(i32), intent(in) :: first
  !> array to read from file
  complex(dp), contiguous, intent(out) :: array(:, :, :, first:)
  !> If present, these descriptors must match those stored in the file
  real(dp), contiguous, optional, intent(in) :: descriptors(:, first:)
  !> Descriptors name (used to build HDF5 dataset name for descriptors)
  character(len=*), optional, intent(in) :: descriptors_name
  !> MPI environment
  type(mpiinfo), intent(in):: mpi_env
  
  integer(i32), allocatable :: dims(:, :)
  real(dp), allocatable :: descriptors_in_file(:, :)
  real(dp), parameter :: tol = 1.0e-7_dp
  type(xhdf5_type) :: h5

  call h5%initialize( h5file, mpi_env )
  call h5%read( h5path, dataset_name, array(:, :, :, first:), [1, 1, 1, first] )
  if( present( descriptors ) ) then
    call assert( present( descriptors_name ), 'descriptors_name must be present when descriptors is present')
    call assert( ubound( descriptors, 2 ) == ubound( array, 4 ), "incompatible descriptors and array")
    allocate( descriptors_in_file, mold=descriptors )
    call h5%read( h5path, dataset_name // descriptors_name, descriptors_in_file, [1, first] )
    allocate( dims(3, first:ubound(array, 4)) )
    call h5%read( h5path, dataset_name // dims_file_name, dims, [1, first] )
    call h5%finalize()
    call terminate_if_false( all_close(descriptors_in_file, descriptors, tol), &
      "Descriptors is incongruent with what is stored in file " // h5file // h5path // dataset_name // descriptors_name )
    call terminate_if_false( all( dims == spread([size(array, 1), size(array, 2), size(array, 3)], dim=2, ncopies=size(array, 4) ) ), &
      "Dims is incongruent with what is stored in file " // h5file // h5path // dataset_name // dims_file_name )
  else
    call h5%finalize()
  end if
end subroutine

!> Remap an array of rank 5 to an array of rank 4 using a pointer
!> and then use [[read_array_rank4]]
subroutine read_array_rank5( h5file, h5path, dataset_name, first, array, mpi_env )
  !> Name of the HDF5 file to read from.
  character(len=*), intent(in) :: h5file
  !> Path in the HDF5 file to read from.
  character(len=*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> first index along 5th dim (needed to determine offsets)
  integer(i32), intent(in) :: first
  !> array to be read from file
  complex(dp), contiguous, target, intent(out) :: array(:, :, :, :, first:)
  !> MPI environment
  type(mpiinfo), intent(in) :: mpi_env
  
  complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

  call map_array_to_pointer( first, array, ptr_rank4 )
  call read_array_rank4( h5file, h5path, dataset_name, lbound(ptr_rank4, 4), ptr_rank4, mpi_env=mpi_env )
end subroutine

!> Write an array of rank=4 in HDF5 format
subroutine write_array_rank4( h5file, h5path, dataset_name, array, first, global_size, descriptors, &
    descriptors_name, mpi_env )
  !> Name of the HDF5 file to write in.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to write in.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> first index along 4th dim (needed to determine offsets)
  integer(i32), intent(in) :: first
  !> array to write to file
  complex(dp), contiguous, intent(in) :: array(:, :, :, first:)
  !> global size of `array` along 4th dim.
  integer(i32), intent(in) :: global_size
  !> If present, these descriptors must match those stored in the file
  real(dp), contiguous, optional, intent(in) :: descriptors(:, first:)
  !> Descriptors name (used to build HDF5 dataset name for descriptors)
  character(len=*), optional, intent(in) :: descriptors_name
  !> MPI environment
  type(mpiinfo), intent(in):: mpi_env
  
  integer(i32), allocatable :: dims(:, :)
  type(xhdf5_type) :: h5

  call h5%initialize(h5file, mpi_env)
  associate( m => size(array, 1), n => size(array, 2), k => size(array, 3))
    call h5%write(h5path, dataset_name, array(:, :, :, first:), [1, 1, 1, first], [m, n, k, global_size] )
    if( present( descriptors ) ) then
      call assert( present( descriptors_name ), 'descriptors_name must be present when descriptors is present')
      call h5%write(h5path, dataset_name // descriptors_name, descriptors(:, first:), [1, first], [size(descriptors, 1), global_size])
      dims = spread([m, n, k], dim=2, ncopies=size(array, 4) )
      call h5%write(h5path, dataset_name // dims_file_name, dims, [1, first], [3, global_size] )
    end if
  end associate
  call h5%finalize()
end subroutine

!> Remap an array of rank 5 to an array of rank 4 using a pointer
!> and then use [[write_array_rank4]]
subroutine write_array_rank5( h5file, h5path, dataset_name, array, first, global_size, mpi_env )
  !> Name of the HDF5 file to write in.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to write in.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> first index along 5th dim (needed to determine offsets)
  integer(i32), intent(in) :: first
  !> Array to be written to file
  complex(dp), contiguous, target, intent(in) :: array(:, :, :, :, first:)
  !> global size of `array` along 4th dim.
  integer(i32), intent(in) :: global_size
  !> MPI environment
  type(mpiinfo), intent(in) :: mpi_env
  ! Local variables
  complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

  call map_array_to_pointer( first, array, ptr_rank4 )
  call write_array_rank4( h5file, h5path, dataset_name, ptr_rank4, &
    lbound(ptr_rank4, 4), global_size*size(ptr_rank4, 4), mpi_env=mpi_env )
end subroutine

end module