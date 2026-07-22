!> Module for reading/writing files in HDF5 format
module rttddft_io_hdf5
#include "asserts.fpp"
  use math_utils, only: all_close
  use mod_mpi_env, only: mpiinfo
  use modmpi, only: terminate_if_false
  use os_utils, only: join_paths, path_exists
  use precision, only: dp, i32
  use rttddft_arrays_utils, only: map_array_to_pointer
  use xhdf5, only: xhdf5_type

  implicit none

  private

  public :: dataset_exists, read_array_hdf5, read_three_arrays_hdf5, write_array_hdf5, write_three_arrays_hdf5

  interface read_array_hdf5
    module procedure :: read_array_real_dp
    module procedure :: read_array_complex_dp
    module procedure :: read_array_rank4
    module procedure :: read_array_rank5
  end interface

  interface read_three_arrays_hdf5
    module procedure :: read_three_arrays_rank3
  end interface

  interface write_array_hdf5
    module procedure :: write_array_complex_dp
    module procedure :: write_array_real_dp
    module procedure :: write_array_rank4
    module procedure :: write_array_rank5
  end interface

  interface write_three_arrays_hdf5
    module procedure :: write_three_arrays_rank3
  end interface

  character(len=*), parameter :: dims_file_name = 'dimensions'

contains
!> Return `.true.` if the dataset exists
logical function dataset_exists( h5file, h5path, dataset_name, mpi_env )
  !> Name of the HDF5 file to read from.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to read from.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> MPI environment
  type(mpiinfo), intent(in):: mpi_env

  type(xhdf5_type) :: h5

  dataset_exists = path_exists( h5file )
  if( dataset_exists ) then
    call h5%initialize( h5file, mpi_env )
    dataset_exists = h5%exists( join_paths( h5path, dataset_name ) )
    call h5%finalize( )
  end if
end function

!> Read a real-dp array. All MPI ranks will read the entire array
subroutine read_array_real_dp( h5file, h5path, dataset_name, array, mpi_env )
  !> Name of the HDF5 file to read from
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to read from
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> Array to read from file
  real(dp), contiguous, intent(out) :: array(..)
  !> MPI environment
  type(mpiinfo), intent(in) :: mpi_env

  type(xhdf5_type) :: h5

  call h5%initialize( h5file, mpi_env )
  call h5%read( h5path, trim(dataset_name), array )
  call h5%finalize( )
end subroutine

!> Read a complex-dp array. All MPI ranks read the entire array
subroutine read_array_complex_dp( h5file, h5path, dataset_name, array, mpi_env )
  !> Name of the HDF5 file to read from
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to read from
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> Array to read from file
  complex(dp), contiguous, intent(out) :: array(..)
  !> MPI environment
  type(mpiinfo), intent(in) :: mpi_env

  type(xhdf5_type) :: h5

  call h5%initialize( h5file, mpi_env )
  call h5%read( h5path, trim(dataset_name), array )
  call h5%finalize( )
end subroutine

!> Read an array of rank=4, with each MPI rank reading its correspondent chunk. 
!> The chunk is determined by the offset `[1, 1, 1, first]`, where `first` is rank-dependent
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
  call h5%read( h5path, trim(dataset_name), array(:, :, :, first:), [1, 1, 1, first] )
  if( present( descriptors ) ) then
    CALL_ASSERT( present( descriptors_name ), 'descriptors_name must be present when descriptors is present')
    CALL_ASSERT( ubound( descriptors, 2 ) == ubound( array, 4 ), "incompatible descriptors and array")
    allocate( descriptors_in_file, mold=descriptors )
    call h5%read( h5path, trim(dataset_name) // trim(descriptors_name), descriptors_in_file, [1, first] )
    allocate( dims(3, first:ubound(array, 4)) )
    call h5%read( h5path, trim(dataset_name) // trim(dims_file_name), dims, [1, first] )
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
  call read_array_rank4( h5file, h5path, trim(dataset_name), lbound(ptr_rank4, 4), ptr_rank4, mpi_env=mpi_env )
end subroutine

!> Read 3 arrays of rank=3, with each MPI rank reading its correspondent arrays chuncks. 
subroutine read_three_arrays_rank3( h5file, h5path, dataset_name, first, array_1, array_2, array_3, mpi_env )
  !> Name of the HDF5 file to read from.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to read from.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> first index along 3rd dim (needed to determine offsets)
  integer(i32), intent(in) :: first
  !> 1st array to read from file
  complex(dp), contiguous, intent(out) :: array_1(:, :, first:)
  !> 2nd array to read from file
  complex(dp), contiguous, intent(out) :: array_2(:, :, first:)
  !> 3rd array to read from file
  complex(dp), contiguous, intent(out) :: array_3(:, :, first:)
  !> MPI environment
  type(mpiinfo), intent(in):: mpi_env
  
  type(xhdf5_type) :: h5

  CALL_ASSERT( all( shape(array_1) == shape(array_2) ), "shape mismatch - array 1 and 2" )
  CALL_ASSERT( all( shape(array_1) == shape(array_3) ), "shape mismatch - array 1 and 3" )
  call h5%initialize( h5file, mpi_env )
  call h5%read( trim(h5path), trim(dataset_name)//'-1', array_1(:, :, first:), [1, 1, first] )
  call h5%read( trim(h5path), trim(dataset_name)//'-2', array_2(:, :, first:), [1, 1, first] )
  call h5%read( trim(h5path), trim(dataset_name)//'-3', array_3(:, :, first:), [1, 1, first] )
  call h5%finalize()
end subroutine

!> Write a complex-dp array (currently only in serial mode, see issue #231)
subroutine write_array_complex_dp( h5file, h5path, dataset_name, array, mpi_env, serial_access )
  !> Name of the HDF5 file to write in.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to write in.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> Array to be written to file
  complex(dp), contiguous, intent(in) :: array(..)
  !> MPI environment
  type(mpiinfo), intent(in) :: mpi_env
  !> If `.true.`, write in serial mode within the current MPI environment
  logical, optional, intent(in) :: serial_access

  logical :: serial_access_local
  type(xhdf5_type) :: h5

  serial_access_local = .false.
  if( present(serial_access) ) serial_access_local = serial_access
  CALL_ASSERT( serial_access_local, "write_array_complex_dp currently only works in serial mode" )
  if( mpi_env%is_root ) then
    call h5%initialize( h5file, mpi_env, serial_access_local )
    call h5%write( h5path, trim( dataset_name ), array )
    call h5%finalize( )
  end if
end subroutine

!> Write a real-dp array (currently only in serial mode, see issue #231)
subroutine write_array_real_dp( h5file, h5path, dataset_name, array, mpi_env, serial_access )
  !> Name of the HDF5 file to write in
  character(len=*), intent(in) :: h5file
  !> Path in the HDF5 file to write in
  character(len=*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> Array to be written to file
  real(dp), contiguous, intent(in) :: array(..)
  !> MPI environment
  type(mpiinfo), intent(in) :: mpi_env
  !> If `.true.`, write in serial mode within the current MPI environment
  logical, optional, intent(in) :: serial_access

  logical :: serial_access_local
  type(xhdf5_type) :: h5

  serial_access_local = .true.
  if( present(serial_access) ) serial_access_local = serial_access
  CALL_ASSERT( serial_access_local, "write_array_real_dp currently only works in serial mode" )
  if( mpi_env%is_root ) then
    call h5%initialize( h5file, mpi_env, serial_access_local )
    call h5%write( h5path, trim( dataset_name ), array )
    call h5%finalize( )
  end if
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
    call h5%write(h5path, trim(dataset_name), array(:, :, :, first:), [1, 1, 1, first], [m, n, k, global_size] )
    if( present( descriptors ) ) then
      CALL_ASSERT( present( descriptors_name ), 'descriptors_name must be present when descriptors is present')
      call h5%write(h5path, trim(dataset_name) // trim(descriptors_name), descriptors(:, first:), [1, first], [size(descriptors, 1), global_size])
      dims = spread([m, n, k], dim=2, ncopies=size(array, 4) )
      call h5%write(h5path, trim(dataset_name) // trim(dims_file_name), dims, [1, first], [3, global_size] )
    end if
  end associate
  call h5%finalize( )
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
  
  complex(dp), contiguous, pointer :: ptr_rank4(:, :, :, :)

  call map_array_to_pointer( first, array, ptr_rank4 )
  call write_array_rank4( h5file, h5path, trim(dataset_name), ptr_rank4, &
    lbound(ptr_rank4, 4), global_size*size(ptr_rank4, 4), mpi_env=mpi_env )
end subroutine

!> Write 3 arrays of rank=3, with each MPI rank reading its correspondent arrays chuncks. 
subroutine write_three_arrays_rank3( h5file, h5path, dataset_name, array_1, array_2, array_3, first, global_size, mpi_env )
  !> Name of the HDF5 file to read from.
  character(*), intent(in) :: h5file
  !> Path in the HDF5 file to read from.
  character(*), intent(in) :: h5path
  !> Dataset name
  character(len=*), intent(in) :: dataset_name
  !> first index along 3rd dim (needed to determine offsets)
  integer(i32), intent(in) :: first
  !> global size of `array` along 4th dim.
  integer(i32), intent(in) :: global_size
  !> 1st array to read from file
  complex(dp), contiguous, intent(in) :: array_1(:, :, first:)
  !> 2nd array to read from file
  complex(dp), contiguous, intent(in) :: array_2(:, :, first:)
  !> 3rd array to read from file
  complex(dp), contiguous, intent(in) :: array_3(:, :, first:)
  !> MPI environment
  type(mpiinfo), intent(in):: mpi_env
  
  type(xhdf5_type) :: h5
  integer(i32) :: m, n

  m = size(array_1, 1); n = size(array_1, 2)
  CALL_ASSERT( all( shape(array_1) == shape(array_2) ), "shape mismatch - array 1 and 2" )
  CALL_ASSERT( all( shape(array_1) == shape(array_3) ), "shape mismatch - array 1 and 3" )
  call h5%initialize( h5file, mpi_env )
  call h5%write( trim(h5path), trim(dataset_name)//'-1', array_1(:, :, first:), [1, 1, first], [m, n, global_size] )
  call h5%write( trim(h5path), trim(dataset_name)//'-2', array_2(:, :, first:), [1, 1, first], [m, n, global_size] )
  call h5%write( trim(h5path), trim(dataset_name)//'-3', array_3(:, :, first:), [1, 1, first], [m, n, global_size] )
  call h5%finalize()
end subroutine

end module