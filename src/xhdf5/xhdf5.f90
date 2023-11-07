!> Module for HDF5 wrappers.
!> Easy to use fortran wrappers for the HDF5 library. 
!> Parallel MPI writing of array strides is supported.
!> Use these routines for writing large binary files that can be easily reused on other mashines.
!>
!> The use of this module requires to build exciting with HDF5. 
!> If exciting is compiled without HDF5 the routines will do nothing. If your feature requires HDF5,
!> Make sure to call `[[abort_if_not_hdf5]]` in the beginning of the code to make sure that exciting gets
!> safely terminated if built without HDF5.
!>
!> As already mentioned, HDF5 is a library to write large binary files, especially writing large arrays. 
!> Therefore it does absolutely make no sense to use them in loops.
!> If you are tempted to do this, you should seriously rethink about your data layout!
!>
!> For examples how to use the library we refer kindly to the unit tests (`[[xhdf5_test.f90]]).
!> These show all examples of usage.
!>
!> Be aware of that the writing routines will overwrite an existing datasets with the same name.
module xhdf5
  use iso_c_binding, only: c_ptr, c_loc

  use precision, only: dp
  use os_utils, only: path_exists, join_paths
  use asserts, only: assert

  use hdf5_utils
  use hdf5_file
  use hdf5_group
  use hdf5_write
  use hdf5_read
  use hdf5_dataset_utils


  implicit none


  private
  public :: xhdf5_type, abort_if_not_hdf5


  !> Value for undefined integers.
  integer(hdf5_id), parameter :: h5id_undefined = -1


  !> Type for handeling an HDF5 file.
  type xhdf5_type
    !> OS file path.
    character(:), allocatable :: path
    !> HDF5 file id.
    integer(hdf5_id) :: h5id = h5id_undefined
    !> MPI communicator handle
    integer :: mpi_comm
    !> Flag to initialize HDF5 file in serial mode, even if compiled with MPI.
    !> By default this is set to false.
    logical :: serial_access = .false.

    contains

    procedure :: initialize
    procedure :: finalize
    procedure :: initialize_group
    procedure :: dataset_shape

    procedure :: exists => link_exists

    procedure :: delete => delete_link

    generic :: write => write_string, &
                        write_integer_rank_0, write_integer_rank_1, write_integer_rank_2, &
                        write_real_dp_rank_0, write_real_dp_rank_1, write_real_dp_rank_2, write_real_dp_rank_3, write_real_dp_rank_4, &
                        write_real_sp_rank_0, write_real_sp_rank_1, write_real_sp_rank_2, write_real_sp_rank_3, write_real_sp_rank_4, &
                        write_complex_dp_rank_1, write_complex_dp_rank_2, write_complex_dp_rank_3, write_complex_dp_rank_4

    procedure :: write_string, &
                 write_integer_rank_0, write_integer_rank_1, write_integer_rank_2, &
                 write_real_dp_rank_0, write_real_dp_rank_1, write_real_dp_rank_2, write_real_dp_rank_3, write_real_dp_rank_4, &
                 write_real_sp_rank_0, write_real_sp_rank_1, write_real_sp_rank_2, write_real_sp_rank_3, write_real_sp_rank_4, &
                 write_complex_dp_rank_1, write_complex_dp_rank_2, write_complex_dp_rank_3, write_complex_dp_rank_4

    generic :: read => read_string, &
                       read_integer_rank_0, read_integer_rank_1, read_integer_rank_2, &
                       read_real_dp_rank_0, read_real_dp_rank_1, read_real_dp_rank_2, read_real_dp_rank_3, read_real_dp_rank_4, &
                       read_real_sp_rank_0, read_real_sp_rank_1, read_real_sp_rank_2, read_real_sp_rank_3, read_real_sp_rank_4, &
                       read_complex_dp_rank_1, read_complex_dp_rank_2, read_complex_dp_rank_3, read_complex_dp_rank_4

    procedure :: read_string, &
                 read_integer_rank_0, read_integer_rank_1, read_integer_rank_2, &
                 read_real_dp_rank_0, read_real_dp_rank_1, read_real_dp_rank_2, read_real_dp_rank_3, read_real_dp_rank_4, &
                 read_real_sp_rank_0, read_real_sp_rank_1, read_real_sp_rank_2, read_real_sp_rank_3, read_real_sp_rank_4, &
                 read_complex_dp_rank_1, read_complex_dp_rank_2, read_complex_dp_rank_3, read_complex_dp_rank_4            

  end type xhdf5_type


  contains


  !> Initialize HDF5 library and, if `path` exists, open the file, else create a file at `path`.
  subroutine initialize(this, path, mpi_comm, serial_access)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Relative path to HDF5 file.
    character(*), intent(in) :: path
    !> MPI communicator handle.
    integer, intent(in) :: mpi_comm
    !> Set to `.true.` if only serial access is possible.
    logical, intent(in), optional :: serial_access

    this%path = trim(path)
    this%mpi_comm = mpi_comm

    if(present(serial_access)) then
      this%serial_access = serial_access
    end if

    call hdf5_initialize()

    if(path_exists(this%path)) then
      call hdf5_open_file(this%path, this%mpi_comm, this%h5id, this%serial_access)
    else 
      call hdf5_create_file(this%path, this%mpi_comm, this%h5id, this%serial_access)
    end if    
  end subroutine


  !> Close file and finalize HDF5 library.
  subroutine finalize(this)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this

    if (this%h5id == h5id_undefined) return
    call hdf5_close_file(this%h5id)
    call hdf5_finalize()
    this%h5id = h5id_undefined
  end subroutine finalize 
  

  !> Create a new group with name `group` at `h5path`. If the group already exists, the routine does nothing.
  subroutine initialize_group(this, h5path, groupname)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path
    !> Group name.
    character(*), intent(in) :: groupname

    if (this%exists(join_paths(h5path, groupname))) return
    call hdf5_create_group(this%h5id, trim(h5path), groupname)
  end subroutine initialize_group
  

  !> Return true or false, whether the link to `h5path` exists or not.
  logical function link_exists(this, h5path)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path

    link_exists = hdf5_link_exists(this%h5id, trim(h5path))
  end function link_exists

  !> Delete link to `h5path`. This action frees the space on the drive allocated by the object `h5path` points to.
  subroutine delete_link(this, h5path)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path

    call hdf5_delete_link(this%h5id, trim(h5path))
  end subroutine delete_link


  !> Get the shape of a dataset.
  subroutine dataset_shape(this, h5path, datasetname, shape, complex_dataset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path
    !> Dataset name.
    character(*), intent(in) :: datasetname
    !> Dataset shape
    integer, allocatable :: shape(:)
    !> Set to `.true.` if the dataset is complex. If not, the first dimension will be 
    !> always 2 for complex datasets.
    logical, optional :: complex_dataset

    logical :: complex_dataset_local
    integer(hdf5_size), allocatable :: shape_local(:)

    complex_dataset_local = .false.
    if(present(complex_dataset)) complex_dataset_local = complex_dataset

    call hdf5_get_dataset_shape(this%h5id, h5path, datasetname, shape_local)
    
    if(complex_dataset_local) then 
      shape = shape_local(2:)
    else 
      shape = shape_local
    end if
  end subroutine dataset_shape


!-----------------------------------------------------------------------------------------------------------------------
! WRITE DATASET


! CHARACTER


  !> Write a string to an HDF5 file.
  subroutine write_string(this, h5path, dataset, string)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> String to write. Must be shorter then `[[hdf5_max_string_len]]`.
    character(*), intent(in), target :: string

    integer, parameter :: dataset_rank = 1


    character(hdf5_max_string_len), target :: string_local
    character(:), allocatable :: string_trim
    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    string_local = repeat(hdf5_blank, hdf5_max_string_len)
    string_trim = trim(adjustl(string))
    
    call assert(len(string_trim) <= hdf5_max_string_len, 'len(string) > hdf5_max_string_len.')

    string_local(1 : len(string_trim)) = string_trim

    data_type_id = hdf5_character()
    dataset_chunk_ptr = c_loc(string_local)
    dataset_chunk_shape = [hdf5_max_string_len]
    dataset_shape_local = [hdf5_max_string_len]
    offset_local = [0]

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_string


! INTEGER

  
  !> Write an integer scalar to an HDF5 file.
  subroutine write_integer_rank_0(this, h5path, dataset, scalar)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Integer to write.
    integer, intent(in), target :: scalar

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_integer()
    dataset_chunk_ptr = c_loc(scalar)
    dataset_chunk_shape = [1]
    dataset_shape_local = [1]
    offset_local = [0]

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_integer_rank_0


  !> Write an integer vector to an HDF5 file.
  subroutine write_integer_rank_1(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Chunk of the vector handled by the current MPI rank.
    integer, intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array.  
    integer, intent(in) :: offset(1)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(1)

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(1), dataset_shape_local(1)
    integer(hdf5_ssize) :: offset_local(1)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_integer()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_integer_rank_1


  !> Write an integer matrix to an HDF5 file.
  subroutine write_integer_rank_2(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    integer, intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(2)

    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_integer()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local,        &
                            this%serial_access)
  end subroutine write_integer_rank_2


! REAL(SP)


  !> Write a single precision scalar to an HDF5 file.
  subroutine write_real_sp_rank_0(this, h5path, dataset, scalar)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Real(sp) to write.
    real(sp), intent(in), target :: scalar

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(scalar)
    dataset_chunk_shape = [1]
    dataset_shape_local = [1]
    offset_local = [0] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_sp_rank_0


  !> Write a single precision vector to an HDF5 file.
  subroutine write_real_sp_rank_1(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(1)

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_sp_rank_1


  !> Write a single precision two rank array to an HDF5 file.
  subroutine write_real_sp_rank_2(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(2)

    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_sp_rank_2


  !> Write a single precision three rank array to an HDF5 file.
  subroutine write_real_sp_rank_3(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(3)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(3)

    integer, parameter :: dataset_rank = 3

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_sp_rank_3


  !> Write a single precision four rank array to an HDF5 file.
  subroutine write_real_sp_rank_4(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:, :, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(4)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(4)

    integer, parameter :: dataset_rank = 4

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_sp_rank_4


! REAL(DP)


  !> Write a double precision scalar to an HDF5 file.
  subroutine write_real_dp_rank_0(this, h5path, dataset, scalar)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Real(dp) to write.
    real(dp), intent(in), target :: scalar

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(scalar)
    dataset_chunk_shape = [1]
    dataset_shape_local = [1]
    offset_local = [0] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_dp_rank_0


  !> Write a double precision vector to an HDF5 file.
  subroutine write_real_dp_rank_1(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(1)

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_dp_rank_1


  !> Write a double precision two rank array to an HDF5 file.
  subroutine write_real_dp_rank_2(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(2)

    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_dp_rank_2


  !> Write a double precision three rank array to an HDF5 file.
  subroutine write_real_dp_rank_3(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(3)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(3)

    integer, parameter :: dataset_rank = 3

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_dp_rank_3


  !> Write a double precision four rank array to an HDF5 file.
  subroutine write_real_dp_rank_4(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:, :, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(4)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(4)

    integer, parameter :: dataset_rank = 4

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    dataset_shape_local = dataset_shape
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_real_dp_rank_4


! COMPLEX(DP)


  !> Write a double vector array to an HDF5 file.
  subroutine write_complex_dp_rank_1(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(1)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    dataset_shape_local = [[2], dataset_shape]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_complex_dp_rank_1


  !> Write a double complex matrix to an HDF5 file.
  subroutine write_complex_dp_rank_2(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(2)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 3

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    dataset_shape_local = [[2], dataset_shape]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_complex_dp_rank_2


  !> Write a double complex three rank array to an HDF5 file.
  subroutine write_complex_dp_rank_3(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(3)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(3)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 4

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank) 
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double() 
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    dataset_shape_local = [[2], dataset_shape]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)
  end subroutine write_complex_dp_rank_3


  !> Write a double complex four rank array to an HDF5 file.
  subroutine write_complex_dp_rank_4(this, h5path, dataset, dataset_chunk, offset, dataset_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to write the dataset to.
    character(*), intent(in) :: h5path
    !> Name of the dataset to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:, :, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(4)
    !> Shape of the whole array to be written.
    integer, intent(in) :: dataset_shape(4)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 5

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank), dataset_shape_local(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank) 
    type(c_ptr) :: dataset_chunk_ptr

    if (this%exists(join_paths(h5path, dataset))) call this%delete(join_paths(h5path, dataset))

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    dataset_shape_local = [[2], dataset_shape]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), &
                'Some elements of offset < 0.')
    call assert(all(offset_local + dataset_chunk_shape <= dataset_shape_local), &
                'Some elements of offset_local + dataset_chunk_shape > dataset_shape_local.')

    call hdf5_write_dataset(this%h5id,           &
                            trim(h5path),        &
                            trim(dataset),       &
                            data_type_id,        &
                            dataset_chunk_ptr,   &
                            dataset_rank,        &
                            dataset_chunk_shape, &
                            dataset_shape_local, &
                            offset_local, &
                            this%serial_access)                        
  end subroutine write_complex_dp_rank_4
  

!-----------------------------------------------------------------------------------------------------------------------
! READ DATASET


! CHARACTER  


  !> Read an character scalar from an HDF5 file.
  subroutine read_string(this, h5path, dataset, string)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> String to read.
    character(:), allocatable, intent(out) :: string

    integer, parameter :: dataset_rank = 1

    character(hdf5_max_string_len), target :: string_local
    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr
    
    data_type_id = hdf5_character()
    dataset_chunk_ptr = c_loc(string_local)
    dataset_chunk_shape = [hdf5_max_string_len]
    offset_local = [0]

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)

    string = trim(adjustl(string_local))
  end subroutine read_string


! INTEGER


  !> Read an integer scalar from an HDF5 file.
  subroutine read_integer_rank_0(this, h5path, dataset, scalar)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> Integer to read
    integer, intent(in), target :: scalar

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_integer()
    dataset_chunk_ptr = c_loc(scalar)
    dataset_chunk_shape = [1]
    offset_local = [0] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_integer_rank_0


  !> Read a double precision two rank array from an HDF5 file.
  subroutine read_integer_rank_1(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    integer, intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_integer()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_integer_rank_1


  !> Read a double precision two rank array from an HDF5 file.
  subroutine read_integer_rank_2(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    integer, intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)

    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_integer()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_integer_rank_2


! REAL(SP)

  
  !> Read a single precision vector from an HDF5 file.
  subroutine read_real_sp_rank_0(this, h5path, dataset, scalar)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> Real(sp) to read.
    real(sp), intent(in), target :: scalar

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(scalar)
    dataset_chunk_shape = [1]
    offset_local = [0] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_sp_rank_0


  !> Read a single precision vector from an HDF5 file.
  subroutine read_real_sp_rank_1(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_sp_rank_1


  !> Read a single precision two rank array from an HDF5 file.
  subroutine read_real_sp_rank_2(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)

    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_sp_rank_2


  !> Read a single precision three rank array from an HDF5 file.
  subroutine read_real_sp_rank_3(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(3)

    integer, parameter :: dataset_rank = 3

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_sp_rank_3


  !> Read a single precision four rank array from an HDF5 file.
  subroutine read_real_sp_rank_4(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(sp), intent(in), target :: dataset_chunk(:, :, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(4)

    integer, parameter :: dataset_rank = 4

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_float()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_sp_rank_4


! REAL(DP)


  !> Read a double precision vector from an HDF5 file.
  subroutine read_real_dp_rank_0(this, h5path, dataset, scalar)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> Real(dp) to read.
    real(dp), intent(in), target :: scalar

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(scalar)
    dataset_chunk_shape = [1]
    offset_local = [0] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_dp_rank_0


  !> Read a double precision vector from an HDF5 file.
  subroutine read_real_dp_rank_1(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)

    integer, parameter :: dataset_rank = 1

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_dp_rank_1


  !> Read a double precision two rank array from an HDF5 file.
  subroutine read_real_dp_rank_2(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)

    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_dp_rank_2


  !> Read a double precision three rank array from an HDF5 file.
  subroutine read_real_dp_rank_3(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(3)

    integer, parameter :: dataset_rank = 3

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_dp_rank_3


  !> Read a double precision four rank array from an HDF5 file.
  subroutine read_real_dp_rank_4(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(dp), intent(in), target :: dataset_chunk(:, :, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(4)

    integer, parameter :: dataset_rank = 4

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = shape(dataset_chunk)
    offset_local = offset - 1 ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_real_dp_rank_4


! COMPLEX(DP)


  !> Read a double complex vector from an HDF5 file.
  subroutine read_complex_dp_rank_1(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(1)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 2

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_complex_dp_rank_1


  !> Read a double complex matrix from an HDF5 file.
  subroutine read_complex_dp_rank_2(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(2)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 3

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank)
    type(c_ptr) :: dataset_chunk_ptr

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_complex_dp_rank_2


  !> Read a double complex three rank array from an HDF5 file.
  subroutine read_complex_dp_rank_3(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(3)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 4

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank) 
    type(c_ptr) :: dataset_chunk_ptr 

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    offset_local = [[0], offset - 1] ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_complex_dp_rank_3


  !> Read a double complex four rank array from an HDF5 file.
  subroutine read_complex_dp_rank_4(this, h5path, dataset, dataset_chunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the HDF5 file to the group to read the dataset from.
    character(*), intent(in) :: h5path
    !> Name of the dataset to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    complex(dp), intent(in), target :: dataset_chunk(:, :, :, :)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer, intent(in) :: offset(4)

    ! HDF5 has no complex type. Real and imaginary part are written to seperate dimensions.
    ! => rank(dataset) = rank(array) + 1
    integer, parameter :: dataset_rank = 5

    integer(hdf5_id) :: data_type_id
    integer(hdf5_size) :: dataset_chunk_shape(dataset_rank)
    integer(hdf5_ssize) :: offset_local(dataset_rank) 
    type(c_ptr) :: dataset_chunk_ptr 

    data_type_id = hdf5_double()
    dataset_chunk_ptr = c_loc(dataset_chunk)
    dataset_chunk_shape = [[2], shape(dataset_chunk)]
    offset_local = [[0], offset - 1]  ! Convert Fortran indexing to C indexing

    call assert(all(offset_local >= 0), 'Some elements of offset < 0.')

    call hdf5_read_dataset(this%h5id,           &
                           trim(h5path),        &
                           trim(dataset),       &
                           data_type_id,        &
                           dataset_chunk_ptr,   &
                           dataset_rank,        &
                           dataset_chunk_shape, &
                           offset_local, &
                            this%serial_access)
  end subroutine read_complex_dp_rank_4

end module xhdf5  
