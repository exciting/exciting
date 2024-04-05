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
!> Be aware of that the writing routines will overwrite an existing data set with the same name.
module xhdf5
  use iso_c_binding, only: c_ptr, c_loc

  use os_utils, only: path_exists, join_paths
  use precision, only: sp, dp, i32
  use modmpi, only: mpiinfo

  use xhdf5_globals
  use xhdf5_error_handling
  use mpi_utils
  use hdf5_file
  use hdf5_group
  use hdf5_write
  use hdf5_read
  use hdf5_dataset_utils
  use hdf5_datatype
  use datainfo_utils


  implicit none


  private
  public :: xhdf5_type, abort_if_not_hdf5, h5group_root


  !> Type for handeling an HDF5 file.
  type xhdf5_type
    !> OS file path.
    character(:), allocatable :: path
    !> HDF5 file id.
    integer(hdf5_id) :: h5id = h5id_undefined
    !> MPI communicator handle
    type(mpi_comm_type) :: mpi_comm
    !> Flag to initialize HDF5 file in serial mode, even if compiled with MPI.
    !> By default this is set to false.
    logical :: serial_access = .false.
    !> Number of overwritten datasets
    integer :: n_overwritten_datasets = 0

    contains

    generic :: initialize => initialize_mpi_comm_type, initialize_mpiinfo
    procedure :: initialize_mpi_comm_type, initialize_mpiinfo

    procedure :: finalize
    procedure :: initialize_group
    procedure :: dataset_shape

    procedure :: exists => link_exists

    procedure :: delete => delete_link

    procedure :: handle_if_dataset_exists
    procedure :: evaluate_overwritten_datasets 

    generic :: write => write_string, write_bool, write_integer, &
                        write_float, write_double, &
                        write_float_complex, write_double_complex

    procedure :: write_string, write_bool, write_integer, &
                 write_float, write_double, &
                 write_float_complex, write_double_complex

    generic :: read => read_string, read_bool, read_integer, &
                       read_float, read_double, &
                       read_float_complex, read_double_complex

    procedure :: read_string, read_bool, read_integer, &
                 read_float, read_double, &
                 read_float_complex, read_double_complex

  end type xhdf5_type


  contains


  !> Initialize HDF5 library and, if `path` exists, open the file, else create a file at `path`.
  subroutine initialize_mpi_comm_type(this, path, mpi_comm, serial_access)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Relative path to HDF5 file.
    character(*), intent(in) :: path
    !> MPI communicator.
    type(mpi_comm_type), intent(in) :: mpi_comm
    !> Set to `.true.` for serial access in an MPI environment. 
    !> This is only allowed if the root process is the caller.
    logical, intent(in), optional :: serial_access

    this%path = trim(path)
    this%mpi_comm = mpi_comm

    if(present(serial_access)) then
      this%serial_access = serial_access
    end if

    if(this%serial_access) then
      call xhdf5_assert(this%mpi_comm, comm_to_rank(this%mpi_comm) == root_rank, &
              'Error(xhdf5%initialize_mpi_comm_type): serial_access is set to .true., only the root process &
              is allowed to call xhdf5 routines.')
    end if

    call hdf5_initialize(this%mpi_comm)

    if(path_exists(this%path)) then
      call hdf5_open_file(this%path, this%mpi_comm, this%h5id, this%serial_access)
    else 
      call hdf5_create_file(this%path, this%mpi_comm, this%h5id, this%serial_access)
    end if
  end subroutine


  !> Initialize HDF5 library and, if `path` exists, open the file, else create a file at `path`.
  subroutine initialize_mpiinfo(this, path, mpi_comm, serial_access)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Relative path to HDF5 file.
    character(*), intent(in) :: path
    !> MPI communicator.
    type(mpiinfo), intent(in) :: mpi_comm
    !> Set to `.true.` for serial access in an MPI environment. 
    !> This is only allowed if the root process is the caller.
    logical, intent(in), optional :: serial_access

    logical :: serial_access_local 

    serial_access_local = .false.
    if(present(serial_access)) serial_access_local = serial_access

    call this%initialize(path, mpi_comm_type(mpi_comm%comm), serial_access=serial_access_local)
  end subroutine


  !> Close file and finalize HDF5 library.
  subroutine finalize(this)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this

    if (this%h5id == h5id_undefined) return
    call hdf5_close_file(this%mpi_comm, this%h5id)
    call hdf5_finalize(this%mpi_comm)
    call this%evaluate_overwritten_datasets()
    this%h5id = h5id_undefined
  end subroutine finalize 
  

  !> Create a new group with name `group` at `h5path`. If the group already exists, the routine does nothing.
  subroutine initialize_group(this, h5path, groupname)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path
    !> Group name.
    character(*), intent(in) :: groupname

    if (this%exists(join_paths(h5path, groupname))) return
    call hdf5_create_group(this%mpi_comm, this%h5id, trim(h5path), groupname)
  end subroutine initialize_group
  

  !> Return true or false, whether the link to `h5path` exists or not.
  logical function link_exists(this, h5path)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path

    link_exists = hdf5_link_exists(this%mpi_comm, this%h5id, trim(h5path))
  end function link_exists

  !> Delete link to `h5path`. This action frees the space on the drive allocated by the object `h5path` points to.
  subroutine delete_link(this, h5path)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path

    call hdf5_delete_link(this%mpi_comm, this%h5id, trim(h5path))
  end subroutine delete_link


  !> Get the shape of a data set.
  subroutine dataset_shape(this, h5path, datasetname, shape, complex_dataset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path
    !> Dataset name.
    character(*), intent(in) :: datasetname
    !> Dataset shape
    integer, allocatable :: shape(:)
    !> Set to `.true.` if the data set is complex. If not, the first dimension will be 
    !> always 2 for complex datasets.
    logical, optional :: complex_dataset

    logical :: complex_dataset_local
    integer(hdf5_size), allocatable :: shape_local(:)

    complex_dataset_local = .false.
    if(present(complex_dataset)) complex_dataset_local = complex_dataset

    call hdf5_get_dataset_shape(this%mpi_comm, this%h5id, h5path, datasetname, shape_local)

    if(complex_dataset_local) then 
      shape = shape_local(2:)
    else 
      shape = shape_local
    end if
  end subroutine dataset_shape


  !> Warn about overwritten datasets.
  subroutine evaluate_overwritten_datasets(this)
    !> HDF5 file handler.
    class(xhdf5_type), intent(in) :: this

#ifdef USE_ASSERT
    if(this%n_overwritten_datasets > 0) then
      write(error_unit, '(1x, A, I3, A, A, A)') &
        "Warning(xhdf5): ", this%n_overwritten_datasets, " dataset(s) were overwritten in file ", this%path, "."
    end if 
#endif 
  end subroutine evaluate_overwritten_datasets


  !> Handles what to do if a data set already exists.
  !> If this is the case, the data set will be overwritten with new data.
  !> In debug mode, a warning is printed to the terminal.
  subroutine handle_if_dataset_exists(this, h5path, datasetname)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the hdf5 file.
    character(*), intent(in) :: h5path
    !> Dataset name.
    character(*), intent(in) :: datasetname

    if (this%exists(join_paths(h5path, datasetname))) then
      this%n_overwritten_datasets = this%n_overwritten_datasets + 1
      call this%delete(join_paths(h5path, datasetname))

#ifdef USE_ASSERT
      write(error_unit, '(1x, A, A, A)') &
          "Warning(xhdf5): Dataset ", datasetname, " already existed and was deleted. &
          It will be overwritten with new data."
#endif
    end if

  end subroutine

!-----------------------------------------------------------------------------------------------------------------------
! WRITE DATASET

  !> Write a string to an HDF5 file.
  subroutine write_string(this, h5path, dataset, string)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> String to write.
    character(*), intent(in), target :: string

    integer(hdf5_id) :: string_type
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo

    call datainfo%initialize([1], [-1], [-1], .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'write_integer')

    call this%handle_if_dataset_exists(h5path, dataset)

    string_type = hdf5_string(this%mpi_comm, len(string, kind=hdf5_size))
    datachunk_ptr = c_loc(string)
    call hdf5_write_dataset(this%mpi_comm, this%h5id, trim(h5path), dataset, string_type, datachunk_ptr, datainfo, this%serial_access)
  end subroutine write_string


  !> Write a boolean to an HDF5 file.
  !> The routine does not really writes a boolean but a string representing the boolean. 
  !> For '.true.' it writes `'true'` and for `.false.` `'false'`.
  subroutine write_bool(this, h5path, dataset, bool)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Boolean to write.
    logical, intent(in), target :: bool

    character(:), allocatable :: bool_string 

    if(bool) then
      bool_string = true_string
    else 
      bool_string = false_string
    end if

    call write_string(this, h5path, dataset, bool_string)
  end subroutine write_bool


  !> Write an integer array to an HDF5 file.
  subroutine write_integer(this, h5path, dataset, datachunk, offset, global_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    integer(i32), intent(in), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)
    !> Shape of the whole array to be written.
    integer(i32), intent(in), optional :: global_shape(:)

    integer(i32), allocatable :: global_shape_local(:), offset_local(:)
    type(datainfo_type) :: datainfo
    type(c_ptr) :: datachunk_ptr
    
    offset_local = [-1]
    if(present(offset)) offset_local = offset

    global_shape_local = [-1]
    if(present(global_shape)) global_shape_local = global_shape

    call datainfo%initialize(shape(datachunk), global_shape_local, offset_local, .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'write_integer')

    call this%handle_if_dataset_exists(h5path, dataset)

    datachunk_ptr = c_loc(datachunk)
    call hdf5_write_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_integer(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine write_integer


  !> Write an float array to an HDF5 file.
  subroutine write_float(this, h5path, dataset, datachunk, offset, global_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(sp), intent(in), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)
    !> Shape of the whole array to be written.
    integer(i32), intent(in), optional :: global_shape(:)

    integer(i32), allocatable :: global_shape_local(:), offset_local(:)
    type(datainfo_type) :: datainfo
    type(c_ptr) :: datachunk_ptr
    
    offset_local = [-1]
    if(present(offset)) offset_local = offset

    global_shape_local = [-1]
    if(present(global_shape)) global_shape_local = global_shape

    call datainfo%initialize(shape(datachunk), global_shape_local, offset_local, .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'write_float')

    call this%handle_if_dataset_exists(h5path, dataset)

    datachunk_ptr = c_loc(datachunk)
    call hdf5_write_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_float(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine write_float


  !> Write a double array to an HDF5 file.
  subroutine write_double(this, h5path, dataset, datachunk, offset, global_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    real(dp), intent(in), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)
    !> Shape of the whole array to be written.
    integer(i32), intent(in), optional :: global_shape(:)

    integer(i32), allocatable :: global_shape_local(:), offset_local(:)
    type(datainfo_type) :: datainfo
    type(c_ptr) :: datachunk_ptr
    
    offset_local = [-1]
    if(present(offset)) offset_local = offset

    global_shape_local = [-1]
    if(present(global_shape)) global_shape_local = global_shape

    call datainfo%initialize(shape(datachunk), global_shape_local, offset_local, .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'write_double')

    call this%handle_if_dataset_exists(h5path, dataset)

    datachunk_ptr = c_loc(datachunk)
    call hdf5_write_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_double(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine write_double


  !> Write a double complex array to an HDF5 file.
  subroutine write_float_complex(this, h5path, dataset, datachunk, offset, global_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    complex(sp), intent(in), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)
    !> Shape of the whole array to be written.
    integer(i32), intent(in), optional :: global_shape(:)

    integer(i32), allocatable :: global_shape_local(:), offset_local(:)
    type(datainfo_type) :: datainfo
    type(c_ptr) :: datachunk_ptr
    
    offset_local = [-1]
    if(present(offset)) offset_local = offset

    global_shape_local = [-1]
    if(present(global_shape)) global_shape_local = global_shape

    call datainfo%initialize(shape(datachunk), global_shape_local, offset_local, .true.)
    call datainfo%sanity_checks(this%mpi_comm, 'write_float_complex')

    call this%handle_if_dataset_exists(h5path, dataset)

    datachunk_ptr = c_loc(datachunk)
    call hdf5_write_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_float(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine write_float_complex  



  !> Write a double complex array to an HDF5 file.
  subroutine write_double_complex(this, h5path, dataset, datachunk, offset, global_shape)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Array chunk handled by the current MPI rank.
    complex(dp), intent(in), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)
    !> Shape of the whole array to be written.
    integer(i32), intent(in), optional :: global_shape(:)

    integer(i32), allocatable :: global_shape_local(:), offset_local(:)
    type(datainfo_type) :: datainfo
    type(c_ptr) :: datachunk_ptr
    
    offset_local = [-1]
    if(present(offset)) offset_local = offset

    global_shape_local = [-1]
    if(present(global_shape)) global_shape_local = global_shape

    call datainfo%initialize(shape(datachunk), global_shape_local, offset_local, .true.)
    call datainfo%sanity_checks(this%mpi_comm, 'write_double_complex')

    call this%handle_if_dataset_exists(h5path, dataset)

    datachunk_ptr = c_loc(datachunk)
    call hdf5_write_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_double(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine write_double_complex
  

!-----------------------------------------------------------------------------------------------------------------------
! READ DATASET


  !> Read a character scalar from an HDF5 file.
  subroutine read_string(this, h5path, dataset, string)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to read the data set from.
    character(*), intent(in) :: h5path
    !> Name of the data set to read from.
    character(*), intent(in) :: dataset
    !> String to read.
    character(:), allocatable, intent(out), target :: string

    integer, parameter :: dataset_rank = 1

    integer(hdf5_size) :: string_len
    integer(hdf5_id) :: dataset_type
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo
    
    dataset_type = hdf5_get_type(this%mpi_comm, this%h5id, h5path, dataset)
    allocate(character(len=hdf5_get_type_size(this%mpi_comm, dataset_type)) :: string)
    call datainfo%initialize([1], [-1], [-1], .false.)
    
    datachunk_ptr = c_loc(string)
    call hdf5_read_dataset(this%mpi_comm, this%h5id, h5path, dataset, dataset_type, datachunk_ptr, datainfo, this%serial_access)
  end subroutine read_string


  !> Read a boolean to an HDF5 file.
  !> The boolean is not saved as a boolean in the file but as a string representing the state
  !> of the boolean. See [[read_bool]].
  subroutine read_bool(this, h5path, dataset, bool)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to write the data set to.
    character(*), intent(in) :: h5path
    !> Name of the data set to write.
    character(*), intent(in) :: dataset
    !> Boolean to read.
    logical, intent(out), target :: bool

    character(:), allocatable :: bool_string

    call read_string(this, h5path, dataset, bool_string)

    if (bool_string == true_string) then
      bool = .true.
    else if(bool_string == false_string) then 
      bool = .false.
    else 
      call xhdf5_assert(this%mpi_comm, .false., 'dataset does not contain a string representing the state of a bool.')
    end if 
  end subroutine read_bool


  !> Read an integer array from an HDF5 file.
  subroutine read_integer(this, h5path, dataset, datachunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to read the data set from.
    character(*), intent(in) :: h5path
    !> Name of the data set to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    integer(i32), intent(in), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)

    integer(i32), allocatable :: offset_local(:), global_shape(:)
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo

    offset_local = [-1]
    if(present(offset)) offset_local = offset

    call this%dataset_shape(h5path, dataset, global_shape)
    call datainfo%initialize(shape(datachunk), global_shape, offset_local, .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'read_integer')

    datachunk_ptr = c_loc(datachunk)
    call hdf5_read_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_integer(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine read_integer


  !> Read a float array from an HDF5 file.
  subroutine read_float(this, h5path, dataset, datachunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to read the data set from.
    character(*), intent(in) :: h5path
    !> Name of the data set to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(sp), intent(out), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)

    integer(i32), allocatable :: offset_local(:), global_shape(:)
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo

    offset_local = [-1]
    if(present(offset)) offset_local = offset

    call this%dataset_shape(h5path, dataset, global_shape)
    call datainfo%initialize(shape(datachunk), global_shape, offset_local, .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'read_float')

    datachunk_ptr = c_loc(datachunk)
    call hdf5_read_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_float(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine read_float


  !> Read a double array from an HDF5 file.
  subroutine read_double(this, h5path, dataset, datachunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to read the data set from.
    character(*), intent(in) :: h5path
    !> Name of the data set to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    real(dp), intent(out), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)

    integer(i32), allocatable :: offset_local(:), global_shape(:)
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo

    offset_local = [-1]
    if(present(offset)) offset_local = offset

    call this%dataset_shape(h5path, dataset, global_shape)
    call datainfo%initialize(shape(datachunk), global_shape, offset_local, .false.)
    call datainfo%sanity_checks(this%mpi_comm, 'read_double')

    datachunk_ptr = c_loc(datachunk)
    call hdf5_read_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_double(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine read_double


  !> Read a double complex array from an HDF5 file.
  subroutine read_float_complex(this, h5path, dataset, datachunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to read the data set from.
    character(*), intent(in) :: h5path
    !> Name of the data set to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    complex(sp), intent(out), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)

    integer(i32), allocatable :: offset_local(:), global_shape(:)
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo

    offset_local = [-1]
    if(present(offset)) offset_local = offset

    call this%dataset_shape(h5path, dataset, global_shape, .true.)
    call datainfo%initialize(shape(datachunk), global_shape, offset_local, .true.)
    call datainfo%sanity_checks(this%mpi_comm, 'read_float')

    datachunk_ptr = c_loc(datachunk)
    call hdf5_read_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_float(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine read_float_complex


  !> Read a double complex array from an HDF5 file.
  subroutine read_double_complex(this, h5path, dataset, datachunk, offset)
    !> HDF5 file handler.
    class(xhdf5_type), intent(inout) :: this
    !> Absolute path in the HDF5 file to the group to read the data set from.
    character(*), intent(in) :: h5path
    !> Name of the data set to read from.
    character(*), intent(in) :: dataset
    !> On output, contains the array chunk, handled by the current MPI rank.
    complex(dp), intent(out), contiguous, target :: datachunk(..)
    !> Offset of array chunk handled by the current MPI rank in the whole array. 
    integer(i32), intent(in), optional :: offset(:)

    integer(i32), allocatable :: offset_local(:), global_shape(:)
    type(c_ptr) :: datachunk_ptr
    type(datainfo_type) :: datainfo

    offset_local = [-1]
    if(present(offset)) offset_local = offset

    call this%dataset_shape(h5path, dataset, global_shape, .true.)
    call datainfo%initialize(shape(datachunk), global_shape, offset_local, .true.)
    call datainfo%sanity_checks(this%mpi_comm, 'read_double_complex')

    datachunk_ptr = c_loc(datachunk)
    call hdf5_read_dataset(this%mpi_comm, this%h5id, h5path, dataset, hdf5_double(), datachunk_ptr, datainfo, this%serial_access)
  end subroutine read_double_complex

end module xhdf5  
