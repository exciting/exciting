module block_data_file
  use iso_fortran_env, only: int64
  use precision, only: dp
  use asserts, only: assert
  use modmpi

  implicit none
  private

  !> Class for managing serial and parallel binary data i/o.
  !> Write and read block data, saved as arrays, of any type and shape to binary files.
  type :: block_data_file_type
    private
    !> MPI file handler or i/o unit
    integer :: fid = 0
    !> path to the files location
    character(:), allocatable :: path
    !> `.true.` if file is currently open
    logical :: file_is_open = .false.
    !> rank of a data block
    integer :: block_rank = 0
    !> shaoe of a data block
    integer, allocatable :: block_shape(:)
    !> record length (in system dependent units)
    !> (a single record consists of the variables `block_rank`, `block_shape` (for consistency checks)
    !> and `product(block_shape)` data entries of the specified data type)
    integer(kind=int64) :: record_length = 0
    !> dummy element for data type determination
    class(*), allocatable :: type_dummy
      
    contains
      private

      procedure, public :: exists      => block_data_file_exists
      procedure, public :: is_open     => block_data_file_is_open
      procedure, public :: delete      => block_data_file_delete
      procedure, public :: open        => block_data_file_open
      procedure, public :: close       => block_data_file_close

      procedure, public :: get_path, get_file_id
      procedure, public :: get_block_rank, get_block_shape

      procedure :: read_integer, read_double, read_double_complex
      generic, public :: read => read_integer, read_double, read_double_complex

      procedure :: write_integer, write_double, write_double_complex
      generic, public :: write => write_integer, write_double, write_double_complex
  end type block_data_file_type
  !> constructor
  interface block_data_file_type
    procedure :: setup_block_data_file_type
  end interface

  public :: block_data_file_type

  contains
    
    !> constructor for block data file type
    function setup_block_data_file_type(file_path, block_shape, type_dummy) result(this)
      !> path to the files location
      character(*), intent(in) :: file_path
      !> dimensions (number of elements per dimension) of data block
      integer, intent(in) :: block_shape(:)
      !> single element of data to determine size
      class(*), intent(in) :: type_dummy
      type(block_data_file_type) :: this

      integer :: ierr, intsize, record_length
      character(:), allocatable :: errmsg

      errmsg = ''

      ! set file path
      this%path = trim(file_path)
      ! set rank and dimensions of data block
      this%block_rank = size(block_shape)
      this%block_shape = block_shape
      ! set data type dummy element
      this%type_dummy = type_dummy
      ! set record length
#ifdef MPI
      ! for parallel i/o the record length is the block size in bytes
      call MPI_type_size( MPI_INTEGER, intsize, ierr)
      call mpi_error_to_string('Failed to retreive MPI_INTEGER size.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(setup_block_data_file_type) '//errmsg)
      this%record_length = intsize + this%block_rank*intsize ! for rank and dimensions
      this%record_length = this%record_length + product(this%block_shape)*(storage_size(this%type_dummy)/8) ! for data block
#else
      ! for serial i/o the record length is the block size in compiler dependent units
      inquire(iolength=this%record_length) this%block_rank, this%block_shape ! for rank and dimensions
      dummytype: select type( t => this%type_dummy)
        type is(integer)
          inquire( iolength=record_length) t
        type is(real(dp))
          inquire( iolength=record_length) t
        type is(complex(dp))
          inquire( iolength=record_length) t
        class default
          call terminate_if_false( .false., &
                 '(setup_block_data_file_type) Unsupported data dype.')
      end select dummytype
      this%record_length = this%record_length + product(this%block_shape)*record_length ! for data block
#endif
    end function setup_block_data_file_type

    !> open block data file for serial or parallel i/o; 
    !> if file doesn't exist, it will be created
    subroutine block_data_file_open(this, mpi_comm, delete_existing, mpi_access_mode)
      use file_utils, only: file_is_open
      use m_getunit
      class(block_data_file_type), intent(inout) :: this
      !> MPI communicator
      type(mpiinfo), intent(inout) :: mpi_comm
      !> if file already exists then delete it and create new one (default: `.false.`)
      logical, optional, intent(in) :: delete_existing
      !> MPI access mode (default: `MPI_MODE_CREATE+MPI_MODE_RDWR`)
      integer, optional, intent(in) :: mpi_access_mode

      logical :: delete
      integer :: ierr, mode
      character(:), allocatable :: errmsg

      errmsg = ''

      delete = .false.
      if (present(delete_existing)) delete = delete_existing

      ! return if file was already opened
      call terminate_if_false( .not. this%is_open(), '(block_data_file_open) &
        The file you try to has already been opened.' )

      ! delete existing file if necessary
      if (delete .and. this%exists()) call this%delete(mpi_comm)

#ifdef MPI
      ! open file for parallel i/o
      mode = MPI_MODE_CREATE+MPI_MODE_RDWR
      if (present(mpi_access_mode)) mode = mpi_access_mode
      call MPI_file_open( mpi_comm%comm, trim(this%path), mode, MPI_INFO_NULL, this%fid, ierr)
#else
      ! open file for serial i/o
      call getunit(this%fid)
      open(this%fid, file=trim(this%path), action='readwrite', form='unformatted', &
           access='direct', recl=this%record_length, iostat=ierr)
#endif

      call mpi_error_to_string('Failed to open file.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(block_data_file_open) '//errmsg)

#ifdef MPI
      this%file_is_open = (this%file_is_open .or. (ierr == 0))
#else
      this%file_is_open = (this%file_is_open .or. file_is_open(this%fid, ierr))
#endif
    end subroutine block_data_file_open

    !> close block data file
    subroutine block_data_file_close(this, mpi_comm)
      use file_utils, only: close_file
      class(block_data_file_type), intent(inout) :: this
      !> MPI communicator
      type(mpiinfo), intent(inout) :: mpi_comm

      integer :: ierr
      character(:), allocatable :: errmsg

      errmsg = ''

      call close_file(this%fid, mpi_comm, ierr)
      call mpi_error_to_string('Failed to close file.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(block_data_file_close) '//errmsg)
      this%file_is_open = (this%file_is_open .and. ierr /= 0)
      this%fid = 0
    end subroutine block_data_file_close

    !> delete block data file
    subroutine block_data_file_delete(this, mpi_comm)
      use file_utils, only: delete_file
      class(block_data_file_type), intent(inout) :: this
      !> MPI communicator
      type(mpiinfo), intent(inout) :: mpi_comm

      integer :: ierr
      character(:), allocatable :: errmsg

      errmsg = ''

      if( this%file_is_open ) call this%close(mpi_comm)
      call delete_file(trim(this%path), mpi_comm, ierr)
      call mpi_error_to_string('Failed to delete file.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(block_data_file_delete) '//errmsg)
    end subroutine block_data_file_delete

    !> check if block data file exists
    function block_data_file_exists(this) result(l)
      use os_utils, only: path_exists
      class(block_data_file_type), intent(inout) :: this

      integer :: ierr
      logical :: l

      l = path_exists(trim(this%path), ierr)
    end function block_data_file_exists

    !> check if block data file is open
    function block_data_file_is_open(this) result(l)
      class(block_data_file_type), intent(in) :: this
      logical :: l
      l = this%file_is_open
    end function block_data_file_is_open

    !> get block data file path
    function get_path(this) result(path)
      class(block_data_file_type), intent(in) :: this
      character(:), allocatable :: path
      path = this%path
    end function get_path

    !> get block data file file id
    function get_file_id(this) result(fid)
      class(block_data_file_type), intent(in) :: this
      integer :: fid
      fid = this%fid
    end function get_file_id

    !> get block data file block rank
    function get_block_rank(this) result(block_rank)
      class(block_data_file_type), intent(in) :: this
      integer :: block_rank
      block_rank = this%block_rank
    end function get_block_rank

    !> get block data file block shape
    function get_block_shape(this) result(block_shape)
      class(block_data_file_type), intent(in) :: this
      integer, allocatable :: block_shape(:)
      block_shape = this%block_shape
    end function get_block_shape

    !> write integer data block to file
    subroutine write_integer(this, record, data_block)
      use iso_c_binding
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> data block to write
      integer, target, intent(in) :: data_block(..)

      integer :: ierr, n
      integer, allocatable :: shp(:)
      integer, pointer :: buffer(:)
      character(:), allocatable :: errmsg

      ! check if file object hosts integers
      dummytype: select type( t => this%type_dummy)
        type is(integer)
          exit dummytype
        class default
          call assert( .false., 'block_data_file_type object not set up to handle data of type integer.')
      end select dummytype

      n = product(this%block_shape)

      ! check block size
      allocate( shp, source=shape( data_block ) )
      call assert( size( shp ) == this%block_rank, &
        'Rank of `data_block` does not match rank defined for `block_data_file_type` object.' )
      call assert( all( shp == this%block_shape ), &
        'Shape of `data_block` does not match shape defined for `block_data_file_type` object.' )

      ! write rank and shape
      call write_rank_and_shape(this, record)

      ! write data
      errmsg = 'Error in writing integer data block to file'
      call c_f_pointer( c_loc( data_block ), buffer, [n] )
#ifdef MPI
      call MPI_file_write_at(this%fid, get_MPI_data_offset(this, record), buffer, n, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#else
      write(this%fid, rec=record, iostat=ierr) this%block_rank, this%block_shape, buffer
#endif
      call mpi_error_to_string('', ierr, errmsg)

      call terminate_if_false( ierr == 0, '(write_integer) '//errmsg )
    end subroutine write_integer

    !> write double real data block to file
    subroutine write_double(this, record, data_block)
      use iso_c_binding
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> data block to read
      real(dp), target, intent(in) :: data_block(..)

      integer :: ierr, n
      integer, allocatable :: shp(:)
      real(dp), pointer :: buffer(:)
      character(:), allocatable :: errmsg

      ! check if file object hosts integers
      dummytype: select type( t => this%type_dummy)
        type is(real(dp))
          exit dummytype
        class default
          call assert( .false., 'block_data_file_type object not set up to handle data of type double.')
      end select dummytype

      n = product(this%block_shape)

      ! check block size
      allocate( shp, source=shape( data_block ) )
      call assert( size( shp ) == this%block_rank, &
        'Rank of `data_block` does not match rank defined for `block_data_file_type` object.' )
      call assert( all( shp == this%block_shape ), &
        'Shape of `data_block` does not match shape defined for `block_data_file_type` object.' )

      ! write rank and shape
      call write_rank_and_shape(this, record)

      ! write data
      errmsg = 'Error in writing real data block to file'
      call c_f_pointer( c_loc( data_block ), buffer, [n] )
#ifdef MPI
      call MPI_file_write_at(this%fid, get_MPI_data_offset(this, record), buffer, n, MPI_DOUBLE, MPI_STATUS_IGNORE, ierr)
#else
      write(this%fid, rec=record, iostat=ierr) this%block_rank, this%block_shape, buffer
#endif
      call mpi_error_to_string('', ierr, errmsg)

      call terminate_if_false( ierr == 0, '(write_double) '//errmsg )
    end subroutine write_double

    !> write double complex data block to file
    subroutine write_double_complex(this, record, data_block)
      use iso_c_binding
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> data block to read
      complex(dp), target, intent(in) :: data_block(..)

      integer :: ierr, n
      integer, allocatable :: shp(:)
      complex(dp), pointer :: buffer(:)
      character(:), allocatable :: errmsg

      ! check if file object hosts integers
      dummytype: select type( t => this%type_dummy)
        type is(complex(dp))
          exit dummytype
        class default
          call assert( .false., 'block_data_file_type object not set up to handle data of type double complex.')
      end select dummytype

      ! check block size
      allocate( shp, source=shape( data_block ) )
      call assert( size( shp ) == this%block_rank, &
        'Rank of `data_block` does not match rank defined for `block_data_file_type` object.' )
      call assert( all( shp == this%block_shape ), &
        'Shape of `data_block` does not match shape defined for `block_data_file_type` object.' )

      n = product(this%block_shape)

      ! write rank and shape
      call write_rank_and_shape(this, record)

      ! write data
      errmsg = 'Error in writing complex data block to file'
      call c_f_pointer( c_loc( data_block ), buffer, [n] )
#ifdef MPI
      call MPI_file_write_at(this%fid, get_MPI_data_offset(this, record), buffer, n, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
#else
      write(this%fid, rec=record, iostat=ierr) this%block_rank, this%block_shape, buffer
#endif

      call mpi_error_to_string('', ierr, errmsg)

      call terminate_if_false( ierr == 0, '(write_double_complex) '//errmsg )
    end subroutine write_double_complex

    !> read integer data block from file
    subroutine read_integer(this, record, data_block)
      use iso_c_binding
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> data block to read
      integer, target, intent(out) :: data_block(..)

      integer :: ierr, n, block_rank, block_shape(this%block_rank)
      integer, allocatable :: shp(:)
      integer, pointer :: buffer(:)
      character(:), allocatable :: errmsg

      ! check if file object hosts integers
      dummytype: select type( t => this%type_dummy)
        type is(integer)
          exit dummytype
        class default
          call terminate_if_false( .false., '(read_integer) Block data file not initialized for integers.')
      end select dummytype

      n = product(this%block_shape)

      ! check block size
      allocate( shp, source=shape( data_block ) )
      call assert( size( shp ) == this%block_rank, &
        'Rank of `data_block` does not match rank defined for `block_data_file_type` object.' )
      call assert( all( shp == this%block_shape ), &
        'Shape of `data_block` does not match shape defined for `block_data_file_type` object.' )

      ! read block rank and shape
      call read_rank_and_shape(this, record, block_rank, block_shape)

      ! read data
      errmsg = 'Error in reading integer data block from file'
      call c_f_pointer( c_loc( data_block ), buffer, [n] )
#ifdef MPI
      call MPI_file_read_at(this%fid, get_MPI_data_offset(this, record), buffer, n, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
#else
      read(this%fid, rec=record, iostat=ierr) block_rank, block_shape, buffer
#endif

      call mpi_error_to_string('', ierr, errmsg)

      call terminate_if_false( ierr == 0, '(read_integer) '//errmsg )
    end subroutine read_integer

    !> read double real data block from file
    subroutine read_double(this, record, data_block)
      use iso_c_binding
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> data block to read
      real(dp), target, intent(out) :: data_block(..)

      integer :: ierr, n, block_rank, block_shape(this%block_rank)
      integer, allocatable :: shp(:)
      real(dp), pointer :: buffer(:)
      character(:), allocatable :: errmsg

      ! check if file object hosts integers
      dummytype: select type( t => this%type_dummy)
        type is(real(dp))
          exit dummytype
        class default
          call terminate_if_false( .false., '(read_double) Block data file not initialized for double reals.')
      end select dummytype

      n = product(this%block_shape)

      ! check block size
      allocate( shp, source=shape( data_block ) )
      call assert( size( shp ) == this%block_rank, &
        'Rank of `data_block` does not match rank defined for `block_data_file_type` object.' )
      call assert( all( shp == this%block_shape ), &
        'Shape of `data_block` does not match shape defined for `block_data_file_type` object.' )

      ! read block rank and shape
      call read_rank_and_shape(this, record, block_rank, block_shape)

      ! read data
      errmsg = 'Error in reading real data block from file'
      call c_f_pointer( c_loc( data_block ), buffer, [n] )
#ifdef MPI
      call MPI_file_read_at(this%fid, get_MPI_data_offset(this, record), buffer, n, MPI_DOUBLE, MPI_STATUS_IGNORE, ierr)
#else
      read(this%fid, rec=record, iostat=ierr) block_rank, block_shape, buffer
#endif

      call mpi_error_to_string('', ierr, errmsg)

      call terminate_if_false( ierr == 0, '(read_duble) '//errmsg )
    end subroutine read_double

    !> read double complex data block from file
    subroutine read_double_complex(this, record, data_block)
      use iso_c_binding
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> data block to read
      complex(dp), target, intent(out) :: data_block(..)

      integer :: ierr, n, block_rank, block_shape(this%block_rank)
      integer, allocatable :: shp(:)
      complex(dp), pointer :: buffer(:)
      character(:), allocatable :: errmsg

      ! check if file object hosts integers
      dummytype: select type( t => this%type_dummy)
        type is(complex(dp))
          exit dummytype
        class default
          call terminate_if_false( .false., '(read_double_complex) Block data file not initialized for double reals.')
      end select dummytype

      n = product(this%block_shape)

      ! check block size
      allocate( shp, source=shape( data_block ) )
      call assert( size( shp ) == this%block_rank, &
        'Rank of `data_block` does not match rank defined for `block_data_file_type` object.' )
      call assert( all( shp == this%block_shape ), &
        'Shape of `data_block` does not match shape defined for `block_data_file_type` object.' )

      ! read block rank and shape
      call read_rank_and_shape(this, record, block_rank, block_shape)

      ! read data
      errmsg = 'Error in reading complex data block from file'
      call c_f_pointer( c_loc( data_block ), buffer, [n] )
#ifdef MPI
      call MPI_file_read_at(this%fid, get_MPI_data_offset(this, record), buffer, n, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
#else
      read(this%fid, rec=record, iostat=ierr) block_rank, block_shape, buffer
#endif

      call mpi_error_to_string('', ierr, errmsg)

      call terminate_if_false( ierr == 0, '(read_double_complex) '//errmsg )
    end subroutine read_double_complex

    !> get file offset for record data for parallel i/o
    function get_MPI_data_offset(this, record) result(file_offset)
      !> block data file object
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> file offset
#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: file_offset
#else
      integer :: file_offset
#endif

      integer :: ierr, intsize

#ifdef MPI
      call MPI_type_size(MPI_INTEGER, intsize, ierr)
      file_offset = (record-1)*this%record_length + intsize + this%block_rank*intsize
#else
      ! to be used with MPI only
      file_offset = 0
#endif
    end function

    !> write block rank and shape for a given record
    subroutine write_rank_and_shape(this, record)
      !> block data file object
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record

#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: file_offset
#endif
      integer :: ierr, intsize
      character(:), allocatable :: errmsg

      errmsg = ''

      ! read rank and shape
#ifdef MPI
      call MPI_type_size(MPI_INTEGER, intsize, ierr)
      file_offset = (record-1)*this%record_length
      call MPI_file_write_at(this%fid, file_offset, this%block_rank, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      call mpi_error_to_string('Failed to write record rank.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(write_rank_and_shape) '//errmsg)
      file_offset = file_offset + intsize
      call MPI_file_write_at(this%fid, file_offset, this%block_shape, this%block_rank, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      call mpi_error_to_string('Failed to write record shape.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(write_rank_and_shape) '//errmsg)
#else
      write(this%fid, rec=record, iostat=ierr) this%block_rank, this%block_shape
      call terminate_if_false( ierr == 0, '(write_rank_and_shape) Failed to write record rank and shape.')
#endif
    end subroutine write_rank_and_shape

    !> read block rank and shape for a given record
    subroutine read_rank_and_shape(this, record, block_rank, block_shape)
      use precision, only: str_32, str_1024
      !> block data file object
      class(block_data_file_type), intent(inout) :: this
      !> record number / index of data block
      integer, intent(in) :: record
      !> block rank
      integer, intent(out) :: block_rank
      !> block shape
      integer, intent(out) :: block_shape(this%block_rank)

#ifdef MPI
      integer(kind=MPI_OFFSET_KIND) :: file_offset
#endif
      integer :: ierr, intsize
      character(len=str_32) :: frmt
      character(len=str_1024) :: msg
      character(:), allocatable :: errmsg

      errmsg = ''

      ! read rank and shape
#ifdef MPI
      call MPI_type_size(MPI_INTEGER, intsize, ierr)
      file_offset = (record-1)*this%record_length
      call MPI_file_read_at(this%fid, file_offset, block_rank, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      call mpi_error_to_string('Failed to read record rank.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(read_rank_and_shape) '//errmsg)
      file_offset = file_offset + intsize
      call MPI_file_read_at(this%fid, file_offset, block_shape, block_rank, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      call mpi_error_to_string('Failed to read record shape.', ierr, errmsg)
      call terminate_if_false( ierr == 0, '(read_rank_and_shape) '//errmsg)
#else
      read(this%fid, rec=record, iostat=ierr) block_rank, block_shape
      call terminate_if_false( ierr == 0, '(read_rank_and_shape) Failed to read record rank and shape.')
#endif

      ! sanity check
      write(frmt, '("(a,a,a,i6,a,",i4,"i4,""/"",",i4,"i4)")') 1, 1
      write(msg, trim(frmt)) &
        'Data blocks in file and in block_data_file_type object have different ranks.', &
        new_line('a')//' file: '//trim(this%path), &
        new_line('a')//' record: ', record, &
        new_line('a')//' rank (file object / file): ', this%block_rank, block_rank
      call terminate_if_false( this%block_rank == block_rank, &
        '(read_rank_and_shape) '//new_line('a')//trim(msg))
      
      write( frmt, '("(a,a,a,i6,a,",i4,"i9,""/"",",i4,"i9)")') this%block_rank, block_rank
      write( msg, trim(frmt)) &
        'Data blocks in file and in block_data_file_type object have different shapes.', &
        new_line('a')//' file: '//trim(this%path), &
        new_line('a')//' record: ', record, &
        new_line('a')//' shape (file object / file): ', this%block_shape, block_shape
      call terminate_if_false( all(this%block_shape == block_shape), &
        '(read_rank_and_shape) '//new_line('a')//trim(msg))
    end subroutine read_rank_and_shape

end module block_data_file
