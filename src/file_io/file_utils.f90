module file_utils
  use mod_misc, only: filext
  use modmpi
  use precision, only: i32, str_256
  use exciting_mpi, only: xmpi_bcast

  implicit none
  private

  public :: add_default_extension, &
            close_file, &
            copy_text_file, &
            delete_file, &
            file_is_open, &
            read_last_and_penultimate_lines_from_file

  interface file_is_open
    procedure :: file_is_open_by_name, file_is_open_by_unit
  end interface file_is_open

  interface close_file
    procedure :: close_file_serial, close_file_parallel
  end interface close_file

  interface delete_file
    procedure :: delete_file_serial, delete_file_parallel
  end interface delete_file

  contains

    !> close a file for serial i/o
    subroutine close_file_serial(fid, ierr)
      !> i/o unit
      integer, intent(in) :: fid
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr

      ! return if file is not open
      if (.not. file_is_open(fid, ierr)) return

      ! close file
      close(fid, iostat=ierr)
    end subroutine close_file_serial

    !> close a file for parallel i/o
    subroutine close_file_parallel(fid, comm, ierr)
      !> MPI file handler
      integer, intent(inout) :: fid
      !> MPI communicator
      type(mpiinfo), intent(inout) :: comm
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr

#ifndef MPI
      ! fall back to serial routine if MPI is not in use
      call close_file_serial(fid, ierr)
      comm%ierr = ierr
      return
#else
      ! ommitting the check if file is open here, because I
      ! didn't find out how to check with MPI
      ! feel free do add

      ! close file
      call MPI_file_close(fid, ierr)
      comm%ierr = ierr
#endif
    end subroutine close_file_parallel

    !> delete file for serial i/o
    subroutine delete_file_serial(fname, ierr)
      use os_utils, only: path_exists
      !> file name
      character(*), intent(in) :: fname
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr

      integer :: id

      ! return if file does not exist
      if (.not. path_exists(fname, ierr)) return

      ! open file if it isn't
      if (.not. file_is_open_by_name(fname, id, ierr)) then
        open(newunit=id, file=trim(fname), status='old', iostat=ierr)
      end if
      if (ierr /= 0) return

      ! delete file
      close(id, status='delete', iostat=ierr)
    end subroutine delete_file_serial

    !> delete file for parallel i/o
    subroutine delete_file_parallel(fname, comm, ierr)
      use os_utils, only: path_exists
      !> file name
      character(*), intent(in) :: fname
      !> MPI communicator
      type(mpiinfo), intent(inout) :: comm
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr

      logical :: exists, all_exists

#ifndef MPI
      ! fall back to serial routine if MPI is not in use
      call delete_file_serial(fname, ierr)
      comm%ierr = ierr
      return
#else
      ! return if file does not exist
      exists = path_exists(fname, comm%ierr)
      call MPI_Allreduce(exists, all_exists, 1, MPI_LOGICAL, MPI_LAND, comm%comm, comm%ierr)
      call terminate_if_false( all_exists .eqv. exists, '(delete_file_parallel) &
        File exists for one but nor for all processes.' )
      if( .not. exists ) return

      ! delete file
      if (comm%rank == 0) call MPI_file_delete(trim(fname), MPI_INFO_NULL, ierr)
      call xmpi_bcast( comm, ierr )
#endif
    end subroutine delete_file_parallel

    !> check if file is open by file name
    function file_is_open_by_name(fname, fid, ierr) result(isopen)
      !> file name
      character(*), intent(in) :: fname
      !> i/o unit of file; `-1` if file is not open
      integer, intent(out) :: fid
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr
      logical :: isopen

      isopen = .false.

      ! check if file is connected to a unit
      inquire(file=trim(fname), number=fid, iostat=ierr)

      ! check if unit is open
      if (fid >= 0 .and. ierr == 0) &
        isopen = file_is_open_by_unit(fid, ierr)
    end function file_is_open_by_name

    !> check if file is open by unit
    function file_is_open_by_unit(fid, ierr) result(isopen)
      !> i/o unit of file; `-1` if file is not open
      integer, intent(in) :: fid
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr
      logical :: isopen

      inquire(unit=fid, opened=isopen, iostat=ierr)
    end function file_is_open_by_unit

    !> Copy a file, from source to destination, line by line. It only works for text files
    subroutine copy_text_file( source_name, destination_name )
      !> Name of the source file
      character(len=*), intent(in) :: source_name
      !> Name of the destination file
      character(len=*), intent(in) :: destination_name

      integer(i32) :: unit_source, unit_dest, ios
      character(len=str_256) :: line
      logical :: file_exists

      inquire( file=trim(source_name), exist=file_exists )
      call terminate_if_false( file_exists, "Error: File " // trim(source_name) // " not found" )
      open( newunit=unit_source, file=source_name, status="old", position="rewind", action="read" )
      open( newunit=unit_dest, file=destination_name, status="replace", position="rewind", action="write" )
      do 
        read( unit_source, '(A)', iostat=ios ) line
        if ( ios /= 0 ) exit
        write( unit_dest, '(A)' ) trim( line )
      end do
      close( unit_source )
      close( unit_dest )
    end subroutine

    !> Read a file and store the content of the last and penultiname lines
    subroutine read_last_and_penultimate_lines_from_file( file_name, last_line, penultimate_line )
      !> File name
      character(len=*), intent(in) :: file_name
      !> Last line
      character(len=*), intent(out) :: last_line
      !> Penultimate line
      character(len=*), intent(out) :: penultimate_line
    
      character(len=str_256) :: line
      integer(i32) :: unit, ios
      logical :: file_exists
    
      last_line = "empty"; penultimate_line = "empty"
      inquire( file=trim( file_name ), exist=file_exists )
      call terminate_if_false( file_exists, "Error: File " // trim(file_name) // " not found" )
      open( newunit=unit, file=file_name, status="old", action="read")
      do 
        read( unit, '(A)', iostat=ios ) line
        if ( ios /= 0 ) exit
        penultimate_line = last_line
        last_line = line
      end do
      if( trim(penultimate_line) == "empty" ) call terminate( "Error: file " // file_name // " has less than two lines" )
      close( unit )
    end subroutine  

    !> Add the default extension `filext` from the [[mod_misc]] (usually `.OUT`) to a base file name
    pure function add_default_extension( file_name ) result( name )
      !> base file name
      character(len=*), intent(in)  :: file_name
      !> file name with default extension
      character(len=:), allocatable :: name

      name = trim( file_name )//trim( filext )
    end function

end module file_utils
