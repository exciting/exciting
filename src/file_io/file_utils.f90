module file_utils
  use modmpi

  implicit none
  private

  public :: file_is_open, close_file, delete_file

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
      use m_getunit
      !> file name
      character(*), intent(in) :: fname
      !> system dependent error code; `0` on success
      integer, intent(out) :: ierr

      integer :: id

      ! return if file does not exist
      if (.not. path_exists(fname, ierr)) return

      ! open file if it isn't
      if (.not. file_is_open_by_name(fname, id, ierr)) then
        call getunit(id)
        open(unit=id, file=trim(fname), status='old', iostat=ierr)
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
      call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, comm%comm, comm%ierr)
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

end module file_utils
