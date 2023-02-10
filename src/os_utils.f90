!> module with utility functions for simple os commands
module os_utils
  use modmpi

  implicit none
  private

  public :: system_cmd, &
            make_directory_command, remove_directory_command, &
            make_directory, remove_directory, &
            path_exists

contains

  !> Creates a directory in the run directory (only if the directory
  !> is not yet present).
  function make_directory_command(directory_name) result(command)
    !> Name of directory
    character(*), intent(in) :: directory_name

    character(:), allocatable :: command

    command = 'test ! -e '//trim(adjustl(directory_name))//&
                &' && mkdir '//trim(adjustl(directory_name))
  end function make_directory_command

  !> Removes a directory from the run directory (only if the directory
  !> is not yet present).
  function remove_directory_command(directory_name) result(command)
    !> Name of directory
    character(*), intent(in) :: directory_name

    character(:), allocatable :: command

    command = 'test ! -e '//trim(adjustl(directory_name))//&
                &' && rm -r '//trim(adjustl(directory_name))
  end function remove_directory_command

  !> execute any system command
  function system_cmd(cmd) result(ierr)
    !> command string
    character(*), intent(in) :: cmd
    !> system specific error code, 0 on success
    integer :: ierr

#ifdef IFORT
    integer, external :: system
    ierr = system(trim(cmd))
#else
    call system(trim(cmd), status=ierr)
#endif
  end function system_cmd

  !> create new directory (if not yet existent)
  function make_directory(dir, comm) result(ierr)
    !> name of / relative path to directory
    character(*), intent(in) :: dir
    !> MPI communicator
    type(mpiinfo), intent(inout) :: comm
    !> system specific error code, 0 on success
    integer :: ierr

    ierr = system_cmd('mkdir -p '//trim(adjustl(dir)))

#ifdef MPI
    call MPI_Allreduce( MPI_IN_PLACE, ierr, 1, MPI_INT, MPI_MAX, comm%comm, comm%ierr )
    comm%ierr = ierr
#endif
  end function make_directory

  !> remove directory (and all of its content)
  function remove_directory(dir, comm) result(ierr)
    !> name of / relative path to directory
    character(*), intent(in) :: dir
    !> MPI communicator
    type(mpiinfo), intent(inout) :: comm
    !> system specific error code, 0 on success
    integer :: ierr

    if (comm%rank == 0) then 
      ierr = system_cmd('rm -rf '//trim(adjustl(dir)))
      comm%ierr = ierr
    end if
    call barrier(mpicom=comm)
#ifdef MPI
    call MPI_Bcast(ierr, 1, MPI_INTEGER, 0, comm%comm, comm%ierr)
#endif
  end function remove_directory

  !> check if path exists
  function path_exists(path, ierr) result(exists)
    !> path name
    character(*), intent(in) :: path
    !> system dependent error code; `0` on success
    integer, intent(out) :: ierr
    logical :: exists

    character(:), allocatable :: dir

#ifdef IFORT
    dir = trim(path)
    if( dir(len(dir):len(dir)) == '/' ) then
      inquire(directory=trim(path), exist=exists, iostat=ierr)
    else
      inquire(file=trim(path), exist=exists, iostat=ierr)
    end if
#else
    inquire(file=trim(path), exist=exists, iostat=ierr)
#endif
  end function path_exists

end module os_utils
