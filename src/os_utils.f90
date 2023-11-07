!> Module with utility functions for simple os commands.
module os_utils
  use modmpi

  implicit none
  private

  public :: system_cmd, &
            make_directory_command, remove_directory_command, &
            make_directory, remove_directory, &
            path_exists, join_paths

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

#ifdef __INTEL_COMPILER
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
    integer, intent(out), optional :: ierr
    logical :: exists

    integer :: ierr_local
    character(:), allocatable :: dir

#ifdef __INTEL_COMPILER
    dir = trim(path)
    if( dir(len(dir):len(dir)) == '/' ) then
      inquire(directory=trim(path), exist=exists, iostat=ierr_local)
    else
      inquire(file=trim(path), exist=exists, iostat=ierr_local)
    end if
#else
    inquire(file=trim(path), exist=exists, iostat=ierr_local)
#endif
  if (present(ierr)) ierr = ierr_local
  end function path_exists

  !> Join two paths following UNIX convention:
  !>
  !> - `join_paths('path/one/', 'path/two') --> 'path/one/path/two'
  !>
  !> - `join_paths('path/one', 'path/two') --> 'path/one/path/two'
  !>
  !> - `join_paths('path/one', '/path/two') --> 'path/one/path/two'
  !>
  !> - `join_paths('path/one/', '/path/two') --> 'path/one/path/two'
  !>
  !> - `join_paths('path/one////', '///path/two') --> 'path/one/path/two'
  function join_paths(path1, path2) result(joined_path)
    !> Paths to join.
    character(*), intent(in) :: path1, path2 
    character(:), allocatable :: joined_path
    
    character(:), allocatable :: path1_cleaned, path2_cleaned

    character(1), parameter :: delimiter = '/'

    integer :: last_char_index

    ! Remove delimiter from the end of path1
    path1_cleaned = adjustl(trim(path1))
    last_char_index = len(path1_cleaned)
    if (path1_cleaned == delimiter) then
      path1_cleaned = ""
    else if (path1_cleaned(last_char_index:last_char_index) == delimiter) then
      path1_cleaned = path1_cleaned(:last_char_index - 1)
    end if

    ! Remove delimiter from the beginning of path2
    path2_cleaned = adjustl(trim(path2))
    if (path2_cleaned == delimiter) then
      path2_cleaned = ""
    else if (path2_cleaned(1:1) == delimiter) then
      path2_cleaned = path2_cleaned(2:)
    end if

    ! Join paths
    joined_path = path1_cleaned // delimiter // path2_cleaned
  end function join_paths

end module os_utils