module os_utils

    implicit none

    private 
    public :: create_directory, remove_directory

contains

    !> Creates a directory in the run directory (only if the directory
    !> is not yet present).
    subroutine create_directory(directory_name)
        !> Name of directory
        character(*), intent(in) :: directory_name

        character(:), allocatable :: command

        command = 'test ! -e '//trim(adjustl(directory_name))//&
                    &' && mkdir '//trim(adjustl(directory_name))
        call system(command)

    end subroutine


    !> Removes a directory in the run directory recursively 
    !>  (only if the directory is present).
    subroutine remove_directory(directory_name)
        !> Name of directory
        character(*), intent(in) :: directory_name

        character(:), allocatable :: command

        command = 'test ! -e '//trim(adjustl(directory_name))//&
                    &' && rm -r '//trim(adjustl(directory_name))
        call system(command)

    end subroutine

end module
