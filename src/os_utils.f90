module os_utils

    implicit none

    private 
    public :: make_directory_command, remove_directory_command

contains

    !> Creates a directory in the run directory (only if the directory
    !> is not yet present).
    function make_directory_command(directory_name) result(command)
        !> Name of directory
        character(*), intent(in) :: directory_name

        character(:), allocatable :: command

        command = 'test ! -e '//trim(adjustl(directory_name))//&
                    &' && mkdir '//trim(adjustl(directory_name))
    end 


   !> Creates a directory in the run directory (only if the directory
    !> is not yet present).
    function remove_directory_command(directory_name) result(command)
        !> Name of directory
        character(*), intent(in) :: directory_name

        character(:), allocatable :: command

        command = 'test ! -e '//trim(adjustl(directory_name))//&
        &' && rm -r '//trim(adjustl(directory_name))
    end 

end module
