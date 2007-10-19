
module m_filedel
  implicit none
contains

  subroutine filedel(fnam)
    use m_getunit
    implicit none
    ! arguments
    character(*), intent(in) :: fnam
    ! local variables
    character(*), parameter :: thisnam = 'filedel'
    integer :: fu
    logical :: existent, opened

    ! check if file exists
    inquire(file=trim(fnam),exist=existent)
    if (.not.existent) then
!!$       write(*,*) 'Warning('//thisnam//'): attempted to delete &
!!$            non-existent file:'//trim(fnam)
       return
    end if
    ! check if file is opened
    inquire(file=trim(fnam),opened=opened,number=fu)
    if (opened) then
       close(fu)
    end if
    call getunit(fu)
    open(fu,file=trim(fnam),action='write')
    ! delete file
    close(fu,status='delete')

  end subroutine filedel

end module m_filedel
