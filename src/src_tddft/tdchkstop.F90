
subroutine tdchkstop
  use modtddft
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tdchkstop'
  integer :: un
  logical :: exist

  inquire(file='STOP',exist=exist)
  if (exist) then
     call getunit(un)
     write(unitout,'(a)') 'Warning('//thisnam//') stopped with message: '// &
          trim(msg)
     open(un,file='STOP')
     close(un,status='delete')
     call terminate()
  end if

end subroutine tdchkstop
