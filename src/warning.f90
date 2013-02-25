subroutine warning(message)
use modmpi
    
    implicit none
    character :: message*(*)
   
    if (rank==0) then
      open(100,File='WARNINGS.OUT',Action='WRITE',Position='APPEND')
      write(100,*) trim(message)
      close(100)
    end if

end subroutine
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine delete_warnings()
use modmpi

    implicit none
    logical :: exist
    
    if (rank==0) then
      inquire (File='WARNINGS.OUT', Exist=exist)
      if (exist) Then
        open(100, File='WARNINGS.OUT')
        close(100, Status='DELETE')
      end if
    end if

end subroutine
