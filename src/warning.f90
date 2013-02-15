subroutine warning(message)
    
    implicit none
    character :: message*(*)
   
    open(100,File='WARNINGS.OUT',Action='WRITE',Position='APPEND')
    write(100,*) trim(message)
    close(100)

end subroutine
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine delete_warnings()

    implicit none
    logical :: exist

    inquire (File='WARNINGS.OUT', Exist=exist)
    if (exist) Then
        open(100, File='WARNINGS.OUT')
        close(100, Status='DELETE')
    end if

end subroutine
