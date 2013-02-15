subroutine warning(message)
    
    character, optional :: message*(*)
    logical :: exist
   
    if (message.ne.'') then
       open(100,File='WARNINGS.OUT',Action='WRITE',Position='APPEND')
       write(100,*) trim(message)
       close(100)
    else
        inquire (File='WARNINGS.OUT', Exist=exist)
        if (exist) Then
            open(100, File='WARNINGS.OUT')
            close(100, Status='DELETE')
        end if   
   end if

end
