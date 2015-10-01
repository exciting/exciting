Subroutine delinveps
      Use modmain
      use modmpi
      Implicit None
      Logical :: exist
      character(128)::sbuffer
      exist = .False.
#ifdef MPI
      write(sbuffer,*)rank
      Inquire (File='INVEPS'//trim(adjustl(sbuffer))//'.OUT', Exist=Exist)
      If (exist) Then
          open(44,file='INVEPS'//trim(adjustl(sbuffer))//'.OUT')
          Close (44, Status='DELETE')
      End if
   If (rank==0) Then   
#else
      Inquire (File='INVEPS0.OUT', Exist=Exist)
      If (exist) Then
          open(44,file='INVEPS0.OUT')
          Close (44, Status='DELETE')
      End if
#endif
      Inquire (File='INVHEAD.OUT', Exist=Exist)
      If (exist) Then
              open(42,file='INVHEAD.OUT')
              close(42,Status='DELETE')
      End If 
      Inquire (File='INVWING1.OUT', Exist=Exist)
      If (exist) Then
              open(42,file='INVWING1.OUT')
              close(42,Status='DELETE')
      End If 
      Inquire (File='INVWING2.OUT', Exist=Exist)
      If (exist) Then
              open(42,file='INVWING2.OUT')
              close(42,Status='DELETE')
      End If 
#ifdef MPI
   End If
#endif  
      Return
End Subroutine
