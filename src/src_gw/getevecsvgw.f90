
Subroutine getevecsvgw(ik, evecsv)

      Use modmain
      Use modinput
      Use modmpi
      Implicit None

! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: recl, nstsv_
      Real (8) :: vkl_ (3)
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
      Logical :: exist

      Inquire (IoLength=Recl) vkl_, nstsv_, evecsv
      filetag = trim (filetag_evecsv)
      Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
      If (exist) Then
        Open (70, File=outfilenamestring(filetag, ik), Action='READ', &
       &  Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Else
        Write (*,*) '(getevecsvgw): Error when reading EVECSV.OUT file!'
        Stop
      End If
      Read (70, Rec=ik) vkl_, nstsv_, evecsv
      Close (70)
      
Return
End Subroutine
