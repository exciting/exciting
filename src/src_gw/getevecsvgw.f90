
Subroutine getevecsvgw(ik, evecsv)

      Use modmain
      Use modinput
      Use modmpi
      Implicit None

! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
! local variables
      Integer :: nstsv_
      Integer(4) :: recl
      Real (8) :: vkl_ (3)
      character(256) :: filename
      Logical :: exist

      if (input%gw%skipgnd) then 
         filename = "EVECSV.OUT"
      else     
         filename = "EVECSV_GW.OUT"
      end if 
      
      Inquire (File=trim(filename), Exist=Exist)
      If (exist) Then
        Inquire (IoLength=Recl) vkl_, nstsv_, evecsv
        Open (70, File=trim(filename), Action='READ', &
       &  Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Else
        Write (*,*) '(getevecsvgw): Cannot open EVECSV_GW.OUT file!'
        Stop
      End If
      Read (70, Rec=ik) vkl_, nstsv_, evecsv
      Close (70)
      
Return
End Subroutine
