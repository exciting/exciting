
Subroutine getevecfvgw(ik, evecfv)

      Use modmain
      Use modinput
      Use modmpi

      Implicit None
  ! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv, nspnfv)
  ! local variables
      Integer(4) :: recl 
      Integer(4) :: nmatmax_, nstfv_, nspnfv_
      Real(8) :: vkl_(3)
      character(256) :: filename
      logical :: exist

      if (input%gw%skipgnd) then 
         filename = "EVECFV.OUT"
      else     
         filename = "EVECFV_GW.OUT"
      end if 
      
      Inquire (File=trim(filename), Exist=exist)
      If (exist) Then
        Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
        Open (70, File=trim(filename), Action='READ', &
        &  Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Else
        Write (*,*) '(getevecfvgw): Cannot open EVECFV_GW.OUT file!'
        Stop
      End If
      Read(70, Rec=ik) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
      Close(70)
      
      Return
End Subroutine

