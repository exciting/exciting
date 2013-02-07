
Subroutine getevecfvgw(ik, evecfv)

      Use modmain
      Use modinput
      Use modmpi

      Implicit None
  ! arguments
      Integer, Intent (In) :: ik
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv, nspnfv)
  ! local variables
      Integer(8) :: recl 
      Integer(4) :: nmatmax_, nstfv_, nspnfv_
      Real(8) :: vkl_(3)
      Character (256) :: filetag
      Character (256), External :: outfilenamestring
      Logical :: Exist

      Inquire (IoLength=Recl) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
      filetag = trim (filetag_evecfv)
      Inquire (File=outfilenamestring(filetag, ik), Exist=Exist)
      If (exist) Then
        Open (70, File=outfilenamestring(filetag, ik), Action='READ&
       &', Form='UNFORMATTED', Access='DIRECT', Recl=Recl)
      Else
        Write (*,*) '(getevecfvgw): Error when reading",outfilenamestring(filetag, ik)," file!'
        Stop
      End If
      Read (70, Rec=ik) vkl_, nmatmax_, nstfv_, nspnfv_, evecfv
      Close (70)
      
Return
End Subroutine
!EOC
