!
!
!check parameterd for validity
!
!
Subroutine check_msecparameters ()
      Use modmixermsec
      Implicit None
      Logical :: usererror
      usererror = .False.
      If (qmx_input .Ge. 0.6) Then
         Write (60,*) ':WARNING: Mixing parameter may be too large and &
        &greedy'
         usererror = .True.
      Else If (qmx_input .Le. 0.025) Then
         Write (60,*) ':WARNING: Mixing parameter may be too small'
         usererror = .True.
      End If
!
!        if(noldstepsmax .lt. 1)then
!                write(60,*)':WARNING: Number of memory steps too small, using 8'
!                noldstepsmax=8
!                usererror=.true.
!        else if(noldstepsmax .lt. 4)then
!                write(60,*)':WARNING: Number of memory steps many be too small'
!                usererror=.true.
!        endif
!
End Subroutine
