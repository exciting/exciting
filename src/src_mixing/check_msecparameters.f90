

!check parameterd for validity


subroutine check_msecparameters()
 use  modmixermsec
 implicit none
 logical::usererror
 usererror=.false.
 if(qmx_input .ge. 0.6)then
		write(60, *)':WARNING: Mixing parameter may be too large and greedy'
		usererror=.true.
	else if(qmx_input .le. 0.025)then
		write(60, *)':WARNING: Mixing parameter may be too small'
		usererror=.true.
	endif
!
!        if(noldstepsmax .lt. 1)then
!                write(60,*)':WARNING: Number of memory steps too small, using 8'
!                noldstepsmax=8
!                usererror=.true.
!        else if(noldstepsmax .lt. 4)then
!                write(60,*)':WARNING: Number of memory steps many be too small'
!                usererror=.true.
!        endif

 end subroutine
