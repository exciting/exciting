subroutine nlomain

  use optica
  implicit none
  
  open(166,file='NLOINFO.OUT',status='UNKNOWN',form='FORMATTED')
  call boxmsg(166,'=',' Non-linear optical susceptibility')

! get and validate all needed input data
  call readinp
  
! calculate NLO
  call nlinopt 
  
  call boxmsg(166,'=',' End of the program NLO')
  close(166)

end subroutine
