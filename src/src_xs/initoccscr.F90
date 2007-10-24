
subroutine initoccscr
  use modmain
  use modtddft
  implicit none

  ! initialize number of occupied and unoccupied states
  if (nstoccscr == -1) nstoccscr=nint(chgval/2.d0)
  if (nstuoccscr == -1) nstuoccscr=nstfv-nint(chgval/2.d0)

end subroutine initoccscr
