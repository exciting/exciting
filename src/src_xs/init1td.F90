
subroutine init1td
  use modxs, only: skipallocs1
  implicit none
  skipallocs1=.true.
  ! call init1 without (re-)allocation of radial functions
  call init1
  skipallocs1=.false.
end subroutine init1td
