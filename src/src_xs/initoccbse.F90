
subroutine initoccbse(nempty_)
  use modmain
  use modxs

  implicit none
  ! arguments
  integer, intent(inout) :: nempty_
  ! local variables
  character(*), parameter :: thisnam = 'initoccbse'
  integer :: nemptyt, nvalel

  ! number of valence electrons
  nvalel=nint(chgval/2.d0)

  ! number of states below Fermi energy
  if (nstbef == -1) nstbef=nvalel
  if (nstbef > nvalel) then
     write(unitout,'(a,2i6)') 'Error('//trim(thisnam)//'): number of &
          &states below Fermi energy for BSE too large (proposed/max):',&
          nstbef,nvalel
     call terminate
  end if

  ! number of states above Fermi energy
  if (nstabf==-1) then
     ! if "nstabf" is not specified define it using "nempty"
     nstabf=nempty_+1
  else
     ! replace "nempty" using "nstabf"
     nemptyt=nempty_
     nempty_=nstabf-1
     write(unitout,'("Info(init0): nempty has been adjusted from ",I6," to ",&
          &I6)') nemptyt,nempty_
  end if

end subroutine initoccbse
