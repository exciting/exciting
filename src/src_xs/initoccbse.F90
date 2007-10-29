
subroutine initoccbse
  use modmain
  use modtddft
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'initoccbse'
  integer :: nvalel

  ! number of valence electrons
  nvalel=nint(chgval/2.d0)
  ! initialize number of occupied states
  if (nstoccbse == -1) nstoccbse=nvalel
  if (nstoccbse > nvalel) then
     write(unitout,'(a,2i6)') 'Warning('//trim(thisnam)//'): number of &
          occupied states for BSE too large - resetting (proposed/reset):',&
          nstoccbse,nvalel
     nstoccbse=nvalel
  end if

end subroutine initoccbse
