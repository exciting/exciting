! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
subroutine writedielt(filtag, nw, w, dt, switch)
  use mod_misc, only: filext
  use m_getunit

  implicit none

  ! arguments
  character(*), intent(in) :: filtag
  integer, intent(in) :: nw
  real(8), intent(in) :: w(nw)
  complex(8), intent(in) :: dt(3, 3, nw)
  integer, intent(in) :: switch

  ! local variables
  integer :: un, oct, iw

  call getunit(un)
  open(un, file=trim(filtag)//trim(filext), form='formatted',&
    & action='write', status='replace')
  write(un,*)

  if(switch .eq. 0) then
    write(un, '(" (dielectric tensor, independent particle approximation)")')
  else
    write(un, '(" (dielectric tensor, random phase approximation)")')
  end if
  write(un,*)

  do iw = 1, nw
    write(un, '(" frequency index and value: ",i6,f14.8)') iw, w(iw)
    write(un, '(" real part, imaginary part below")')
    write(un, '(3f14.8,5x,3f14.8)')&
      & (dble(dt(oct, :, iw)), aimag(dt(oct, :, iw)), oct=1, 3)
    write(un,*)
  end do

  close(un)
end subroutine writedielt
