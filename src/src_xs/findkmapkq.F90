
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findkmapkq(iq,vq,voff,map)
  use modmain
  use modxs
  use modmpi
  implicit none
  ! arguments
  integer, intent(in) :: iq
  real(8), intent(in) :: vq(3),voff(3)
  integer, intent(out) :: map(nkpt)
  ! local variables
  integer :: ik,iv(3),ivt(3)
  real(8) :: vkq(3)
  real(8), external :: r3taxi
  character(256) :: filnam
  do ik=1,nkpt
     vkq(:)=vkl(:,ik)+vq(:)
     call r3frac(epslat,vkq,ivt)
     iv(:)=nint(vkq(:)*ngridk(:)-voff(:))
     map(ik)=ikmap(iv(1),iv(2),iv(3))
  end do
end subroutine findkmapkq
