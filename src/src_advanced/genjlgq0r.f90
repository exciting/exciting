
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genjlgq0r(gq0,jlgq0r)
use modmain
implicit none
! arguments
real(8), intent(in) :: gq0
real(8), intent(out) :: jlgq0r(0:lmaxvr,nrcmtmax,nspecies)
! local variables
integer is,irc
real(8) t1
do is=1,nspecies
  do irc=1,nrcmt(is)
    t1=gq0*rcmt(irc,is)
    call sbessel(lmaxvr,t1,jlgq0r(:,irc,is))
  end do
end do
return
end subroutine

