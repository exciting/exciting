
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genffacg(is,ffacg)
use modmain
implicit none
! arguments
integer, intent(in) :: is
real(8), intent(out) :: ffacg(ngvec)
! local variables
integer ig
real(8) t1,t2,t3,t4
t1=fourpi/omega
t2=cfdamp/gmaxvr
do ig=1,ngvec
  if (gc(ig).gt.epslat) then
    if (cfdamp.ne.0.d0) then
! use damping if required
      t3=exp(-(t2*gc(ig))**2)
    else
      t3=1.d0
    end if
    t4=gc(ig)*rmt(is)
    ffacg(ig)=t1*t3*(sin(t4)-t4*cos(t4))/(gc(ig)**3)
  else
    ffacg(ig)=(t1/3.d0)*rmt(is)**3
  end if
end do
return
end subroutine

