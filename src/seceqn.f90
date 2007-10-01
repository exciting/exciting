
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine seceqn(ik,evalfv,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: evalfv(nstfv,nspnfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
do ispn=1,nspnfv
! find the matching coefficients
  call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
   sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
! solve the first-variational secular equation
  call seceqnfv(ik,ispn,apwalm(1,1,1,1,ispn),evalfv(1,ispn),evecfv(1,1,ispn))
end do
if (spinsprl) then
! solve spin-spiral secular equation
  call seceqnss(ik,apwalm,evalfv,evecfv,evecsv)
else
! solve second-variational secular equation
  call seceqnsv(ik,apwalm,evalfv,evecfv,evecsv)
end if
! compute the spin characters
call spinchar(ik,evecsv)
deallocate(apwalm)
return
end subroutine

