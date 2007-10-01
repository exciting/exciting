
! Copyright (C) 2002-2005 S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine moke
use modmain
implicit none
! local variables
integer iw
real(8), parameter :: eps=1.d-8
complex(8) zt1,zt2
! allocatable arrays
real(8), allocatable :: w(:)
real(8), allocatable :: eps1(:)
real(8), allocatable :: eps2(:)
real(8), allocatable :: sigma1(:,:)
real(8), allocatable :: sigma2(:,:)
complex(8), allocatable :: kerr(:)
! calculate dielectric function
noptcomp=1
intraband=.true.
optcomp(:,1)=1
call linopt
intraband=.false.
optcomp(1,1)=1
optcomp(2,1)=2
call linopt
! allocate local arrays
allocate(w(nwdos))
allocate(eps1(nwdos),eps2(nwdos))
allocate(sigma1(nwdos,2),sigma2(nwdos,2))
allocate(kerr(nwdos))
! read dielectric tensor from file
open(50,file='EPSILON_11.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),eps1(iw)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),eps2(iw)
end do
close(50)
! read diagonal contribution to optical conductivity
open(50,file='SIGMA_11.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigma1(iw,1)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigma2(iw,1)
end do
close(50)
! read off-diagonal contribution to optical conductivity
open(50,file='SIGMA_12.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigma1(iw,2)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigma2(iw,2)
end do
close(50)
! calculate the MOKE
do iw=1,nwdos
  zt1=sigma1(iw,2)+zi*sigma2(iw,2)
  zt2=(sigma1(iw,1)+zi*sigma2(iw,1))*sqrt(eps1(iw)+zi*eps2(iw))
  if (abs(zt2).gt.eps) then
    kerr(iw)=-zt1/zt2
  else
    kerr(iw)=0.d0
  end if
end do
open(50,file='KERR.OUT',action='WRITE',form='FORMATTED')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),dble(kerr(iw))*180.d0/pi
end do
write(50,'("     ")')
do iw=1,nwdos
  write(50,'(2G18.10)') w(iw),aimag(kerr(iw))*180.d0/pi
end do
close(50)
write(*,*)
write(*,'("Info(moke):")')
write(*,'(" Kerr rotation in degrees written to KERR.OUT")')
write(*,*)
deallocate(eps1,eps2,sigma1,sigma2)
deallocate(w,kerr)
return
end subroutine

