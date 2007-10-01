
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
real(8), allocatable :: sigint1(:)
real(8), allocatable :: sigint2(:)
real(8), allocatable :: epsint1(:)
real(8), allocatable :: epsint2(:)
complex(8), allocatable :: kerr(:)
noptcomp=2
optcomp(:,1)=1
optcomp(1,2)=1
optcomp(2,2)=2
optcomp(3,2)=1
! calculate dielectric function
call linopt
allocate(w(nwdos))
allocate(eps1(nwdos),eps2(nwdos))
allocate(epsint1(nwdos),epsint2(nwdos))
allocate(sigma1(nwdos,2),sigma2(nwdos,2))
allocate(sigint1(nwdos),sigint2(nwdos))
allocate(kerr(nwdos))
! read interband contribution to dielectric tensor
open(50,file='EPSILON_11.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),eps1(iw)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),eps2(iw)
end do
close(50)
! read intraband contribution to dielectric tensor
open(50,file='EPSINTRA_11.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),epsint1(iw)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),epsint2(iw)
end do
close(50)
! read diagonal interband contribution to optical conductivity
open(50,file='SIGMA_11.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigma1(iw,1)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigma2(iw,1)
end do
close(50)
! read diagonal intraband contribution to optical conductivity
open(50,file='SIGINTRA_11.OUT',action='READ',form='FORMATTED')
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigint1(iw)
end do
read(50,*)
do iw=1,nwdos
  read(50,'(2G18.10)') w(iw),sigint2(iw)
end do
close(50)
! read off-diagonal interband contribution to optical conductivity
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
  zt2=((sigma1(iw,1)+sigint1(iw))+zi*(sigma2(iw,1)+sigint2(iw)))* &
   sqrt((eps1(iw)+epsint1(iw))+zi*(eps2(iw)+epsint2(iw)))
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
deallocate(eps1,eps2,epsint1,epsint2)
deallocate(sigma1,sigma2,sigint1,sigint2)
deallocate(w,kerr)
return
end subroutine

