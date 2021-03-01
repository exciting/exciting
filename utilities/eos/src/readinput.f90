
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readinput
use modmain
implicit none
! local variables
integer ipt
open(50,file='eos.in',action='READ',status='OLD',form='FORMATTED')
read(50,*) cname
read(50,*) natoms
if (natoms.le.0) then
  write(*,*)
  write(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms
  write(*,*)
  stop
end if
read(50,*) etype
read(50,*) vplt1,vplt2,nvplt
read(50,*) nevpt
if (nevpt.le.0) then
  write(*,*)
  write(*,'("Error(readinput): nevpt <= 0 : ",I8)') nevpt
  write(*,*)
  stop
end if
allocate(vpt(nevpt),ept(nevpt))
do ipt=1,nevpt
  read(50,*) vpt(ipt),ept(ipt)
end do
close(50)
return
end subroutine

