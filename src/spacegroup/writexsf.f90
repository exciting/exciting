
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writexsf
use modmain
implicit none
! local variables
integer is,ia
real(8) v(3)
open(50,file='crystal.xsf',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("CRYSTAL")')
write(50,*)
write(50,'("PRIMVEC")')
write(50,'(3G18.10)') avec(:,1)*au2ang
write(50,'(3G18.10)') avec(:,2)*au2ang
write(50,'(3G18.10)') avec(:,3)*au2ang
write(50,*)
write(50,'("PRIMCOORD")')
write(50,'(2I8)') natmtot,1
do is=1,nspecies
  do ia=1,natoms(is)
    call r3mv(avec,atposl(1,ia,is),v)
    write(50,'(A,3G18.10)') trim(spsymb(is)),v*au2ang
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writexsf):")')
write(*,'(" XCrysDen file written to crystal.xsf")')
return
end subroutine

