
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geomplot
use modmain
implicit none
! local variables
integer is,ia
! Bohr to Angstroms (CODATA 2002)
real(8), parameter :: au_to_ang=0.5291772108d0
real(8) v1(3),v2(3),v3(3),v4(3),t1
real(8) dxx,dyx,dyy,dzx,dzy,dzz
! initialise universal variables
call init0
!------------------------------------------------!
!     write the XCrysden file to crystal.xsf     !
!------------------------------------------------!
open(50,file='crystal.xsf',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("CRYSTAL")')
write(50,*)
write(50,'("PRIMVEC")')
write(50,'(3G18.10)') avec(:,1)*au_to_ang
write(50,'(3G18.10)') avec(:,2)*au_to_ang
write(50,'(3G18.10)') avec(:,3)*au_to_ang
write(50,*)
write(50,'("PRIMCOORD")')
write(50,'(2I8)') natmtot,1
do is=1,nspecies
  do ia=1,natoms(is)
    call r3mv(avec,atposl(:,ia,is),v1)
    write(50,'(A,3G18.10)') trim(spsymb(is)),v1(:)*au_to_ang
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writexsf):")')
write(*,'(" XCrysDen file written to crystal.xsf")')
!-----------------------------------------------!
!     write the V_Sim file to crystal.ascii     !
!-----------------------------------------------!
! determine coordinate system vectors
t1=sqrt(avec(1,1)**2+avec(2,1)**2+avec(3,1)**2)
v1(:)=avec(:,1)/t1
t1=sqrt(avec(1,2)**2+avec(2,2)**2+avec(3,2)**2)
v2(:)=avec(:,2)/t1
call r3cross(v1,v2,v3)
t1=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
v3(:)=v3(:)/t1
call r3cross(v3,v1,v2)
t1=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
v2(:)=v2(:)/t1
dxx=dot_product(avec(:,1),v1(:))
dyx=dot_product(avec(:,2),v1(:))
dyy=dot_product(avec(:,2),v2(:))
dzx=dot_product(avec(:,3),v1(:))
dzy=dot_product(avec(:,3),v2(:))
dzz=dot_product(avec(:,3),v3(:))
open(50,file='crystal.ascii',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'(3G18.10)') dxx,dyx,dyy
write(50,'(3G18.10)') dzx,dzy,dzz
write(50,*)
do is=1,nspecies
  do ia=1,natoms(is)
    v4(1)=dot_product(atposc(:,ia,is),v1(:))
    v4(2)=dot_product(atposc(:,ia,is),v2(:))
    v4(3)=dot_product(atposc(:,ia,is),v3(:))
    write(50,'(3G18.10," ",A)') v4,trim(spsymb(is))
  end do
end do
close(50)
write(*,*)
write(*,'("Info(writevsim):")')
write(*,'(" V_Sim file written to crystal.ascii")')
return
end subroutine

