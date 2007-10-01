
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symvect
! !INTERFACE:
subroutine symvect(vc)
! !INPUT/OUTPUT PARAMETERS:
!   vc : vectors in Cartesian coordinates for all atoms (in,real(3,natmtot))
! !DESCRIPTION:
!   Symmetrises a set of input vectors by rotating and averaging over equivalent
!   atoms. The vectors could be atomic forces for example.
!
! !REVISION HISTORY:
!   Created June 2004 (JKD)
!EOP
!BOC
use modmain
implicit none
! arguments
real(8), intent(inout) :: vc(3,natmtot)
! local variables
integer is,ia1,ia2,ias1,ias2,n,i,isym
real(8) s(3,3),v2(3),t1
! allocatable arrays
real(8), allocatable :: v1(:,:)
allocate(v1(3,natmtot))
do is=1,nspecies
  do ia1=1,natoms(is)
    ias1=idxas(ia1,is)
    v1(:,ias1)=0.d0
    n=0
    do ia2=1,natoms(is)
      ias2=idxas(ia2,is)
      do i=1,nsymeqat(ia2,ia1,is)
        isym=symeqat(i,ia2,ia1,is)
        s(:,:)=dble(symlat(:,:,isym))
! convert symmetry to Cartesian coordinates
        call r3mm(s,ainv,s)
        call r3mm(avec,s,s)
! apply to input vector
        call r3mv(s,vc(1,ias2),v2)
        v1(:,ias1)=v1(:,ias1)+v2(:)
        n=n+1
      end do
    end do
    t1=1.d0/dble(n)
    v1(:,ias1)=t1*v1(:,ias1)
  end do
! copy symmetrised vectors to input array
  do ia1=1,natoms(is)
    ias1=idxas(ia1,is)
    vc(:,ias1)=v1(:,ias1)
! underflow check
    t1=vc(1,ias1)**2+vc(2,ias1)**2+vc(3,ias1)**2
    if (t1.lt.1.d-40) vc(:,ias1)=0.d0
  end do
end do
deallocate(v1)
return
end subroutine
!EOC

