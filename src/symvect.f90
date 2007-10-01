
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symvect
! !INTERFACE:
subroutine symvect(vc)
! !USES:
use modmain
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
implicit none
! arguments
real(8), intent(inout) :: vc(3,natmtot)
! local variables
integer is,ia,ja,ias,jas
integer isym,lspl
real(8) s(3,3),v(3),t1
! automatic arrays
real(8) vl(3,natmtot),vs(3,natmtot)
! convert vectors to lattice coordinates
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call r3mv(ainv,vc(1,ias),vl(1,ias))
  end do
end do
! make symmetric average
vs(:,:)=0.d0
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ja=ieqatom(ia,is,isym)
      jas=idxas(ja,is)
      call r3mv(s,vl(1,jas),v)
      vs(:,ias)=vs(:,ias)+v(:)
    end do
  end do
end do
! normalise and convert to Cartesian coordinates
t1=1.d0/dble(nsymcrys)
vs(:,:)=t1*vs(:,:)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call r3mv(avec,vs(1,ias),vc(1,ias))
  end do
end do
return
end subroutine
!EOC

