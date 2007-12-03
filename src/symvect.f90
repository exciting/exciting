
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
real(8) v(3),t1
! automatic arrays
real(8) vs(3,natmtot)
! make symmetric average
vs(:,:)=0.d0
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ja=ieqatom(ia,is,isym)
      jas=idxas(ja,is)
      call r3mv(symlatc(1,1,lspl),vc(1,jas),v)
      vs(:,ias)=vs(:,ias)+v(:)
    end do
  end do
end do
! normalise
t1=1.d0/dble(nsymcrys)
vc(:,:)=t1*vs(:,:)
return
end subroutine
!EOC

