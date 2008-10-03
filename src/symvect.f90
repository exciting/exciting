
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symvect
! !INTERFACE:
subroutine symvect(tspin,vc)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tspin : .true. if the global spin rotations should be used instead of the
!           spatial rotations (in,logical)
!   vc    : vectors in Cartesian coordinates for all atoms (in,real(3,natmtot))
! !DESCRIPTION:
!   Symmetrises a 3-vector at each atomic site by rotating and averaging over
!   equivalent atoms. Depending on {\tt tspin}, either the spatial or spin part
!   of the crystal symmetries are used.
!
! !REVISION HISTORY:
!   Created June 2004 (JKD)
!   Modified for spin rotations, May 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tspin
real(8), intent(inout) :: vc(3,natmtot)
! local variables
integer is,ia,ja,ias,jas
integer isym,lsp
real(8) v(3),t1
! automatic arrays
real(8) vs(3,natmtot)
! make symmetric average
vs(:,:)=0.d0
do isym=1,nsymcrys
  lsp=lsplsymc(isym)
  if (tspin) lsp=lspnsymc(isym)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ja=ieqatom(ia,is,isym)
      jas=idxas(ja,is)
      call r3mv(symlatc(:,:,lsp),vc(:,jas),v)
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

