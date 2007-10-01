
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: moment
! !INTERFACE:
subroutine moment
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total moments by integrating the
!   magnetisation.
!
! !REVISION HISTORY:
!   Created January 2005 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir,idm
real(8) sum
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
if (.not.spinpol) then
  mommt(:,:)=0.d0
  mommttot(:)=0.d0
  momir(:)=0.d0
  momtot(:)=0.d0
  return
end if
! find the muffin-tin moments
mommttot(:)=0.d0
do idm=1,ndmag
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is)
        fr(ir)=fourpi*magmt(1,ir,ias,idm)*y00*spr(ir,is)**2
      end do
      call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
      mommt(idm,ias)=gr(nrmt(is))
      mommttot(idm)=mommttot(idm)+mommt(idm,ias)
    end do
  end do
end do
! find the interstitial moments
do idm=1,ndmag
  sum=0.d0
  do ir=1,ngrtot
    sum=sum+magir(ir,idm)*cfunir(ir)
  end do
  momir(idm)=sum*omega/dble(ngrtot)
end do
momtot(:)=mommttot(:)+momir(:)
return
end subroutine
!EOC

