
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: symtr2
! !INTERFACE:
subroutine symtr2(t2)
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Symmetrizes a rank-2 tensor wrt.the rotational part of the crystal
!   symmetries:
!   $$ t_{ij}^{\rm sym} = \frac{1}{N_{\alpha}}\sum_{\alpha} \sum_{k,l}
!     \alpha_{ik}\alpha{jl}t_{kl}. $$
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(inout) :: t2(3,3)
  ! local variables
  integer :: iop1,iop2,i,j
  real(8) :: s2(3,3)
  s2(:,:)=0.d0
  do iop1=1,3
     do iop2=1,3
        do i=1,3
           do j=1,3
              s2(iop1,iop2)=s2(iop1,iop2)+symt2(iop1,iop2,i,j)*t2(i,j)
           end do
        end do
     end do
  end do
  t2(:,:)=s2(:,:)
end subroutine symtr2
