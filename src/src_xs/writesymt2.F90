
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writesymt2
! !INTERFACE:
subroutine writesymt2
! !USES:
  use modmain
  use modxs
! !DESCRIPTION:
!   Outputs the symmetrization matrices for the tensor components of a rank-2
!   tensor. The tensor(-field) in reciprocal space must be invariant under
!   coordinate transforms of the system wrt. the rotational part of the crystal
!   symmetries, so we can average:
!   $$ t_{ij}^{\rm sym} = \frac{1}{N_{\alpha}}\sum_{\alpha} \sum_{k,l}
!     \alpha_{ik}\alpha{jl}t_{kl}. $$
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! local variables
  integer :: iop1,iop2,i
  ! output the symmetrization matrices
  open(50,file='SYMT2'//trim(filext),action='WRITE',form='FORMATTED')
  write(50,*)
  write(50,'("(symmetrization matrices are in Cartesian coordinates)")')
  write(50,*)
  do iop1=1,3
     do iop2=1,3
        write(50,'("(",i1,",",i2,")-component")') iop1,iop2
        write(50,'(3f12.8)') (symt2(iop1,iop2,i,:),i=1,3)
        write(50,*)
     end do
  end do
  close(50)
end subroutine writesymt2
