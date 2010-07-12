
! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writesymt2
! !INTERFACE:
Subroutine writesymt2
! !USES:
      Use modmain
      Use modxs
! !DESCRIPTION:
!   Outputs the symmetrization matrices for the tensor components of a rank-2
!   tensor. The tensor(-field) $t_{ij}$ in reciprocal space must be invariant under
!   coordinate transforms of the system wrt. the rotational part of the crystal
!   symmetries, so we can average:
!   $$ t_{ij}^{\rm sym} = \frac{1}{N_{\alpha}}\sum_{\alpha} \sum_{k,l}
!     \alpha_{ik}\alpha_{jl}t_{kl}. $$
!   The symmetrized tensor $t_{ij}^{\rm sym}$ can then be written as
!   $$ t_{ij}^{\rm sym} = \sum_{k,l} T_{ij,kl} t_{kl}, $$
!   with the symmetrization tensor
!   $$ T_{ij,kl} = \frac{1}{N_{\alpha}}\sum_{\alpha}\alpha_{ik}\alpha_{jl} $$
!   where $N_{\alpha}$ is the number of symmetry operations in the space group.
!   For each component $ij$ the symmetrization tensor $T_{ij,kl}$ is written as
!   a matrix in the components $kl$ to the file {\tt SYMT2.OUT}.
!
! !REVISION HISTORY:
!   Created October 2008 (Sagmeister)
!EOP
!BOC
      Implicit None
  ! local variables
      real(8), parameter :: eps=1.d-8
      character(256) :: str
      Integer :: iop1, iop2, i
      real(8) :: t1
  ! output the symmetrization matrices
      Open (50, File='SYMT2'//trim(filext), Action='WRITE', Form='FORMATTED')
      Write (50,*)
      Write (50, '("(symmetrization matrices are in Cartesian coordinates)")')
      Write (50,*)
      Do iop1 = 1, 3
         Do iop2 = 1, 3
            t1=sum(abs(symt2(iop1,iop2,:,:)))
            str=""
            if (t1 .lt. eps) str='," ; zero contribution"'
            Write (50, '("(", i1, ", ", i2, ")-component"'//trim(str)//')') iop1, iop2
            Write (50, '(3f12.8)') (symt2(iop1, iop2, i, :), i=1, 3)
            Write (50,*)
         End Do
      End Do
      Close (50)
End Subroutine writesymt2
