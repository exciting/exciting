!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gensfacgp
! !INTERFACE:
!
!
Subroutine gensfacgp (ngp, vgpc, ld, sfacgp)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,*))
!   ld     : leading dimension (in,integer)
!   sfacgp : structure factors of G+p-vectors (out,complex(ld,natmtot))
! !DESCRIPTION:
!   Generates the atomic structure factors for a set of ${\bf G+p}$-vectors:
!   $$ S_{\alpha}({\bf G+p})=\exp(i({\bf G+p})\cdot{\bf r}_{\alpha}), $$
!   where ${\bf r}_{\alpha}$ is the position of atom $\alpha$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ngp
      Real (8), Intent (In) :: vgpc (3, ngp)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: sfacgp (ld, natmtot)
! local variables
      Integer :: is, ia, ias, igp
      Real (8) :: t1
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do igp = 1, ngp
               t1 = dot_product (vgpc(:, igp), atposc(:, ia, is))
               sfacgp (igp, ias) = cmplx (Cos(t1), Sin(t1), 8)
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
