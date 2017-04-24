! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gensfacgp
! !INTERFACE:
!
!
Subroutine gensfacgp(ngp, vgpc, ld, sfacgp)
! !USES:
      Use mod_atoms, only: nspecies, natoms, natmtot, idxas, atposc
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
!   Added OMP, removed modmpi 2017 (BA)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: ngp
      Real (8), Intent (In) :: vgpc (3, ngp)
      Integer, Intent (In) :: ld
      Complex (8), Intent (Out) :: sfacgp (ld, natmtot)
! local variables
      Integer :: is, ia, igp
      Real (8) :: t1
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igp,t1,is,ia)
#endif    
      Do is = 1, nspecies
         Do ia = 1, natoms(is)
#ifdef USEOMP
!$OMP DO 
#endif    
            Do igp = 1, ngp
               t1 = dot_product (vgpc(:, igp), atposc(:, ia, is))
               sfacgp (igp, idxas(ia, is)) = cmplx (Cos(t1), Sin(t1), 8)
            End Do
#ifdef USEOMP
!$OMP END DO
#endif    
         End Do
      End Do
#ifdef USEOMP
!$OMP END PARALLEL
#endif    
      Return
End Subroutine
!EOC
