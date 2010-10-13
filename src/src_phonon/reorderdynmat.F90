
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine reorderdynmat(dynq,dynmatq)
  Use mod_atoms
  use mod_qpoint
  implicit none
! arguments
  Complex (8), intent(in) :: dynq (3*natmtot,3*natmtot,nqpt)
  Complex (8), intent(out) :: dynmatq(3,3,natmmax,nspecies,natmmax,nspecies,nqpt)
! local variables
  Integer :: n, iq, i, j, is, ia, ip, js, ja, jp
  n = 3 * natmtot
! store dynamical matrices in a different way
  Do iq = 1, nqpt
     i = 0
     Do is = 1, nspecies
        Do ia = 1, natoms (is)
           Do ip = 1, 3
              i = i + 1
              j = 0
              Do js = 1, nspecies
                 Do ja = 1, natoms (js)
                    Do jp = 1, 3
                       j = j + 1
                       dynmatq(ip,jp,ia,is,ja,js,iq)=dynq (i, j, iq)
                    End Do
                 End Do
              End Do
           End Do
! end loops over atoms and species
        End Do
     End Do
! end loop over q-vectors
  End Do
end subroutine
