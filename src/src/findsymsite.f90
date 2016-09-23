!
!
!
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findsymsite
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, js, ia, ja, ias
      Real (8) :: apl (3, maxatoms, maxspecies)
! automatic arrays
      Real (8) :: iea (natmmax, nspecies, 48)
! allocate the site symmetry arrays
      If (allocated(nsymsite)) deallocate (nsymsite)
      Allocate (nsymsite(natmtot))
      If (allocated(lsplsyms)) deallocate (lsplsyms)
      Allocate (lsplsyms(48, natmtot))
      If (allocated(lspnsyms)) deallocate (lspnsyms)
      Allocate (lspnsyms(48, natmtot))
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do js = 1, nspecies
               Do ja = 1, natoms (js)
                  apl (:, ja, js) = input%structure%speciesarray(js)%species%atomarray(ja)%atom%coord(:) - &
                 & input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
               End Do
            End Do
            Call findsym (apl, apl, nsymsite(ias), lsplsyms(:, ias), &
           & lspnsyms(:, ias), iea)
         End Do
      End Do
      Return
End Subroutine
