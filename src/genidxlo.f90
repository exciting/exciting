!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genidxlo
! !INTERFACE:
!
!
Subroutine genidxlo
! !USES:
      Use modmain
! !DESCRIPTION:
!   Generates an index array which maps the local-orbitals in each atom to their
!   locations in the overlap or Hamiltonian matrices. Also finds the total
!   number of local-orbitals.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ias, i, ilo, l, m, lm
! allocate global local-orbital index
      If (allocated(idxlo)) deallocate (idxlo)
      Allocate (idxlo(lolmmax, nlomax, natmtot))
      i = 0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ilo = 1, nlorb (is)
               l = lorbl (ilo, is)
               Do m = - l, l
                  i = i + 1
                  lm = idxlm (l, m)
                  idxlo (lm, ilo, ias) = i
               End Do
            End Do
         End Do
      End Do
      nlotot = i
      Return
End Subroutine
!EOC
