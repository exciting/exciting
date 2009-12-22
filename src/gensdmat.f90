!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gensdmat
! !INTERFACE:
!
!
Subroutine gensdmat (evecsv, sdmat)
! !USES:
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   sdmat  : spin density matrices (out,complex(nspinor,nspinor,nstsv))
! !DESCRIPTION:
!   Computes the spin density matrices for a set of second-variational states.
!
! !REVISION HISTORY:
!   Created September 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Complex (8), Intent (In) :: evecsv (nstsv, nstsv)
      Complex (8), Intent (Out) :: sdmat (nspinor, nspinor, nstsv)
! local variables
      Integer :: ispn, jspn, ist, j
      Complex (8) zt1, zt2
      sdmat (:, :, :) = 0.d0
      Do j = 1, nstsv
         Do ispn = 1, nspinor
            Do jspn = 1, nspinor
               Do ist = 1, nstfv
                  zt1 = evecsv (ist+nstfv*(ispn-1), j)
                  zt2 = evecsv (ist+nstfv*(jspn-1), j)
                  sdmat (ispn, jspn, j) = sdmat (ispn, jspn, j) + zt1 * &
                 & conjg (zt2)
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
