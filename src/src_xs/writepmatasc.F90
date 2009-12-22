!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writepmatasc
      Use modmain
      Use modxs
      Use m_getunit
      Use m_getpmat
      Implicit None
      Complex (8), Allocatable :: pmat (:, :, :)
      Integer :: un, ik, ist1, ist2, oct
  ! initialize global variables
      Call init0
      Call init1
      Call init2
      Allocate (pmat(3, nstsv, nstsv))
      Call getunit (un)
      Open (un, File='PMAT_XS_ASC.OUT', Action='write')
      Do ik = 1, nkpt
     ! read momentum matrix elements if required
         Call getpmat (ik, vkl, 1, nstsv, 1, nstsv, .True., 'PMAT_XS.OU&
        &T', pmat)
         Do ist1 = 1, nstsv
            Do ist2 = 1, nstsv
               Do oct = 1, 3
                  Write (un, '(3i8, i4, 4g18.10)') ik, ist1, ist2, oct, &
                 & pmat (oct, ist1, ist2), Abs (pmat(oct, ist1, ist2)), &
                 & Abs (pmat(oct, ist2, ist1)-conjg(pmat(oct, ist1, &
                 & ist2)))
               End Do
            End Do
         End Do
      End Do
      Close (un)
      Deallocate (pmat)
End Subroutine writepmatasc
