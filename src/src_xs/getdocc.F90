!
!
!
! Copyright (C) 2010 W. Olovsson, S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getdocc (iq, ik, ikq, l1, u1, l2, u2, docc)
  ! xssave0 has to be called in advance.
      Use modinput
      Use modmain
      Use modxs
      Use m_genfilname
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik, ikq, l1, u1, l2, u2
      Real (8), Intent (Out) :: docc (u1-l1+1, u2-l2+1)
  ! local variables
      Integer :: ist, jst, iqt
      Real (8), Allocatable :: o0 (:), o (:)
      iqt = iq
      Allocate (o0(nstsv), o(nstsv))
  ! eigenvalues and occupancies for k+q-point
      Call getoccsv (vkl(1, ikq), o)
  ! eigenvalues and occupancies for k-point
      Call getoccsv0 (vkl0(1, ik), o0)
  ! loop over band ranges    
      Do ist = l1, u1
         Do jst = l2, u2
         	docc (ist-l1+1, jst-l2+1) = o0 (ist) - o (jst)
         End Do
      End Do
      Deallocate (o0, o)
End Subroutine getdocc
