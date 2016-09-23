!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writekmapkq (iq)
      Use modmain
      Use modxs
      Use m_getunit
      Implicit None
      Integer, Intent (In) :: iq
      Integer :: un, ik
      Call getunit (un)
      Open (un, File='KMAPKQ'//trim(filext), Form='formatted', Action='&
     &write', Status='replace')
      Write (un, '(i9, a)') nkpt, ' : nkpt; k-point, ikmapikq below'
      Do ik = 1, nkpt
         Write (un, '(2i9)') ik, ikmapikq (ik, iq)
      End Do
      Close (un)
End Subroutine writekmapkq
