!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gauntyry
! !INTERFACE:
Complex (8) Function gauntyry (l1, l2, l3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
!   l1, l2, l3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the complex Gaunt-like coefficient given by
!   $\langle Y^{l_1}_{m_1}|R^{l_2}_{m_2}|Y^{l_3}_{m_3}\rangle$, where $Y_{lm}$
!   and $R_{lm}$ are the complex and real spherical harmonics, respectively.
!   Suitable for $l_i$ less than 50. See routine {\tt genrlm}.
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: l1
      Integer, Intent (In) :: l2
      Integer, Intent (In) :: l3
      Integer, Intent (In) :: m1
      Integer, Intent (In) :: m2
      Integer, Intent (In) :: m3
! local variables
! real constant sqrt(2)/2
      Real (8), Parameter :: c1 = 0.7071067811865475244d0
      Real (8) :: t1, t2
! external functions
      Real (8) :: gaunt
      External gaunt
      If (m2 .Gt. 0) Then
         If (Mod(m2, 2) .Eq. 0) Then
            t1 = 1.d0
         Else
            t1 = - 1.d0
         End If
         t2 = c1 * (gaunt(l1, l2, l3, m1, m2, m3)+t1*gaunt(l1, l2, l3, &
        & m1,-m2, m3))
         gauntyry = cmplx (t2, 0.d0, 8)
      Else If (m2 .Lt. 0) Then
         If (Mod(m2, 2) .Eq. 0) Then
            t1 = 1.d0
         Else
            t1 = - 1.d0
         End If
         t2 = c1 * (gaunt(l1, l2, l3, m1, m2, m3)-t1*gaunt(l1, l2, l3, &
        & m1,-m2, m3))
         gauntyry = cmplx (0.d0,-t2, 8)
      Else
         gauntyry = cmplx (gaunt(l1, l2, l3, m1, m2, m3), 0.d0, 8)
      End If
      Return
End Function
!EOC
