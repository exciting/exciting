!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: gaunt
! !INTERFACE:
Real (8) Function gaunt (l1, l2, l3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
!   l1, l2, l3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Gaunt coefficient given by
!   $$  \langle Y^{l_1}_{m_1}|Y^{l_2}_{m_2}|Y^{l_3}_{m_3} \rangle
!    = (-1)^{m_1}\left[\frac{(2l_1+1)(2l_2+1)(2l_3+1)}{4\pi} \right]
!    ^{\frac{1}{2}}
!    \begin{pmatrix} l_1 & l_2 & l_3 \\  0   & 0   & 0   \end{pmatrix}
!    \begin{pmatrix} l_1 & l_2 & l_3 \\ -m_1 & m_2 & m_3 \end{pmatrix}. $$
!   Suitable for $l_i$ less than 50.
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
      Integer :: j, j1, j2, j3, jh
      Real (8) :: t1
! real constant 1/sqrt(4*pi)
      Real (8), Parameter :: c1 = 0.28209479177387814347d0
! external functions
      Real (8) :: wigner3j, factr, factnm
      External wigner3j, factr, factnm
      If ((l1 .Lt. 0) .Or. (l2 .Lt. 0) .Or. (l3 .Lt. 0) .Or. (Abs(m1) &
     & .Gt. l1) .Or. (Abs(m2) .Gt. l2) .Or. (Abs(m3) .Gt. l3)) Then
         Write (*,*)
         Write (*, '("Error(gaunt): non-physical arguments :")')
         Write (*, '("l1 = ", I8, " l2 = ", I8, " l3 = ", I8)') l1, l2, &
        & l3
         Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2, &
        & m3
         Write (*,*)
         Stop
      End If
      If ((l1 .Gt. 50) .Or. (l2 .Gt. 50) .Or. (l3 .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(gaunt): angular momenta out of range : ", 3&
        &I8)') l1, l2, l3
         Write (*,*)
         Stop
      End If
      If (m1-m2-m3 .Ne. 0) Then
         gaunt = 0.d0
         Return
      End If
      j1 = l2 - l1 + l3
      j2 = l1 - l2 + l3
      j3 = l1 + l2 - l3
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0)) Then
         gaunt = 0.d0
         Return
      End If
      j = l1 + l2 + l3
      If (Mod(j, 2) .Ne. 0) Then
         gaunt = 0.d0
         Return
      End If
      jh = j / 2
      t1 = Sqrt (dble((2*l1+1)*(2*l2+1)*(2*l3+1))*factr(j1, &
     & j+1)*factnm(j2, 1)*factnm(j3, 1))
      t1 = t1 * factr (jh, jh-l1) / (factnm(jh-l2, 1)*factnm(jh-l3, 1))
      gaunt = t1 * c1 * wigner3j (l1, l2, l3,-m1, m2, m3)
      If (Mod(m1+jh, 2) .Ne. 0) gaunt = - gaunt
      Return
End Function
!EOC
