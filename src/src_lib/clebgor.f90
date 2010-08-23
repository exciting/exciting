!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: clebgor
! !INTERFACE:
Real (8) Function clebgor (j1, j2, j3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Clebsch-Gordan coefficients using the Wigner $3j$-symbols
!   $$ C(J_1 J_2 J_3 | m_1 m_2 m_3)=(-1)^{J_1-J_2+m_3}\sqrt{2J_3+1}
!    \begin{pmatrix} J_1 & J_2 & J_3 \\ m_1 & m_2 & -m_3 \end{pmatrix}. $$
!   Suitable for $J_i\le 50$.
!
! !REVISION HISTORY:
!   Created September 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: j1
      Integer, Intent (In) :: j2
      Integer, Intent (In) :: j3
      Integer, Intent (In) :: m1
      Integer, Intent (In) :: m2
      Integer, Intent (In) :: m3
! external functions
      Real (8) :: wigner3j
      External wigner3j
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0) .Or. (Abs(m1) &
     & .Gt. j1) .Or. (Abs(m2) .Gt. j2) .Or. (Abs(m2) .Gt. j2) .Or. &
     & (Abs(m3) .Gt. j3)) Then
         Write (*,*)
         Write (*, '("Error(clebgor): non-physical arguments :")')
         Write (*, '("j1 = ", I8, " j2 = ", I8, " j3 = ", I8)') j1, j2, &
        & j3
         Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2, &
        & m3
         Write (*,*)
         Stop
      End If
      If ((j1 .Eq. 0) .And. (j2 .Eq. 0) .And. (j3 .Eq. 0)) Then
         clebgor = 1.d0
         Return
      End If
      If ((j1 .Gt. 50) .Or. (j2 .Gt. 50) .Or. (j3 .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(clebgor): angular momenta out of range : ",&
        & 3I8)') j1, j2, j3
         Write (*,*)
         Stop
      End If
      If ((m1+m2-m3 .Ne. 0) .Or. (j2+j3 .Lt. j1) .Or. (j1+j3 .Lt. j2) &
     & .Or. (j1+j2 .Lt. j3)) Then
         clebgor = 0.d0
         Return
      End If
      clebgor = Sqrt (dble(2*j3+1)) * wigner3j (j1, j2, j3, m1, m2,-m3)
      If (Mod(j1-j2+m3, 2) .Ne. 0) clebgor = - clebgor
      Return
End Function
!EOC
