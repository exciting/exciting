!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: wigner3j
! !INTERFACE:
Real (8) Function wigner3j (j1, j2, j3, m1, m2, m3)
! !INPUT/OUTPUT PARAMETERS:
!   j1, j2, j3 : angular momentum quantum numbers (in,integer)
!   m1, m2, m3 : magnetic quantum numbers (in,integer)
! !DESCRIPTION:
!   Returns the Wigner $3j$-symbol. There are many equivalent definitions for
!   the $3j$-symbols, the following provides high accuracy for $j\le 50$
!   \begin{align*}
!    &\begin{pmatrix} j_1 & j_2 & j_3 \\ m_1 & m_2 & m_3 \end{pmatrix}=(-1)^
!    {j1+j2+m3}\\
!    &\times\sqrt{\frac{(j_1+m_1)!(j_2+m_2)!(j_3+m_3)!(j_3-m_3)!(j_1-m_1)!
!    (j_2-m_2)!}{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!(1+j_1+j_2+j_3)!}}
!    \times\sum_{\max(0,j_2-j_3-m_1,j_1-j_3+m_2)}^
!    {\min(j_1+j_2-j_3,j_1-m_1,j_2+m_2)}\\
!    &(-1)^k\frac{(j_2-j_1+j_3)!(j_1-j_2+j_3)!(j_1+j_2-j_3)!}
!    {(j_3-j_1-m_2+k)!(j_3-j_2+m_1+k)!(j_1+j_2-j_3-k)!k!(j_1-m_1-k)!
!    (j_2+m_2-k)}.
!   \end{align*}
!
! !REVISION HISTORY:
!   Created November 2002 (JKD)
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
! local variables
      Integer :: k, k1, k2, l1, l2, l3, n1, n2
      Real (8) :: sgn, sum, t1
! external functions
      Real (8) :: factnm, factr
      External factnm, factr
! check input variables
      If ((j1 .Lt. 0) .Or. (j2 .Lt. 0) .Or. (j3 .Lt. 0) .Or. (Abs(m1) &
     & .Gt. j1) .Or. (Abs(m2) .Gt. j2) .Or. (Abs(m3) .Gt. j3)) Then
         Write (*,*)
         Write (*, '("Error(wigner3j): invalid arguments :")')
         Write (*, '("j1 = ", I8, " j2 = ", I8, " j3 = ", I8)') j1, j2, &
        & j3
         Write (*, '("m1 = ", I8, " m2 = ", I8, " m3 = ", I8)') m1, m2, &
        & m3
         Write (*,*)
         Stop
      End If
      If ((j1 .Eq. 0) .And. (j2 .Eq. 0) .And. (j3 .Eq. 0)) Then
         wigner3j = 1.d0
         Return
      End If
      If ((j1 .Gt. 50) .Or. (j2 .Gt. 50) .Or. (j3 .Gt. 50)) Then
         Write (*,*)
         Write (*, '("Error(wigner3j): angular momenta out of range : "&
        &, 3I8)') j1, j2, j3
         Write (*,*)
         Stop
      End If
      l1 = j2 - j1 + j3
      l2 = j1 - j2 + j3
      l3 = j1 + j2 - j3
      If ((m1+m2+m3 .Ne. 0) .Or. (l1 .Lt. 0) .Or. (l2 .Lt. 0) .Or. (l3 &
     & .Lt. 0)) Then
         wigner3j = 0.d0
         Return
      End If
      n1 = j1 - m1
      n2 = j2 + m2
      k1 = Max (0, j2-j3-m1, j1-j3+m2)
      k2 = Min (l3, n1, n2)
      If (Mod(k1+j1+j2+m3, 2) .Ne. 0) Then
         sgn = - 1.d0
      Else
         sgn = 1.d0
      End If
      sum = 0.d0
      Do k = k1, k2
         t1 = sgn * factr (l1, l1-n2+k) * factr (l2, l2-n1+k) * factr &
        & (l3, l3-k)
         sum = sum + t1 / (factnm(k, 1)*factnm(n1-k, 1)*factnm(n2-k, &
        & 1))
         sgn = - sgn
      End Do
      t1 = factr (j1+m1, l1) * factr (j2+m2, l2) * factr (j3+m3, l3)
      t1 = t1 * factr (j3-m3, 1+j1+j2+j3) * factnm (j1-m1, 1) * factnm &
     & (j2-m2, 1)
      wigner3j = sum * Sqrt (t1)
      Return
End Function
!EOC
