!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: factnm
! !INTERFACE:
Real (8) Function factnm (n, m)
! !INPUT/OUTPUT PARAMETERS:
!   n : input (in,integer)
!   m : order of multifactorial (in,integer)
! !DESCRIPTION:
!   Returns the multifactorial
!   $$ n\underbrace{!!\,...\,!}_{m\,{\rm times}}=\prod_{i\ge 0,\,n-im>0}n-im $$
!   for $n,\,m \ge 0$. $n$ should be less than 150.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: m
! local variables
      Integer :: i, j
      Real (8) :: f1 (24), f2 (38)
      Data f1 / 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0, 5040.d0, &
     & 40320.d0, 362880.d0, 3628800.d0, 39916800.d0, 479001600.d0, &
     & 6227020800.d0, 87178291200.d0, 1307674368000.d0, &
     & 20922789888000.d0, 355687428096000.d0, 6402373705728000.d0, &
     & 121645100408832000.d0, 2432902008176640000.d0, &
     & 51090942171709440000.d0, 1124000727777607680000.d0, &
     & 25852016738884976640000.d0, 620448401733239439360000.d0 /
      Data f2 / 1.d0, 2.d0, 3.d0, 8.d0, 15.d0, 48.d0, 105.d0, 384.d0, &
     & 945.d0, 3840.d0, 10395.d0, 46080.d0, 135135.d0, 645120.d0, &
     & 2027025.d0, 10321920.d0, 34459425.d0, 185794560.d0, &
     & 654729075.d0, 3715891200.d0, 13749310575.d0, 81749606400.d0, &
     & 316234143225.d0, 1961990553600.d0, 7905853580625.d0, &
     & 51011754393600.d0, 213458046676875.d0, 1428329123020800.d0, &
     & 6190283353629375.d0, 42849873690624000.d0, &
     & 191898783962510625.d0, 1371195958099968000.d0, &
     & 6332659870762850625.d0, 46620662575398912000.d0, &
     & 221643095476699771875.d0, 1678343852714360832000.d0, &
     & 8200794532637891559375.d0, 63777066403145711616000.d0 /
! fast return if possible
      If (n .Eq. 0) Then
         factnm = 1.d0
         Return
      End If
      If (m .Eq. 1) Then
         If ((n .Ge. 1) .And. (n .Le. 24)) Then
            factnm = f1 (n)
            Return
         End If
      End If
      If (m .Eq. 2) Then
         If ((n .Ge. 1) .And. (n .Le. 38)) Then
            factnm = f2 (n)
            Return
         End If
      End If
      If (n .Lt. 0) Then
         Write (*,*)
         Write (*, '("Error(factnm): n < 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (m .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(factnm): m <= 0 : ", I8)') m
         Write (*,*)
         Stop
      End If
      If (n .Gt. 150) Then
         Write (*,*)
         Write (*, '("Error(factnm): n out of range : ", I8)') n
         Write (*,*)
         Stop
      End If
      If (m .Eq. 1) Then
         factnm = f1 (24)
         Do i = 25, n
            factnm = factnm * dble (i)
         End Do
      Else
         j = n / m
         If (Mod(n, m) .Eq. 0) j = j - 1
         factnm = dble (n)
         Do i = 1, j
            factnm = factnm * dble (n-i*m)
         End Do
      End If
      Return
End Function
!EOC
