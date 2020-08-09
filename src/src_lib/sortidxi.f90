!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!BOP
! !ROUTINE: sortidx
! !INTERFACE:
!
!
Subroutine sortidxi (n, a, idx)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of elements in array (in,integer)
!   idx : permutation index (out,integer(n))
!   a   : integer array (in,integer(n))
! !DESCRIPTION:
!   Integer version of sortidx.
!
! !REVISION HISTORY:
!   Created August 2020 (SeTi)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Integer, Intent (In) :: a (n)
      Integer, Intent (Out) :: idx (n)
! local variables
      Integer :: i, j, k, l, m
      If (n .Le. 0) Then
         Write (*,*)
         Write (*, '("Error(sortidx): n <= 0 : ", I8)') n
         Write (*,*)
         Stop
      End If
      Do i = 1, n
         idx (i) = i
      End Do
      If (n .Eq. 1) Return
      l = n / 2 + 1
      k = n
10    Continue
      If (l .Gt. 1) Then
         l = l - 1
         m = idx (l)
      Else
         m = idx (k)
         idx (k) = idx (1)
         k = k - 1
         If (k .Eq. 1) Then
            idx (1) = m
            Return
         End If
      End If
      i = l
      j = l + l
20    Continue
      If (j .Le. k) Then
         If (j .Lt. k) Then
            If (a(idx(j)) .Lt. a(idx(j+1))) j = j + 1
         End If
         If (a(m) .Lt. a(idx(j))) Then
            idx (i) = idx (j)
            i = j
            j = j + j
         Else
            j = k + 1
         End If
         Go To 20
      End If
      idx (i) = m
      Go To 10
End Subroutine
!EOC
