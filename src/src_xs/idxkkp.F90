!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
Integer Function idxkkp (ik, ikp, n)
  use modmpi
      Implicit None
  ! arguments
      Integer, Intent (In) :: ik, ikp, n
  ! local variables
      If ((ik .Le. 0) .Or. (ikp .Le. 0) .Or. (n .Le. 0)) Then
         Write (*,*)
         Write (*, '("Error(idxkkp): negative indices or number of poin&
        &ts")')
         Write (*,*)
         Call terminate
      End If
      If (ik .Gt. ikp) Then
         Write (*,*)
         Write (*, '("Error(idxkkp): ik > ikp")')
         Write (*,*)
         Call terminate
      End If
  ! (i,j) -> (i-1)i/2 + (N-i+1)(i-1) + j-i+1
      idxkkp = - (ik*(ik-1)) / 2 + n * (ik-1) + ikp
End Function idxkkp
