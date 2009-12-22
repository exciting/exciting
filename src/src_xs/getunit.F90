!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getunit
      Implicit None
Contains
!
!
      Subroutine getunit (un)
         Implicit None
    ! parameters
         Integer, Intent (Out) :: un
    ! local variables
         Character (*), Parameter :: thisnam = 'getunit'
         Integer :: u, u_lo, u_hi
         Logical :: connected
    ! lower value for units
         u_lo = 100
    ! upper value for units
         u_hi = 5000
         Do u = u_lo, u_hi
            Inquire (u, Opened=connected)
            If ( .Not. connected) Then
               un = u
               Return
            End If
         End Do
         Write (*, '("Error(", a, "): no free file unit available betwe&
        &en", i6, "and", i6)') thisnam, u_lo, u_hi
         Stop
      End Subroutine getunit
!
End Module m_getunit
