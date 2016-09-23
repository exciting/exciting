!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine terminateqry (str)
      Use m_getunit
      Implicit None
  ! arguments
      Character (*) :: str
  ! local variables
      Logical :: exis
      Integer :: un
  ! check for terminating
      Inquire (File='TERMINATE', Exist=exis)
      If (exis) Then
         Write (*, '(a)') 'Error: user termination of program in routin&
        &e: ' // trim (str)
         Call getunit (un)
         Open (un, File='TERMINATE', Action='write')
         Close (un, Status='delete')
         Call terminate
      End If
End Subroutine terminateqry
