!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dyntask (iq, is, ia, ip, status)
      Use modmain
      Implicit None
! arguments
      Integer, Intent (Out) :: iq
      Integer, Intent (Out) :: is
      Integer, Intent (Out) :: ia
      Integer, Intent (Out) :: ip
      character(*), intent(out) :: status
! local variables
      Logical :: exist
      character(256) :: chdummy
      ip = 1
      is = 1
      ia = 1
      iq = 1
      status="unfinished"
      Do ip = 1, 3
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               Do iq = 1, nqpt
                  Call phfext (iq, is, ia, ip, 0, 1, chdummy, filextdyn, chdummy)
                  Inquire (File='DYN'//trim(filextdyn), Exist=Exist)
                  If ( .Not. exist) Then
                     Return
                  End If
               End Do
            End Do
         End Do
      End Do
      Write (*,*)
      Write (*, '("Info(dyntask): Nothing more to do")')
      Write (*,*)
      status="finished"
End Subroutine
