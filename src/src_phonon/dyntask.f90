!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine dyntask (fnum, iq, is, ia, ip, status)
      Use modmain
      Implicit None
! arguments
      Integer, Intent (In) :: fnum
      Integer, Intent (Out) :: iq
      Integer, Intent (Out) :: is
      Integer, Intent (Out) :: ia
      Integer, Intent (Out) :: ip
      character(*), intent(out) :: status
! local variables
      Logical :: exist
      ip = 1
      is = 1
      ia = 1
      iq = 1
      status="unfinished"
      Do ip = 1, 3
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               Do iq = 1, nqpt
                  Call phfext (iq, is, ia, ip, filext)
                  Inquire (File='DYN'//trim(filext), Exist=Exist)
                  If ( .Not. exist) Then
                     Open (fnum, File='DYN'//trim(filext), Action='WRIT&
                    &E', Form='FORMATTED')
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
