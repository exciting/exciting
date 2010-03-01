
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine deldynmat
		use modmain
        implicit none
        integer :: ip,is,ia,iq
        logical :: exist
        Do ip = 1, 3
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               Do iq = 1, nqpt
                  Call phfext (iq, is, ia, ip, filext)
                  Inquire (File='DYN'//trim(filext), Exist=Exist)
                  If ( exist) Then
                     Open (50, File='DYN'//trim(filext))
                     close(50, status='delete')
                  End If
               End Do
            End Do
         End Do
      End Do
end subroutine
