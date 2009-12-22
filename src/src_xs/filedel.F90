!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_filedel
      Implicit None
Contains
!
!
      Subroutine filedel (fnam)
         Use m_getunit
         Implicit None
    ! arguments
         Character (*), Intent (In) :: fnam
    ! local variables
         Integer, Parameter :: verb = 0
         Integer :: un
         Logical :: existent, opened
    ! check if file exists
         Inquire (File=trim(fnam), Exist=existent)
         If ((verb .Gt. 0) .And. ( .Not. existent)) Then
            Write (*, '("Warning(filedel): attempted to delete non-exis&
           &tent file: ", a)') trim (fnam)
            Return
         End If
    ! check if file is opened
         Inquire (File=trim(fnam), Opened=Opened, Number=un)
    ! close file if opened
         If (opened) Then
            Close (un)
         End If
    ! open file for writing
         Call getunit (un)
         Open (un, File=trim(fnam), Action='write')
    ! delete file
         Close (un, Status='delete')
      End Subroutine filedel
!
End Module m_filedel
