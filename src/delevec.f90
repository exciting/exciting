! Copyright (C) 2007 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine delevec
      Use modmain
      Implicit None
      Logical :: exist
      exist = .False.
! delete the first-variational eigenvector file
      Inquire (File=trim(scrpath)//'EVECFV'//trim(filext), Exist=Exist)
      If (exist) Then
         Open (71, File=trim(scrpath)//'EVECFV'//trim(filext))
         Close (71, Status='DELETE')
      End If
      Inquire (File=trim(scrpath)//'EVECSV'//trim(filext), Exist=Exist)
      If (exist) Then
! delete the second-variational eigenvector file
         Open (71, File=trim(scrpath)//'EVECSV'//trim(filext))
         Close (71, Status='DELETE')
      End If
      Return
End Subroutine
