!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phdelete
      Use modmain
      Implicit None
! delete the eigenvector files
      Call delevec
! delete the eigenvalue files
      Open (70, File=trim(scrpath)//'EVALFV'//trim(filext))
      Close (70, Status='DELETE')
      Open (70, File=trim(scrpath)//'EVALSV'//trim(filext))
      Close (70, Status='DELETE')
! delete the occupancy file
      Open (70, File=trim(scrpath)//'OCCSV'//trim(filext))
      Close (70, Status='DELETE')
! delete the STATE.OUT file
      Open (50, File='STATE'//trim(filext))
      Close (50, Status='DELETE')
      Return
End Subroutine
