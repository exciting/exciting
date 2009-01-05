
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phdelete
use modmain
implicit none
! delete the eigenvector files
call delevec
! delete the eigenvalue files
open(70,file=trim(scrpath)//'EVALFV'//trim(filext))
close(70,status='DELETE')
open(70,file=trim(scrpath)//'EVALSV'//trim(filext))
close(70,status='DELETE')
! delete the occupancy file
open(70,file=trim(scrpath)//'OCCSV'//trim(filext))
close(70,status='DELETE')
! delete the STATE.OUT file
open(50,file='STATE'//trim(filext))
close(50,status='DELETE')
return
end subroutine

