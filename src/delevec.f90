
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine delevec
use modmain
implicit none
! delete the first-variational eigenvector file
open(70,file=trim(scrpath)//'EVECFV'//trim(filext))
close(70,status='DELETE')
! delete the second-variational eigenvector file
open(70,file=trim(scrpath)//'EVECSV'//trim(filext))
close(70,status='DELETE')
return
end subroutine

