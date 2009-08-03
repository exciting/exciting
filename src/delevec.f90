! Copyright (C) 2007 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine delevec
use modmain
implicit none
logical::exist
exist=.false.
! delete the first-variational eigenvector file
 inquire(file=trim(scrpath)//'EVECFV'//trim(filext), exist=exist)
 if (exist) then
open(71, file=trim(scrpath)//'EVECFV'//trim(filext))
close(71, status='DELETE')
endif
 inquire(file=trim(scrpath)//'EVECSV'//trim(filext), exist=exist)
 if (exist) then
! delete the second-variational eigenvector file
open(71, file=trim(scrpath)//'EVECSV'//trim(filext))
close(71, status='DELETE')
endif
return
end subroutine
