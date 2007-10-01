
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putoccsv(ik,occsvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: occsvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,occsvp
!$OMP CRITICAL
open(70,file=trim(scrpath)//'OCCSV'//trim(filext),action='WRITE', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
write(70,rec=ik) vkl(:,ik),nstsv,occsvp
close(70)
!$OMP END CRITICAL
return
end subroutine

