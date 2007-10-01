
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getoccsv(vpl,occsvp)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: occsvp(nstsv)
! local variables
integer isym,ik
integer recl,nstsv_
real(8) vkl_(3)
! external functions
real(8) r3taxi
external r3taxi
! find the k-point number
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nstsv_,occsvp
!$OMP CRITICAL
open(70,file=trim(scrpath)//'OCCSV'//trim(filext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nstsv_,occsvp
close(70)
!$OMP END CRITICAL
if (r3taxi(vkl(1,ik),vkl_).gt.epslat) then
  write(*,*)
  write(*,'("Error(getoccsv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" OCCSV.OUT  : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getoccsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" OCCSV.OUT  : ",I8)') nstsv_
  write(*,*)
  stop
end if
return
end subroutine

