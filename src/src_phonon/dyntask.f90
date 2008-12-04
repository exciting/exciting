
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dyntask(fnum,iq,is,ia,ip)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
integer, intent(out) :: iq
integer, intent(out) :: is
integer, intent(out) :: ia
integer, intent(out) :: ip
! local variables
logical exist
ip=1
is=1
ia=1
iq=1
do ip=1,3
  do is=1,nspecies
    do ia=1,natoms(is)
      do iq=1,nqpt
        call phfext(iq,is,ia,ip,filext)
        inquire(file='DYN'//trim(filext),exist=exist)
        if (.not.exist) then
          open(fnum,file='DYN'//trim(filext),action='WRITE',form='FORMATTED')
          return
        end if
      end do
    end do
  end do
end do
write(*,*)
write(*,'("Info(dyntask): Nothing more to do")')
write(*,*)
stop
return
end subroutine
