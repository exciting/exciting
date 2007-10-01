
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dyntask(fnum,iq,is,ia)
use modmain
implicit none
! arguments
integer, intent(in) :: fnum
integer, intent(out) :: iq
integer, intent(out) :: is
integer, intent(out) :: ia
! local variables
logical exist
integer i,j,m(3),n(3)
! external functions
integer gcd
external gcd
iq=1
is=1
ia=1
do iq=1,nqpt
  do i=1,3
    if (ivq(i,iq).ne.0) then
      j=gcd(ivq(i,iq),ngridq(i))
      m(i)=ivq(i,iq)/j
      n(i)=ngridq(i)/j
    else
      m(i)=0
      n(i)=0
    end if
  end do
  do is=1,nspecies
    do ia=1,natoms(is)
      write(filext,'("_Q",2I2.2,"_",2I2.2,"_",2I2.2,"_S",I2.2,"_A",I3.3,&
       &".OUT")') m(1),n(1),m(2),n(2),m(3),n(3),is,ia
      inquire(file='DYN'//trim(filext),exist=exist)
      if (.not.exist) then
        open(fnum,file='DYN'//trim(filext),action='WRITE',form='FORMATTED')
        return
      end if
    end do
  end do
end do
write(*,*)
write(*,'("Info(dyntask): Nothing more to do")')
write(*,*)
stop
return
end subroutine
