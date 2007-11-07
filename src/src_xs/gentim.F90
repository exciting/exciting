
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentim(sec_in,hrs,days,hours,minutes,seconds)
  implicit none
  ! arguments
  real(8), intent(in) :: sec_in
  real(8), intent(out) :: hrs
  integer, intent(out) :: days,hours,minutes,seconds
  ! local variables
  integer :: isecs

  hrs=sec_in/3600.d0
  isecs=nint(sec_in)
  days=isecs/(60*60*24)
  hours=mod(isecs/(60*60),24)
  minutes=mod(isecs/60,60)
  seconds=mod(isecs,60)
end subroutine gentim

character(256) function stringtim(sec,hrs,d,h,m,s)
  use a2str
  implicit none
  ! arguments
  real(8), intent(in) :: sec,hrs
  integer, intent(in) :: d,h,m,s

  stringtim=&
       trim(r2str(sec,'(f12.2)'))//' sec; '//&
       trim(r2str(hrs,'(f12.2)'))//' hrs; ( '//&
       trim(i2str(d,'(i4)'))//' d, '//&
       trim(i2str(h,'(i3.2)'))//' h, '//&
       trim(i2str(m,'(i3.2)'))//' m, '//&
       trim(i2str(s,'(i3.2)'))//' s )'
end function stringtim
