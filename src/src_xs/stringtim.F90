
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

character(256) function stringtim(sec,hrs,d,h,m,s)
  use modmain,only:filext
  implicit none
  ! arguments
  real(8), intent(in) :: sec,hrs
  integer, intent(in) :: d,h,m,s
  ! functions
  character(256), external :: r2str,i2str
  stringtim=&
       trim(r2str(sec,'(F12.2)'))//' sec; '//&
       trim(r2str(hrs,'(F12.2)'))//' hrs; ( '//&
       trim(i2str(d,'(I4)'))//' d, '//&
       trim(i2str(h,'(I3.2)'))//' h, '//&
       trim(i2str(m,'(I3.2)'))//' m, '//&
       trim(i2str(s,'(I3.2)'))//' s )'
end function stringtim
