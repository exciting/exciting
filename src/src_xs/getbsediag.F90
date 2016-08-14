! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine getbsediag
  use modxs, only: bsed, bsedl, bsedu, bsedd
  use m_getunit

  implicit none

  ! Local variables
  integer :: un
  real(8) :: re, im

  call getunit(un)
  open(un, file='BSEDIAG.OUT', action='read',&
    & form='formatted', status='old')

  read(un,*) re, im
  bsed = cmplx(re, im, 8)
  read(un,*) re, im
  bsedl = cmplx(re, im, 8)
  read(un,*) re, im
  bsedu = cmplx(re, im, 8)
  read(un,*) re, im
  bsedd = cmplx(re, im, 8)

  close(un)
end subroutine getbsediag
