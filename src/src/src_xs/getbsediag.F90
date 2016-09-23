!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine getbsediag
      Use modxs
      Use m_getunit
      Implicit None
  ! local variables
      Integer :: un
      Real (8) :: re, im
      Call getunit (un)
      Open (un, File='BSEDIAG.OUT', Action='read', Form='formatted', &
     & Status='old')
      Read (un,*) re, im
      bsed = cmplx (re, im, 8)
      Read (un,*) re, im
      bsedl = cmplx (re, im, 8)
      Read (un,*) re, im
      bsedu = cmplx (re, im, 8)
      Read (un,*) re, im
      bsedd = cmplx (re, im, 8)
      Close (un)
End Subroutine getbsediag
