

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine scrgeneigvec
  use modxs
  use m_genfilname
  implicit none
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  ! generate eigenvectors, eigenvalues, occupancies and APW MT coefficients
  call xsgeneigvec
  write(unitout, '("Info(scrgeneigvec): eigenvectors for screening finished")')
end subroutine scrgeneigvec
