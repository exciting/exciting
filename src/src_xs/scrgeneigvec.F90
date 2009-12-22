!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine scrgeneigvec
      Use modxs
      Use m_genfilname
      Implicit None
      Call genfilname (dotext='_SCR.OUT', setfilext=.True.)
  ! generate eigenvectors, eigenvalues, occupancies and APW MT coefficients
      Call xsgeneigvec
      Write (unitout, '("Info(scrgeneigvec): eigenvectors for screening&
     & finished")')
End Subroutine scrgeneigvec
