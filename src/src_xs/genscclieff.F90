! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscclieff(iqr, nmax, n, scieff)

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, n, nmax
  complex(8), intent(out) :: scieff(nmax, nmax)

  ! Local variables
  logical :: tq0
  complex(8), allocatable :: scrn(:, :), scrnw(:, :, :), scrnh(:, :)

  ! External functions
  logical, external :: tqgamma

  allocate(scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3))

  ! Read screening from file
  call getscreen(iqr, n, scrnh, scrnw, scrn)
  tq0 = tqgamma(iqr)

  if(tq0) then
    ! Averaging using Lebedev-Laikov spherical grids
    call angavsc0(n, nmax, scrnh, scrnw, scrn, scieff)
  else
    ! Averaging using numerical method and extrapolation
    call avscq(iqr, n, nmax, scrn, scieff)
  end if

  deallocate(scrn, scrnw, scrnh)

end subroutine genscclieff
