! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine avscq(iqr, n, nmax, scrn, scieff)
  use modinput, only: input
  use mod_qpoint, only: iqmap
  use modxs, only: ivqr
  use invert

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, n, nmax
  complex(8), intent(in) :: scrn(n, n)
  complex(8), intent(out) :: scieff(nmax, nmax)

  ! Local variables
  integer :: iqrnr, j1, j2, flg
  real(8) :: clwt

  ! Find reduced q-point in non-reduced set
  iqrnr = iqmap(ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))

  ! Invert dielectric matrix
  call zinvert_hermitian(input%xs%bse%scrherm, scrn, scieff(:n, :n))

  do j1 = 1, n
    do j2 = 1, j1
      if((input%xs%bse%sciavqhd .and. (j1 .eq. 1) .and. (j2 .eq. 1))&
        & .or. (input%xs%bse%sciavqwg .and. (j1 .ne. 1)&
        & .and. (j2 .eq. 1)) .or. (input%xs%bse%sciavqwg .and. (j1 .eq. 1)&
        & .and. (j2 .ne. 1)) .or. (input%xs%bse%sciavqbd .and. (j1 .ne. 1)&
        & .and. (j2 .ne. 1))) then
        ! Numerical averaging on grids with extrapolation to continuum
        flg = 2
      else
        ! Analytic expression, no averaging
        flg = 0
      end if
      ! Generate the (averaged) symmetrized coulomb potential
      call genwiqggp(flg, iqrnr, j1, j2, clwt)
      ! Multiply with averaged coulomb potential
      scieff(j1, j2) = scieff(j1, j2) * clwt
      ! Set upper triangle
      scieff(j2, j1) = conjg(scieff(j1, j2))
    end do
  end do

end subroutine avscq
