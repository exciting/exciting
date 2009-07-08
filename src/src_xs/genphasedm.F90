

! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine genphasedm(iq, jsym, nmax, n, phfdm, tphf)
  use modmain
use modinput
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iq, jsym, nmax, n
  complex(8), intent(out) :: phfdm(nmax, nmax)
  ! true if non-trivial phase appears at least for one (G,Gp) component
  logical, intent(out) :: tphf
  ! local variables
  real(8), parameter :: epsortho=1.d-12
  real(8) :: vtl(3), t1, t2, t3
  integer :: igq1, igq2, ivg1(3), ivg2(3), iv(3)
  do igq1=1, n
     ivg1(:)=ivg(:, igqig(igq1, iq))
     do igq2=igq1, n
        ! G-vector difference
	ivg2(:)=ivg1(:)-ivg(:, igqig(igq2, iq))
        ! translation vector vtl(s)
	vtl=vtlsymc(:, jsym)
	call r3frac(input%structure%epslat, vtl, iv)
	t1=twopi*dot_product(dble(ivg2), vtl)
	t2=cos(t1)
	t3=sin(t1)
	if (abs(t2).lt.epsortho) t2=0.d0
	if (abs(t3).lt.epsortho) t3=0.d0
        ! phase factor for dielectric matrix (due to translations)
	phfdm(igq1, igq2)=cmplx(t2, t3, 8)
	phfdm(igq2, igq1)=conjg(phfdm(igq1, igq2))
	if (input%xs%dbglev.gt.2) then
	   write(40, '(a, i5, 2x, 2i5, 2x, i5, 2g18.10)') 'q, g, gp, jsym, phf', &
		iq, igq1, igq2, jsym, phfdm(igq1, igq2)
	end if
        ! end loop over (G,Gp)-vectors
     end do
  end do
  ! occurrance of non-trivial phase for q-point
  tphf=.false.
  if (any(abs(phfdm(1:n, 1:n)-1.d0).gt.input%structure%epslat)) tphf=.true.
end subroutine genphasedm
