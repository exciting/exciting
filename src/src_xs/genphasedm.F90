! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genphasedm
! !INTERFACE:
Subroutine genphasedm(iq, jsym, nmax, n, phfdm, tphf)
! !USES:
  use modinput
  use modxs, only: igqig
  use mod_Gvector, only: ivg
  use mod_symmetry, only: vtlsymc
  use mod_constants, only: twopi
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   integer(4), iq : Index of non-reduced q-point
!   integer(4), jsym : Index of crystal symmetry operation that maps qr onto q
!   integer(4), nmax : Maximal number of G+q vectors over all q
!   integer(4), n    : Number of G+q vectors for current q
! OUT:
!   complex(8), phfdm(nmax, nmax) : Phase factor for the dielectric matrix 
!                                   due to translations
!   logical, tphf : True if non-trivial phase appears at least for one (G,Gp) component
!
! !DESCRIPTION:
!   The routine computed the phase factors needed to reduce the dielectric matrix
!   at any q-point to corresponding matrix at a reduces q-point. The phases occur 
!   due to the translational vector $\vec{t}$ of the crystal symmetry operation 
!   chosen to reduce $\vec{q}$ to $\vec{q}_\text{r}$. The rotaional effects of the
!   symmetry operation are accounted for differently by remapping 
!   the set of $\left\{\vec{G}+\vec{q}\right\}$ vectors to the set of 
!   $\left\{\vec{G}+\vec{q}_\text{r}\right\}$ vectors. It is
!   $\epsilon_{\vec{G},\vec{G}'}(\vec{q}) = 
!     e^{i 2*pi (\vec{G}-\vec{G}') \cdot \vec{t}}
!     \epsilon_{\tilde{\vec{G}},\tilde{\vec{G}}'}(\vec{q}_r)$
!   
! !REVISION HISTORY:
!   Added to documentation scheme. Note: Description may be buggy. (Aurich 2016)
!EOP
!BOC

  implicit none

  ! Arguments
  integer(4), intent(in) :: iq, jsym, nmax, n
  complex(8), intent(out) :: phfdm(nmax, nmax)
  logical, intent(out) :: tphf

  ! Local variables
  real(8), parameter :: epsortho = 1.d-12
  real(8) :: vtl(3), t1, t2, t3
  integer(4) :: igq1, igq2, ivg1(3), ivg2(3), iv(3)

  ! Loop over G+q
  do igq1 = 1, n

    ! Get 3d integer coordinates of G corresponding to G+q
    ivg1(:) = ivg(:, igqig(igq1, iq))

    ! Loop over G'+q
    do igq2 = igq1, n

      ! G-vector difference G2 = G-G'
      ivg2(:) = ivg1(:) - ivg(:, igqig(igq2, iq))
      
      ! Translation vector t of the crystal symmetry operation vtl(s)
      vtl = vtlsymc(:, jsym)

      ! Map back to [0,1)
      call r3frac(input%structure%epslat, vtl, iv)

      ! 2pi * G2*t
      t1 = twopi * dot_product(dble(ivg2), vtl)
      t2 = cos(t1)
      t3 = sin(t1)

      if(abs(t2) .lt. epsortho) t2 = 0.d0
      if(abs(t3) .lt. epsortho) t3 = 0.d0

      ! Phase factor for dielectric matrix (due to translations)
      ! e^{i 2*pi (G-G')*t}
      phfdm(igq1, igq2) = cmplx(t2, t3, 8)
      phfdm(igq2, igq1) = conjg(phfdm(igq1, igq2))

      if(input%xs%dbglev .gt. 2) then
        write(40, '(a, i5, 2x, 2i5, 2x, i5, 2g18.10)') 'q, g, gp, jsym, phf',&
          & iq, igq1, igq2, jsym, phfdm(igq1, igq2)
      end if

    ! End loop over (G,Gp)-vectors
    end do

  end do

  ! Occurrence of non-trivial phase for q-point
  tphf = .false.
  if(any(abs(phfdm(1:n, 1:n)-1.d0) .gt. input%structure%epslat)) tphf = .true.

end subroutine genphasedm
!EOC
