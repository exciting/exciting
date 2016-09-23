
! Copyright (C) 2010 S. Sagmeister and  C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpchgs
! !INTERFACE:
subroutine genpchgs(ik,evecfv,evecsv)
! !USES:
      use modinput
      use modmain
! !DESCRIPTION:
!  Generate partial charges for each state $j$, atom $\alpha$ and for each
!  $(lm)$ combination.
!
! !REVISION HISTORY:
!   Created 2010 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: ik
  complex(8), intent(in) :: evecfv(nmatmax, nstfv, nspnfv)
  complex(8), intent(in) :: evecsv(nstsv, nstsv)
  ! local variables
  integer :: lmax,lmmax,ispn,is,ia,ias,ist,l,m,lm
  real(8) :: t1,t2
  Complex (8), Allocatable :: dmat (:, :, :, :, :)
  Complex (8), Allocatable :: apwalm (:, :, :, :, :)
  ! set lmax for maximum accuracy
  lmax=input%groundstate%lmaxvr
  lmmax = (lmax+1) ** 2
  Allocate (dmat(lmmax, lmmax, nspinor, nspinor, nstsv))
  Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  ! find the matching coefficients
  Do ispn = 1, nspnfv
    Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :,ispn, ik), &
     & sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
  End Do
  ! average band character over spin for all atoms
  Do is = 1, nspecies
    Do ia = 1, natoms (is)
      ias = idxas (ia, is)
      ! generate the diagonal of the density matrix
      Call gendmat (.True., .True., 0, lmax, is, ia, ngk(:,ik), apwalm, &
       & evecfv, evecsv, lmmax, dmat)
      Do ist = 1, nstsv
        t1 = wkpt(ik) * occsv(ist,ik)
        If (abs(t1) .gt. input%groundstate%epsocc) Then
          Do l = 0, lmax
            Do m = - l, l
              t2 = 0.d0
              lm = idxlm (l, m)
              Do ispn = 1, nspinor
                t2 = t2 + dble (dmat(lm, lm, ispn, ispn, ist))
              End Do
              ! add for current k-point
              chgpart (lm, ias, ist) = chgpart (lm, ias, ist) + dble(t2)*t1
            End Do
          End Do
        end if
      ! end loop over states
      End Do
    End Do
  End Do
  Deallocate (dmat, apwalm)
end subroutine
!EOC
