


! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqn


subroutine seceqn(ik, evalfv, evecfv, evecsv)
  ! !USES:
use modinput
  use modmain
  use modmpi
  use sclcontroll
  use diisinterfaces

  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  !   evecsv : second-variational eigenvectors (out,complex(nstsv,nstsv))
  ! !DESCRIPTION:
  !   Solves the first- and second-variational secular equations. See routines
  !   {\tt match}, {\tt seceqnfv}, {\tt seceqnss} and {\tt seceqnsv}.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !EOP
  !BOC
  implicit none
  ! arguments
  integer, intent(in) :: ik
  real(8), intent(out) :: evalfv(nstfv, nspnfv)
  complex(8), intent(out) :: evecfv(nmatmax, nstfv, nspnfv)
  complex(8), intent(out) :: evecsv(nstsv, nstsv)
  ! local variables
  integer::ispn


  ! allocatable arrays
  complex(8), allocatable :: apwalm(:, :, :, :, :)
  allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)
#ifdef KSMP
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO
#endif
  !
  !-IMPORTANT: the first-variational spinor index and the k-point index have been
  ! swapped in the following arrays: ngk, igkig, vgkl, vgkc, gkc, tpgkc, sfacgk
  !
  do ispn=1, nspnfv
     ! find the matching coefficients
     call match(ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, ik), &
	  sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
     ! solve the first-variational secular equation
     if(doDIIScycle()) then
	call DIISseceqnfv(ik, ispn, apwalm(:, :, :, :, ispn), vgkc(:, :, ispn, ik), evalfv, evecfv)

	if (ik.eq.lastk(rank)) diiscounter=diiscounter+1

     else if(doARPACKiteration()) then
	call  iterativearpacksecequn(ik, ispn, apwalm(1, 1, 1, 1, ispn), &
	     vgkc(1, 1, ispn, ik), evalfv, evecfv)
     else if (doLAPACKsolver()) then
	if (tseqit) then
           ! iteratively
	   call seceqnit(nmat(ispn, ik), ngk(ispn, ik), igkig(:, ispn, ik), vkl(:, ik), &
		vgkl(:, :, ispn, ik), vgkc(:, :, ispn, ik), apwalm(:, :, :, :, ispn), evalfv(:, ispn), &
		evecfv(:, :, ispn))
	else
           ! directly
	   call seceqnfv(nmat(ispn, ik), ngk(ispn, ik), igkig(:, ispn, ik), &
		vgkc(:, :, ispn, ik), apwalm(:, :, :, :, ispn), evalfv(:, ispn), evecfv(:, :, ispn))
	end if
     else if(.true.) then
	write(*, *)"error in solverselect secequn.F90"

     endif
  end do
#ifdef KSMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
  if (isspinspiral()) then
     ! solve the spin-spiral second-variational secular equation
     call seceqnss(ik, apwalm, evalfv, evecfv, evecsv)
  else
     ! solve the second-variational secular equation
     call seceqnsv(ik, apwalm, evalfv, evecfv, evecsv)
  end if

  deallocate(apwalm)
  return
end subroutine seceqn
!EOC
