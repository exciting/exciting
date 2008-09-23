
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqn
subroutine seceqn(ik,evalfv,evecfv,evecsv)
! !USES:
use modmain
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
real(8), intent(out) :: evalfv(nstfv,nspnfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! loop over first-variational spins (nspnfv=2 for spin-spirals only)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do ispn=1,nspnfv
! find the matching coefficients
  call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
   sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
! solve the first-variational secular equation
  if (tseqit) then
! iteratively
    call seceqnit(nmat(ispn,ik),ngk(ispn,ik),igkig(:,ispn,ik),vkl(:,ik), &
     vgkl(:,:,ispn,ik),vgkc(:,:,ispn,ik),apwalm(:,:,:,:,ispn),evalfv(:,ispn), &
     evecfv(:,:,ispn))
  else
! directly
    call seceqnfv(nmat(ispn,ik),ngk(ispn,ik),igkig(:,ispn,ik), &
     vgkc(:,:,ispn,ik),apwalm(:,:,:,:,ispn),evalfv(:,ispn),evecfv(:,:,ispn))
  end if
end do
!$OMP END DO
!$OMP END PARALLEL
if (spinsprl) then
! solve the spin-spiral second-variational secular equation
  call seceqnss(ik,apwalm,evalfv,evecfv,evecsv)
else
! solve the second-variational secular equation
  call seceqnsv(ik,apwalm,evalfv,evecfv,evecsv)
end if
deallocate(apwalm)
return
end subroutine
!EOC

