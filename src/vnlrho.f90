
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: vnlrho
! !INTERFACE:
subroutine vnlrho(wfmt1,wfmt2,wfir1,wfir2,zrhomt,zrhoir)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   wfmt1  : muffin-tin part of wavefunction 1 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfmt2  : muffin-tin part of wavefunction 2 in spherical coordinates
!            (in,complex(lmmaxvr,nrcmtmax,natmtot,nspinor))
!   wfir1  : interstitial wavefunction 1 (in,complex(ngrtot))
!   wfir2  : interstitial wavefunction 2 (in,complex(ngrtot))
!   zrhomt : muffin-tin charge density in spherical harmonics
!            (out,complex(lmmaxvr,nrcmtmax,natmtot))
!   zrhoir : interstitial charge density (out,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the complex overlap charge density from two input wavefunctions:
!   $$ \rho({\bf r})\equiv\Psi_1^*({\bf r})\Psi_2({\bf r}). $$
!   Note that the muffin-tin wavefunctions are provided in spherical coordinates
!   and the returned density is in terms of spherical harmonic coefficients. See
!   also routine {\tt vnlrhomt}.
!
! !REVISION HISTORY:
!   Created November 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) ::  wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor)
complex(8), intent(in) ::  wfir1(ngrtot,nspinor)
complex(8), intent(in) ::  wfir2(ngrtot,nspinor)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(out) :: zrhoir(ngrtot)
! local variables
integer is,ia,ias,nr,ir
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
if (spinpol) allocate(zfmt(lmmaxvr,nrcmtmax))
! muffin-tin part
do is=1,nspecies
  nr=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    call vnlrhomt(is,wfmt1(1,1,ias,1),wfmt2(1,1,ias,1),zrhomt(1,1,ias))
    if (spinpol) then
! spin-polarised
      call vnlrhomt(is,wfmt1(1,1,ias,2),wfmt2(1,1,ias,2),zfmt)
      zrhomt(:,1:nr,ias)=zrhomt(:,1:nr,ias)+zfmt(:,1:nr)
    end if
  end do
end do
! interstitial part
if (spinpol) then
! spin-polarised
  do ir=1,ngrtot
    zrhoir(ir)=conjg(wfir1(ir,1))*wfir2(ir,1)+conjg(wfir1(ir,2))*wfir2(ir,2)
  end do
else
! spin-unpolarised
  do ir=1,ngrtot
    zrhoir(ir)=conjg(wfir1(ir,1))*wfir2(ir,1)
  end do
end if
if (spinpol) deallocate(zfmt)
return
end subroutine
!EOC

