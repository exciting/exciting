
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfinp
! !INTERFACE:
complex(8) function zfinp(tsh,zfmt1,zfmt2,zfir1,zfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the muffin-tin functions are in spherical harmonics
!           (in,logical)
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngrtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. In the intersitial region,
!   the integrand is multiplied with the characteristic function, to remove the
!   contribution from the muffin-tin. See routines {\tt zfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
complex(8), intent(in) :: zfmt1(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) :: zfmt2(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) :: zfir1(ngrtot)
complex(8), intent(in) :: zfir2(ngrtot)
! local variables
integer is,ia,ias,ir
complex(8) zsum
! external functions
complex(8) zfmtinp
external zfmtinp
zsum=0.d0
! interstitial contribution
do ir=1,ngrtot
  zsum=zsum+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
end do
zsum=zsum*omega/dble(ngrtot)
! muffin-tin contribution
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zsum=zsum+zfmtinp(tsh,lmaxvr,nrcmt(is),rcmt(:,is),lmmaxvr,zfmt1(:,:,ias), &
     zfmt2(:,:,ias))
  end do
end do
zfinp=zsum
return
end function
!EOC

