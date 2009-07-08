

! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.

! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: zfinp2
! !INTERFACE:
complex(8) function zfinp2(ngp1, ngp2, igpig, zfmt1, zfmt2, zfir1, zfir2)
! !USES:
use modinput
  use modmain
! !INPUT/OUTPUT PARAMETERS:
!   zfmt1 : first complex function in spherical harmonics for all muffin-tins
!           (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfmt2 : second complex function in spherical harmonics for all muffin-tins
!           (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngrtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. In the intersitial region,
!   the integrand is multiplied with the smooth characteristic function,
!   $\tilde{\Theta}({\bf r})$, to remove the contribution from the muffin-tin.
!   See routines {\tt zfmtinp} and {\tt gencfun}. Based upon the routine
!   {\tt zfinp}.
!
! !REVISION HISTORY:
!   Created January 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: ngp1, ngp2, igpig(ngkmax)
  complex(8), intent(in) :: zfmt1(lmmaxvr, nrcmtmax, natmtot)
  complex(8), intent(in) :: zfmt2(lmmaxvr, nrcmtmax, natmtot)
  complex(8), intent(in) :: zfir1(ngp1)
  complex(8), intent(in) :: zfir2(ngp2)
  ! local variables
  integer::is, ia, ias, ig, igp1, igp2, iv(3)!!$, ir
  complex(8) zsum
  ! external functions
  complex(8) zfmtinp
  external zfmtinp
!!$! interstitial contribution
!!$zsum=0.d0
!!$do ir=1,ngrtot
!!$  zsum=zsum+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
!!$end do
!!$  zsum=zsum*omega/dble(ngrtot)
  do igp1=1, ngp1
     do igp2=1, ngp2
	iv(:) = ivg(:, igpig(igp1)) - ivg(:, igpig(igp2))
	ig = ivgig(iv(1), iv(2), iv(3))
	zsum=zsum+cfunig(ig)*conjg(zfir1(igp1))*zfir2(igp2)
     end do
  end do
  ! muffin-tin contribution
  do is=1, nspecies
     do ia=1, natoms(is)
	ias=idxas(ia, is)
	zsum = zsum + zfmtinp(input%groundstate%lmaxvr, nrcmt(is), rcmt(1, is), lmmaxvr, zfmt1(1, 1, ias), &
	     zfmt2(1, 1, ias))
     end do
  end do
  zfinp2=zsum
end function zfinp2
!EOC
