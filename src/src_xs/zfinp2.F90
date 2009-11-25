!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
!
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: zfinp2
! !INTERFACE:
Complex (8) Function zfinp2 (ngp1, ngp2, igpig, zfmt1, zfmt2, zfir1, &
& zfir2)
! !USES:
      Use modinput
      Use modmain
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
      Implicit None
  ! arguments
      Integer, Intent (In) :: ngp1, ngp2, igpig (ngkmax)
      Complex (8), Intent (In) :: zfmt1 (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zfmt2 (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zfir1 (ngp1)
      Complex (8), Intent (In) :: zfir2 (ngp2)
  ! local variables
      Integer is, ia, ias, ig, igp1, igp2, iv (3)
      Complex (8) zsum
  ! external functions
      Complex (8) zfmtinp
      External zfmtinp
  ! interstitial contribution
      Do igp1 = 1, ngp1
         Do igp2 = 1, ngp2
            iv (:) = ivg (:, igpig(igp1)) - ivg (:, igpig(igp2))
            ig = ivgig (iv(1), iv(2), iv(3))
            zsum = zsum + cfunig (ig) * conjg (zfir1(igp1)) * zfir2 &
           & (igp2)
         End Do
      End Do
  ! muffin-tin contribution
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zsum = zsum + zfmtinp (input%groundstate%lmaxvr, nrcmt(is), &
           & rcmt(1, is), lmmaxvr, zfmt1(1, 1, ias), zfmt2(1, 1, ias))
         End Do
      End Do
      zfinp2 = zsum
End Function zfinp2
!EOC
