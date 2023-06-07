!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
!
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: zfinp
! !INTERFACE:
Complex (8) Function zfinpmt (tsh, zfmt1, zfmt2)
! !USES:
      Use modmain
      Use modinput
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the muffin-tin functions are in spherical harmonics
!           (in,logical)
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))

! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions inside muffin-tin. 
!   The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. 
!   See routines {\tt zfmtinp} and {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!   Modified March 2014 (UW)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Complex (8), Intent (In) :: zfmt1 (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zfmt2 (lmmaxvr, nrcmtmax, natmtot)
! local variables
      Integer :: is, ia, ias
      Complex (8) zsum,zsum0
! external functions
      Complex (8) zfmtinp
      External zfmtinp
      zsum = 0.d0
! interstitial contribution

! muffin-tin contribution
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zsum = zsum + zfmtinp (tsh, input%groundstate%lmaxvr, &
           & nrcmt(is), rcmt(:, is), lmmaxvr, zfmt1(:, :, ias), &
           & zfmt2(:, :, ias))
         End Do
      End Do
      zfinpmt = zsum
      Return
End Function
!EOC
