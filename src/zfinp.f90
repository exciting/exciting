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
Complex (8) Function zfinp (tsh, zfmt1, zfmt2, zfir1, zfir2)
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
!   Modified March 2014 (UW)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tsh
      Complex (8), Intent (In) :: zfmt1 (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zfmt2 (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (In) :: zfir1 (ngrtot)
      Complex (8), Intent (In) :: zfir2 (ngrtot)
! local variables
      Integer :: is, ia, ias, ir
      Complex (8) zsum,zsum0
! external functions
      Complex (8) zfmtinp
      External zfmtinp
      zsum = 0.d0
! interstitial contribution

#ifdef USEOMP
!$OMP PARALLEL PRIVATE (ir,zsum0)
 zsum0= 0.d0
!$OMP DO  
      Do ir = 1, ngrtot
         zsum0 = zsum0 + cfunir (ir) * conjg (zfir1(ir)) * zfir2 (ir)
      End Do
!$OMP END DO
!$OMP ATOMIC
     zsum = zsum +zsum0
!$OMP END PARALLEL  
#else
      Do ir = 1, ngrtot
         zsum = zsum + cfunir (ir) * conjg (zfir1(ir)) * zfir2 (ir)
      End Do
#endif
      zsum = zsum * omega / dble (ngrtot)
! muffin-tin contribution
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zsum = zsum + zfmtinp (tsh, input%groundstate%lmaxvr, &
           & nrcmt(is), rcmt(:, is), lmmaxvr, zfmt1(:, :, ias), &
           & zfmt2(:, :, ias))
         End Do
      End Do
      zfinp = zsum
      Return
End Function
!EOC
