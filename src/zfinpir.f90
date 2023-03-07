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
! !ROUTINE: zfinpir
! !INTERFACE:
Complex (8) Function zfinpir (zfir1, zfir2)
! !USES:
      Use modmain
      Use modinput
! !INPUT/OUTPUT PARAMETERS:
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngrtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions. 
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!   Modified March 2014 (UW)
!EOP
!BOC
      Implicit None
! arguments
      Complex (8), Intent (In) :: zfir1 (ngrtot)
      Complex (8), Intent (In) :: zfir2 (ngrtot)
! local variables
      Complex (8) zsum
      Integer :: ir
! interstitial contribution
      zsum = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DO PRIVATE (ir) REDUCTION(+:zsum)
      Do ir = 1, ngrtot
         zsum = zsum + conjg (zfir1(ir)) * zfir2 (ir)
      End Do
!$OMP END PARALLEL DO
#else
      Do ir = 1, ngrtot
         zsum = zsum + conjg (zfir1(ir)) * zfir2 (ir)
      End Do
#endif
      zfinpir = zsum * omega / dble (ngrtot)
      Return
End Function
!EOC
