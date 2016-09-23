!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genveffig
! !INTERFACE:
!
!
Subroutine genveffig
! !USES:
      Use modmain
! !DESCRIPTION:
!   Generates the Fourier transform of the effective potential in the
!   intersitial region. The potential is first multiplied by the characteristic
!   function which zeros it in the muffin-tins. See routine {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ig, ifg
! allocatable arrays
      Complex (8), Allocatable :: zfft (:)
      Allocate (zfft(ngrtot))
! multiply effective potential with characteristic function
      zfft (:) = veffir (:) * cfunir (:)
! Fourier transform to G-space
      Call zfftifc (3, ngrid,-1, zfft)
      Do ig = 1, ngvec
         ifg = igfft (ig)
         veffig (ig) = zfft (ifg)
      End Do
      Deallocate (zfft)
      Return
End Subroutine
!EOC
