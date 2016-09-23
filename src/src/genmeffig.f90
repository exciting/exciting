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
Subroutine genmeffig
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
      Real(8) :: alpha,a2
      Parameter (alpha=1d0 / 137.03599911d0, a2=0.5d0*alpha**2)
! allocatable arrays
      Complex (8), Allocatable :: zfft (:)
      Allocate (zfft(ngrtot))
      if (input%groundstate%ValenceRelativity.ne."none") then
         if (allocated(meffig)) deallocate(meffig)
         allocate(meffig(ngvec))
         Do ig = 1, ngrtot
            zfft(ig)=cfunir(ig)/(1d0-veffir(ig)*a2)
         End Do
         Call zfftifc (3, ngrid,-1, zfft)
         Do ig = 1, ngvec
           ifg = igfft (ig)
           meffig (ig) = zfft (ifg)
         End Do
      endif
      if (input%groundstate%ValenceRelativity.eq."iora*") then
         if (allocated(m2effig)) deallocate(m2effig)
         allocate(m2effig(ngvec))
         Do ig = 1, ngrtot
            zfft(ig)=cfunir(ig)/((1d0-veffir(ig)*a2)**2)
         End Do
         Call zfftifc (3, ngrid,-1, zfft)
         Do ig = 1, ngvec
           ifg = igfft (ig)
           m2effig (ig) = zfft (ifg)
         End Do
      endif

      Deallocate (zfft)
      Return
End Subroutine
!EOC
