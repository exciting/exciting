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
      Real(8) :: energyref,alpha,a2
      Parameter (alpha=1d0 / 137.03599911d0, a2=0.5d0*alpha**2)
! allocatable arrays
      Complex (8), Allocatable :: zfft (:)
      Allocate (zfft(ngrtot))
      if (input%groundstate%ValenceRelativity.eq."scalar") then
!         write(*,*) 'howdy'
         energyref=input%groundstate%energyref
!         if (.not.allocated(meffig))
         write(*,*) 'allocation' 
         allocate(meffig(ngvec))
!         zfft(:)=cfunir(:)/(1d0+(energyref-veffir(:))*a2)
         Do ig = 1, ngrtot
            zfft(ig)=cfunir(ig)/(1d0+(energyref-veffir(ig))*a2)
         End Do
         Call zfftifc (3, ngrid,-1, zfft)
         Do ig = 1, ngvec
           ifg = igfft (ig)
           meffig (ig) = zfft (ifg)
         End Do
!         write(*,*) 'done'
!         Do ig = 1, ngvec
!            zfft(:)
!         enddo
      endif

      Deallocate (zfft)
      Return
End Subroutine
!EOC
