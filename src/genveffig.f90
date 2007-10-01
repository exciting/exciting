
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genveffig
! !INTERFACE:
subroutine genveffig
! !USES:
use modmain
! !DESCRIPTION:
!   Generates the Fourier transform of the effective potential in the
!   intersitial region. The potential is first multiplied by the characteristic
!   function which zeros it in the muffin-tins. See routine {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ig,ifg
! allocatable arrays
complex(8), allocatable :: zfft(:)
allocate(zfft(ngrtot))
! multiply effective potential with characteristic function
zfft(:)=veffir(:)*cfunir(:)
! Fourier transform to G-space
call zfftifc(3,ngrid,-1,zfft)
do ig=1,ngvec
  ifg=igfft(ig)
  veffig(ig)=zfft(ifg)
end do
deallocate(zfft)
return
end subroutine
!EOC
