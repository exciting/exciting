
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genvxcig
! !INTERFACE:
subroutine genvxcig
! !USES:
    use modmain
    use modgw
! !DESCRIPTION:
!   Generates the Fourier transform of the xc potential in the
!   intersitial region. The potential is first multiplied by the characteristic
!   function which zeros it in the muffin-tins. See routine {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created August 2006 (RGA)
!   Modified September 2015 (Ronaldo) - to include DFT-1/2
!EOP
!BOC
    implicit none
! local variables
    integer ig,ifg
! allocatable arrays
    complex(8), allocatable :: zfft(:)
    allocate(zfft(ngrtot))
    if (allocated(vxcig)) deallocate(vxcig)
    allocate(vxcig(ngvec))
! multiply effective potential with smooth characteristic function
! Check if a DFT-1/2 calculation is required
    if (associated(input%groundstate%dfthalf)) then
      zfft(:)=(vxcir(:)+vhalfir(:))*cfunir(:)
    else
      zfft(:)=vxcir(:)*cfunir(:)
    end if
! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,zfft)
    do ig=1,ngvec
      ifg=igfft(ig)
      vxcig(ig)=zfft(ifg)
    end do
    deallocate(zfft)
    return
end subroutine
!EOC
