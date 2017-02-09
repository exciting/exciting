!BOP
!!ROUTINE: genvxcig
!!INTERFACE:
!
subroutine genvxcig
!
!!USES:
    use modmain, only : ngrid, ngrtot, igfft, &
    &                   vxcir, cfunir
    use modgw,   only : Gset
    use mod_vxc
    
!!DESCRIPTION:
!   Generates the Fourier transform of the xc potential in the
!   intersitial region. The potential is first multiplied by the characteristic
!   function which zeros it in the muffin-tins. See routine {\tt gencfun}.
!
!!REVISION HISTORY:
!   Created August 2006 (RGA) based on the routine genveffig.f90
!   Revisited October 2013 (DIN)
!EOP
!BOC
    implicit none
    ! local variables
    integer ig, ifg
    ! allocatable arrays
    complex(8), allocatable :: zfft(:)
    allocate(zfft(ngrtot))
    allocate(vxcig(Gset%ngvec))

    ! multiply effective potential with smooth characteristic function
    zfft(:) = vxcir(:)*cfunir(:)

    ! Fourier transform to G-space
    call zfftifc(3,ngrid,-1,zfft)

    do ig = 1, Gset%ngvec
      ifg = igfft(ig)
      vxcig(ig) = zfft(ifg)
    end do
    
    deallocate(zfft)
    return
end subroutine
!EOC
