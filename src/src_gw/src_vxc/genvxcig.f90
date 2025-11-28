module genvxcig
  implicit none
  private

  public :: generate_vxcig
contains

!> Generates the Fourier transform of the xc potential in the
!> intersitial region. The potential is first multiplied by the characteristic
!> function which zeros it in the muffin-tins. See routine {\tt gencfun}.
subroutine generate_vxcig(vxcir, vxcig)  
    use mod_Gvector, only: igfft, ngrid, cfunir
    use m_zfftifc, only: zfftifc
    use precision, only: dp, i32
    
    real(dp), intent(in) :: vxcir(:)
    complex(dp), intent(out) :: vxcig(:)

    integer(i32) :: ig, ifg, ng, ngrtot
    complex(dp), allocatable :: zfft(:)

    ngrtot = size( vxcir )
    allocate( zfft(ngrtot) )
    
    ! multiply effective potential with smooth characteristic function
    zfft(:) = vxcir(:)*cfunir(:)

    ! Fourier transform to G-space
    call zfftifc( 3, ngrid, -1, zfft )

    ng = size( vxcig )
    do ig = 1, ng
      ifg = igfft(ig)
      vxcig(ig) = zfft(ifg)
    end do
    
end subroutine
end module