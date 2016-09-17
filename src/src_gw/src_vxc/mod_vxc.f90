
module mod_vxc

! APW-APW exchange-correlation  integrals
      real(8), allocatable :: vxcraa(:,:,:,:,:,:)
      
! local-orbital-APW exchange-correlation  integrals
      real(8), allocatable :: vxcrloa(:,:,:,:,:)
      
! local-orbital-local-orbital exchange-correlation  integrals
      real(8), allocatable :: vxcrlolo(:,:,:,:)

! G-space interstitial exchange-correlation potential
      complex(8), allocatable :: vxcig(:)
      
! diagonal matrix elements of the exchange-correlation potential
      complex(8), allocatable :: vxcnn(:,:)
      
end module
