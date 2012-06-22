!BOP
! !ROUTINE: genvxcapw
! !INTERFACE:
      subroutine genvxcapw(ik,apwalm,vxc)
! !USES:
      use modmain
      use modgw
! !INPUT PARAMETERS:
      implicit none
      integer(4), intent(in)  :: ik      ! k-point number
! APW matching coefficients      
      complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
! !OUTPUT PARAMETERS:
      complex(8), intent(out) :: vxc(*) ! XC-potential matrix
! !LOCAL VARIABLES:
      integer(4) :: ia
      integer(4) :: is
      complex(8) :: v(nmatmax)
! !REVISION HISTORY:
!
! Created August 2006 (RGA)
! Revisited: May 2011 by DIN 
!
!EOP
!BOC            
      vxc(1:npmat(1,ik))=0.d0

!     muffin-tin contributions
      do is=1,nspecies
        do ia=1,natoms(is)
          call vxcaa(.false.,is,ia,ngk(1,ik),apwalm,v,vxc)
          call vxcalo(.false.,is,ia,ngk(1,ik),apwalm,v,vxc)
          call vxclolo(.false.,is,ia,ngk(1,ik),v,vxc)
        end do
      end do

!     interstitial contributions
      call vxcistl(.false.,ngk(1,ik),igkig(1,1,ik),vgkc(1,1,1,ik),v,vxc)
      
      return
      end subroutine genvxcapw
!EOC      
