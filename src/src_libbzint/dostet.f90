!BOP
!
! !ROUTINE: dostet
!
! !INTERFACE:
      real(8) function dostet(nsp,nik,nbd,eband,ntet,tetc,wtet,vt,e)
!     
! !DESCRIPTION:
!
! This subroutine calculates the density of states at an energy e
!
!     
! !USES:
      use order
      
      implicit none     
      
! !INPUT PARAMETERS:

      integer, intent(in) :: nsp ! Number of spin
      integer, intent(in) :: nik ! Number of irreducible k-points
      integer, intent(in) :: nbd ! Maximum number of bands
      real(8), intent(in) :: eband(nbd,nik,nsp) ! Band energies
      integer, intent(in) :: ntet        ! Number of tetrahedra
      integer, intent(in) :: tetc(4,*)! id. numbers of the corners  of the tetrahedra
      integer, intent(in) :: wtet(*)  ! weight of each tetrahedron
      real(8), intent(in) :: vt         ! the volume of the tetrahedra
      real(8), intent(in) :: e      !  energy

! !LOCAL VARIABLES:

      integer :: itet,i,ib,isp
      real(8), dimension(4) :: ee
      real(8), external :: dos1t

! !REVISION HISTORY:
!
! Created:  3th. March 2004. by RGA 
!
!EOP
!
!BOC
      dostet=0.0d0
      do isp=1,nsp
        do itet=1,ntet
          do ib=1,nbd 
            do i=1,4
              ee(i)=eband(ib,tetc(i,itet),isp)
            enddo
            call sort(4,ee)
            dostet=dostet+wtet(itet)*dos1t(ee,e,vt)
          enddo
        enddo
      enddo 
      dostet=dostet*2.d0/nsp
      end function dostet
      
!EOC      
          
        
      
