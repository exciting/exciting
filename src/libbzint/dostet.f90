
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: dostet
!
! !INTERFACE:
      real(8) function dostet(nik,nbd,eband,ntet,tetc,wtet,vt,e)
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

      integer(4), intent(in) :: nik ! Number of irreducible k-points
      
      integer(4), intent(in) :: nbd ! Maximum number of bands
      
      real(8), intent(in)    :: eband(nbd,nik) ! Band energies
      
      integer(4), intent(in) :: ntet        ! Number of tetrahedra
      
      integer(4), intent(in) :: tetc(4,*)! id. numbers of the corners
!                                             of the tetrahedra
  
      integer(4), intent(in) :: wtet(*)  ! weight of each tetrahedron
      
      real(8), intent(in)    :: vt       ! the volume of the tetrahedra

      real(8), intent(in)    :: e        !  energy


! !LOCAL VARIABLES:

      integer(4) :: itet,i,ib
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
      do itet=1,ntet
        do ib=1,nbd 
          do i=1,4
            ee(i)=eband(ib,tetc(i,itet))
          enddo
          call sort(4,ee)
          dostet=dostet+wtet(itet)*dos1t(ee,e,vt)
        enddo
      enddo
      
      end function dostet
      
!EOC      
          
        
      
