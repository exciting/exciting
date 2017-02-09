!BOP
!
! !ROUTINE: idos
!
! !INTERFACE:
      real(8) function idos(nsp,nik,nbd,eband,ntet,tetc,wtet,vt,e)
!      
! !DESCRIPTION:
!
! This function calculates the integrated density of states at an energy
! e
!
     
! !USES:
      use order
      
      implicit none     
      
! !INPUT PARAMETERS:

      integer(4), intent(in):: nsp  !! number of spin channels 
      integer(4), intent(in) :: nik ! Number of irreducible k-points
      integer(4), intent(in) :: nbd ! Maximum number of bands
      real(8), intent(in) :: eband(nbd,nik,nsp) ! Band energies
      integer(4), intent(in) :: ntet        ! Number of tetrahedra
      integer(4), intent(in) :: tetc(4,*)! id. numbers of the corners of the tetrahedra
      integer(4), intent(in) :: wtet(*)  ! weight of each tetrahedron
      real(8), intent(in)    :: vt         ! the volume of the tetrahedra
      real(8), intent(in)    :: e        !  energy

! !LOCAL VARIABLES:

      integer(4) :: itet,i,ib,isp
      real(8) :: idt,sfact
      real(8), dimension(4) :: ee
      real(8), external :: intdos1t

! !REVISION HISTORY:
!
! Created:  3th. March 2004. by RGA 
!
!EOP
!
!BOC
      idos=0.0d0
      do isp=1,nsp
        do itet=1,ntet
          do ib=1,nbd 
            do i=1,4
              ee(i)=eband(ib,tetc(i,itet),isp)
            enddo
            call sort(4,ee)
            idt=intdos1t(ee,e,vt)
            idos=idos+wtet(itet)*idt
          enddo
        enddo
      enddo 
      idos=idos*2.d0/nsp
      return
      
      end function idos
      
!EOC      
          
        
      
