!BOP
!
! !ROUTINE: calcdos
!
! !INTERFACE:
      subroutine calcdos(nik,nb,eb,ntet,tetc,wtet,vt,nsp,emax,emin, &
     &                   np,dos)
!
! !DESCRIPTION:
!  This subroutine calculated the density of states in the given energy range
! !USES:
      implicit none     
! !INPUT PARAMETERS:

      integer, intent(in) :: nik                 ! Number of irreducible k-points
      integer, intent(in) :: nb                  ! Maximum number of bands
      integer, intent(in) :: nsp                 ! =1 / 2 for spin unpolarized/polarized caculations 
      real(8), intent(in) :: eb(nb,nik,nsp)      ! Band energies
      integer, intent(in) :: ntet                ! Number of tetrahedra
      integer, intent(in) :: tetc(4,*)           ! id. numbers of the corners of the tetrahedr
      integer, intent(in) :: wtet(*)             ! weight of each tetrahedron
      real(8), intent(in) :: vt                  ! the volume of the tetrahedra
      real(8), intent(in) :: emin, emax          ! energy range
      integer             :: np                  ! number of points
     
      
! !OUTPUT PARAMETERS:      
      real(8), intent(out) :: dos(np,nsp)         ! the density of states

! !REVISION HISTORY:
!
! Created 10th. March 2004 by RGA
! Modified by Jiang
!
! !LOCAL VARIABLES:

      integer :: i,isp
      real(8) :: ein,sfact
      real(8), external :: dostet

!EOP
!BOC
      if(nsp.eq.1) then 
        sfact=2.d0
      else 
        sfact=1.d0
      endif 

      do isp=1,nsp
        do i=1,np
          ein=emin+real(i-1)*(emax-emin)/(np-1)
          dos(i,isp)=sfact*dostet(nik,nb,eb(:,:,isp),ntet,tetc,wtet,vt,ein)
        enddo 
      enddo 

      end subroutine calcdos
!EOC        
          

      
      
            



      
