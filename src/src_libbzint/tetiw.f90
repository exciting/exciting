!BOP
!
! !ROUTINE: tetiw 
!
! !INTERFACE:
       subroutine tetiw(nik,nt,nb,ebd,tetc,wtet,v,efer,iw)

!     
! !DESCRIPTION:
!
!   This subroutine gives the weight of on one k-point for a certain band
! in an operator integration.
!  

! !USES:
 
       use tetra_internal
       
       implicit none      
       
! !INPUT PARAMETERS:
 
       integer(4), intent(in) :: nik        ! Number of irreducible 
!                                             k-points
       
       integer(4), intent(in) :: nt         ! Number of tetrahedra
       
       integer(4), intent(in) :: nb         ! Number of bands
       
       real(8), target, intent(in) :: ebd(nb,nik)  ! Band energies
       
       integer(4), target, intent(in) :: tetc(4,*)  ! id. numbers of 
!                                                     the corners
!                                                     of the tetrahedra
   
       integer(4), target, intent(in) :: wtet(*)    ! weight of each 
!                                                     tetrahedron
       
       real(8), intent(in)    :: v         ! the volume of the tetrahedra
 
       real(8), intent(in)    :: efer       ! fermi energy
       
       
! !OUTPUT PARAMETERS:
       
       real(8), intent(out)   :: iw(nb,nik)    ! the value of the integral
       
       
! !INTRINSIC ROUTINES:
       
       intrinsic size
       
! !EXTERNAL ROUTINES:

       external intw
       
! !REVISION HISTORY:
!
!   Created: 4th. March 2004 by RGA
!
!EOP
!BOC
 
      nirkp = nik
      ntet  = nt
      nband = nb
      vt = v
      tetcorn => tetc(1:4,1:ntet)
      tetweig => wtet(1:ntet)
      eband   => ebd(1:nband,1:nik)
!
!     intw is called to calculate the weigths
!      
      call intw(efer,iw)
    
      end subroutine tetiw
      
!EOC
