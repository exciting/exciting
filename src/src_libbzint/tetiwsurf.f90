!BOP
!
! !ROUTINE: tetiwsurf
!
! !INTERFACE:
       subroutine tetiwsurf(nik,nt,nb,ebd,tetc,wtet,v,omeg,iwsurf)
!     
! !DESCRIPTION:
!
!   This subroutine calculates the weight on one k-point of a certain band for 
! a normal surface integration which is not q-dependent. 

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
 
       real(8), intent(in)    :: omeg       ! fermi energy
       
       
! !OUTPUT PARAMETERS:
       
       real(8), intent(out)   :: iwsurf(nb,nik)    ! the value of the weight
       
       
! !INTRINSIC ROUTINES:
       
       intrinsic size
       
! !EXTERNAL ROUTINES:
       
       external intwsurf
       
! !REVISION HISTORY:
!
!   Created: 10th. Jan 2005 by XZL
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
!     intwsurf is analogous to intw, difference is that it is for surface integration
!       
      call intwsurf(omeg,iwsurf)

      end subroutine tetiwsurf
      
!EOC
