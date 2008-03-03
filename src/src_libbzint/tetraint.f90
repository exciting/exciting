!BOP
!
! !ROUTINE: tetraint 
!
! !INTERFACE:
       subroutine tetraint(nik,nt,nb,ebd,tetc,wtet,v,efer,opm,intopm)
!     
! !DESCRIPTION:
!
!   This subroutine integrates the values of an operator, which is essentially
! the same as tetiw, with the only difference that the value of this operator 
! is taken as input and the result is as output. We don't have to bother with 
! the weight things. It is all inside of the library.
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
       
       real(8), intent(in)    :: opm(nik,*) ! the values of the operator
!                                             that has to be integrated
       
! !OUTPUT PARAMETERS:
       
       real(8), intent(out)   :: intopm    ! the value of the integral
 
!  
! !LOCAL VARIABLES:
 
       integer(4) :: ik,ib
       
       real(8), allocatable :: kw(:,:)    ! Weight of each k-point
       
! !SYSTEM SUBROUTINES:
       
       intrinsic size
       
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
      
      allocate(kw(nband,nirkp))
!
!     use intw to calculate the weights
!
      call intw(efer,kw)

      intopm = 0

!      
!     calculate the integrated value of the operator with the weights got above
!  essential the same thing we do outside of the library when we call the tetiw
!  instead
!
      do ik=1,nik
        do ib=1,nband
          intopm=intopm+kw(ib,ik)*opm(ik,ib)
        enddo
      enddo
      
      end subroutine tetraint
      
!EOC
