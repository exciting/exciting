!BOP
!
! !ROUTINE: tetraconv 
!
! !INTERFACE:
       subroutine tetraconv(nik,nt,nb,ebd,tetc,linkt,v,efer,opm,intopm)

! !USES:
 
       use tetra_internal
       
       implicit none      

! !INPUT PARAMETERS:
 
       integer(4), intent(in) :: nik        ! Number of irreducible 
!                                             k-points
       
       integer(4), intent(in) :: nt         ! Number of tetrahedra
       
       integer(4), intent(in) :: nb         ! Number of bands
       
       real(8), target, intent(in) :: ebd(nik,*)   ! Band energies
       
       integer(4), target, intent(in) :: tetc(4,*) ! id. numbers of the
!                                                    corners of the 
!                                                    tetrahedra
   
       integer(4), target, intent(in) :: linkt(*)  ! weight of each 
!                                                    tetrahedron
       
       real(8), intent(in)    :: v        ! the volume of the tetrahedra
 
       real(8), intent(in)    :: efer     ! fermi energy
       
       real(8), intent(in)    :: opm(nik,nik,nb,nb) ! the values of the 
!                                                 operator to be 
!                                                 integrated
       
! !OUTPUT PARAMETERS:
       
       real(8), intent(out)   :: intopm    ! the value of the integral
       
! !DESCRIPTION: 
! 
!   This subroutine integrates the values of an operator which diagonal
!  values are given in \verb"opm" using the improved tetrahedron method
!  
! !LOCAL VARIABLES:
 
       integer(4) :: ik,jk,ib,jb
       
       real(8), allocatable :: kw(:,:,:,:)   ! Weight of each k-point
       
! !SYSTEM SUBROUTINES:
       
       intrinsic size
       
       external convw

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
      tetln => linkt(1:ntet)
      eband   => ebd(1:nik,1:nband)
      
      allocate(kw(nirkp,nirkp,nband,nband))
      
      call convw(efer,kw)
      intopm = 0
!      write(6,*)'tetraconv'
!      do ik=1,nik
!        write(6,*)'ik =',ik
!        do jk=1,nik
!          write(6,*)'jk =',jk
!          do ib=1,nband
!            write(6,*)(kw(ik,jk,ib,jb),jb=1,nband)
!            do jb=1,nband
!              intopm=intopm+kw(ik,jk,ib,jb)*opm(ik,jk,ib,jb)
!            enddo
!          enddo
!        enddo
!      enddo
      
      end subroutine tetraconv
      
!EOC
