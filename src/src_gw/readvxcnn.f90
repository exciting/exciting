!BOP
!
! !ROUTINE: writevxcnn
!
! !INTERFACE:
      subroutine readvxcnn

! !DESCRIPTION:
! 
! This subroutine reads the matrix elements $v^{xc}_{nn}(\vec{k})$  to file
!
! !USES:
!
      use modmain
      use modgw
       
! !LOCAL VARIABLES:
      
      integer(4) :: ist, i  !(Counter) Runs over bands
      integer(4) :: ikp     !(Counter) Runs over k-points
      real(8)    :: kvec(3)
      real(8)    :: vxc

! !REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
      open(66,file='VXCNN.OUT',action='READ',form='FORMATTED')
      do ikp=1,nkpt
        read(66,*) kvec
        do ist=ibgw,nbgw
          read(66,*) i, vxc
          vxcnn(ist,ikp)=cmplx(vxc/hev,0.0d0,8)
        enddo
      enddo    
      close(66)    
      
      return
      end subroutine
!EOC      
