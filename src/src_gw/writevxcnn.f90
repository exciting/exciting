!BOP
!
! !ROUTINE: writevxcnn
!
! !INTERFACE:
      subroutine writevxcnn

! !DESCRIPTION:
! 
! This subroutine writes the matrix elements $v^{xc}_{nn}(\vec{k})$  to file
!
! !USES:
!
      use modmain
      use modgw
       
! !LOCAL VARIABLES:
      
      integer(4) :: ist  !(Counter) Runs over bands
      integer(4) :: ikp  !(Counter) Runs over k-points

! !REVISION HISTORY:
!
! Created 16.08.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC
      open(66,file='VXCNN.OUT',action='WRITE',form='FORMATTED')
      do ikp = 1, nkpt
        write(66,1) vkl(1:3,ikp),ikp,ibgw,nbgw
        do ist = ibgw, nbgw
          write(66,2) ist, real(vxcnn(ist,ikp))*hev
        enddo
      enddo    
      close(66)    
    1 format(3e19.12,3i6)
    2 format(i4,2f19.12)     
      
      return
      end subroutine writevxcnn
!EOC      
