!BOP
!
! !ROUTINE: writeselfx
!
! !INTERFACE:
      subroutine writeselfx

! !DESCRIPTION:
! 
! This subroutine writes the self-energy to file
!
! !USES:
!
      use modmain
      use modgw     
       
! !LOCAL VARIABLES:
      
      integer(4) :: ie   !(Counter) Runs over bands
      integer(4) :: ikp  !(Counter) Runs over k-points

! !REVISION HISTORY:
!
! Created 16.08.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC
      open(92,file='SELFX.OUT',action='WRITE',form='FORMATTED')
      write(92,*)'# band       selfex [eV]'
      do ikp = 1, nkpt
        write(92,1) vkl(1:3,ikp), ikp, ibgw, nbgw
        do ie = ibgw, nbgw
          write(92,2) ie, real(selfex(ie,ikp))*hev, aimag(selfex(ie,ikp))*hev
        enddo
      enddo    
      close(92)
    1 format(3e19.12,3i6)
    2 format(i4,4x,2f19.12)     
      
      return
      end subroutine writeselfx
!EOC      
