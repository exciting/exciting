!BOP
!
! !ROUTINE: writeselfe
!
! !INTERFACE: 
      subroutine writeselfc

! !DESCRIPTION:
!
! This subroutine writes the selfenergy to disc
!
! !USES:

      use modmain
      use modgw

! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: ie        ! (Counter) Runs over bands.
      integer(4) :: ikp       ! (Counter) Runs over k-points.
      integer(4) :: iom       ! (Counter) Runs over frequencies.

! !REVISION HISTORY:
!
! Created 20.07.05 by RGA.
! Revisited June 2011 by DIN
!
!EOP
!BOC
!
      open(93,file='SELFC.OUT',form='FORMATTED',status='UNKNOWN')
      do ikp = 1, nkpt
        do ie = ibgw, nbgw
          write(93,*)'# k-vector        ikp        band'
          write(93,*) vkl(:,ikp), ikp, ie
          write(93,*)'# omega           selfec [eV]'
          do iom = 1, nomeg
            write(93,1) freqs(iom), real(selfec(ie,ikp,iom))*hev, &
           &                       aimag(selfec(ie,ikp,iom))*hev
          enddo ! iom
        enddo ! ie
      enddo !ikp
      close(93)
    1 format(f15.5,' ',2f19.12)
      
      return
      end subroutine
!EOC          
            
