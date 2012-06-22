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
      write(93,*)'# band          selfec [eV]        selfcs1 [eV]        selfcs2 [eV]'
      do ikp = 1, nkpt
        write(93,*) 'ikp =', ikp
        do ie = ibgw, nbgw
          write(93,*) 'band nr. =', ie
          do iom = 1, nomeg
            write(93,2) freqs(iom), selfec(ie,ikp,iom)*hev,                   &
     &                  selfcs1(ie,ikp,iom)*hev, selfcs2(ie,ikp,iom)*hev
          enddo ! iom
        enddo ! ie
      enddo !ikp
      close(93)
    2 format(g15.5,' ',2g15.5,' ',2g15.5,' ',2g15.5)
      
      return
      end subroutine
!EOC          
            
