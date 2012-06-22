!BOP
!
! !ROUTINE: readselfc
!
! !INTERFACE: 
      subroutine readselfc

! !DESCRIPTION:
!
! This subroutine writes the selfenergy to file
!
! !USES:

      use modmain
      use modgw

! !LOCAL VARIABLES:
      
      implicit none
      
      integer(4) :: ie, ib    ! (Counter) Runs over bands.
      integer(4) :: ikp       ! (Counter) Runs over k-points.
      integer(4) :: iom       ! (Counter) Runs over frequencies.
      real(8)    :: om
      real(8)    :: scr, sci
      
! !REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
!
      open(93,file='SELFC.OUT',form='FORMATTED',status='UNKNOWN')
      read(93,*)
      do ikp = 1, nkpt
        read(93,*)
        do ie = ibgw, nbgw
          read(93,*) ib
          if ((ib<ibgw).or.(ib>nbgw)) then
            write(6,*)'ERROR(readselfc): Inconsistent input parameters'
            write(6,*)'ibgw=',ibgw,'nbgw=',nbgw
          end if
          do iom = 1, nomeg
            read(93,*) om, scr, sci
            selfec(ie,ikp,iom)=cmplx(scr,sci,8)/hev
          enddo ! iom
        enddo ! ie
      enddo !ikp
      close(93)
     
      return
      end subroutine
!EOC          
            
