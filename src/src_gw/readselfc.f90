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
      integer(4) :: ib, nb, nk, no
      
! !REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
!
      open(93,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
      read(93) ib, nb, nk, no, selfec
      close(93)

      if (nk.ne.nkpt) then
        write(6,*)'ERROR(readselfc): Wrong number of k-points'
        write(6,*)'    nk=', nk, '    nkpt=', nkpt
        stop
      end if

      if ((ib.ne.ibgw).or.(nb.ne.nbgw)) then
        write(6,*)'WARNING(readselfc): Different band bounds'
        write(6,*)'    ib=',   ib, '    nb=', nb
        write(6,*)'  ibgw=', ibgw, '  nbgw=', nbgw
      end if

      if (no.ne.nomeg) then
        write(6,*)'ERROR(readselfc): Wrong number of frequencies'
        write(6,*)'    no=', no, '    nomeg=', nomeg
        stop
      end if
     
      return
      end subroutine
!EOC          
            
