!BOP
!
! !ROUTINE: readselfx
!
! !INTERFACE:
      subroutine readselfx

! !DESCRIPTION:
! 
! This subroutine writes the self-energy to file
!
! !USES:
!
      use modmain
      use modgw     
       
! !LOCAL VARIABLES:
      
      integer(4) :: ib, nb, nk
      
! !REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
      open(92,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
      read(92) ib, nb, nk, selfex
      close(92)

      if (nk.ne.nkpt) then
        write(6,*)'ERROR(readselfx): Wrong number of k-points'
        write(6,*)'    nk=', nk, '    nkpt=', nkpt
        stop
      end if

      if ((ib.ne.ibgw).or.(nb.ne.nbgw)) then
        write(6,*)'WARNING(readselfx): Different band bounds'
        write(6,*)'    ib=',   ib, '    nb=', nb
        write(6,*)'  ibgw=', ibgw, '  nbgw=', nbgw
      end if
      
      return
      end subroutine
!EOC      
