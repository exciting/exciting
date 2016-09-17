!BOP
!
!!ROUTINE: readselfx
!
!!INTERFACE:
!
subroutine readselfx
!
!!DESCRIPTION:
! 
! This subroutine reads the exchange self-energy from file
!
!!USES:
    use modgw,   only : kset, ibgw, nbgw, selfex
    use m_getunit
       
!!LOCAL VARIABLES:
    integer(4) :: ib, nb, nk
    integer(4) :: fid
      
!!REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
    call getunit(fid)

    open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, nk
    close(fid)

    if (nk.ne.kset%nkpt) then
      write(*,*)'ERROR(readselfx): Wrong number of k-points'
      write(*,*)'    nk=', nk, '    nkpt=', kset%nkpt
      stop
    end if

    if ((ib.ne.ibgw).or.(nb.ne.nbgw)) then
      write(*,*)'ERROR(readselfx): Different number of bands'
      write(*,*)'    ib=',   ib, '    nb=', nb
      write(*,*)'  ibgw=', ibgw, '  nbgw=', nbgw
      stop
    end if
    
    open(fid,file='SELFX.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, nk, selfex
    close(fid)
      
    return
end subroutine
!EOC      
