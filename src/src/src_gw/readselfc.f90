!BOP
!
!!ROUTINE: readselfc
!
!!INTERFACE:
!
subroutine readselfc
!
!!DESCRIPTION:
!
! This subroutine reads the correlation self-energy from file
!
!!USES:
    use modgw, only : kset, freq, ibgw, nbgw, selfec
    use m_getunit

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ib, nb, nk, no
    integer(4) :: fid
      
!!REVISION HISTORY:
!
! Created June 2011 by DIN
!
!EOP
!BOC
    call getunit(fid)

    open(fid,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, no, nk
    close(fid)

    if (nk.ne.kset%nkpt) then
      write(*,*)'ERROR(readselfc): Wrong number of k-points'
      write(*,*)'    nk=', nk, '    nkpt=', kset%nkpt
      stop
    end if

    if ((ib.ne.ibgw).or.(nb.ne.nbgw)) then
      write(*,*)'WARNING(readselfc): Different number of bands'
      write(*,*)'    ib=',   ib, '    nb=', nb
      write(*,*)'  ibgw=', ibgw, '  nbgw=', nbgw
      stop
    end if

    if (no.ne.freq%nomeg) then
      write(*,*)'ERROR(readselfc): Wrong number of frequencies'
      write(*,*)'    no=', no, '    freq%nomeg=', freq%nomeg
      stop
    end if
    
    open(fid,file='SELFC.OUT',form='UNFORMATTED',status='UNKNOWN')
    read(fid) ib, nb, no, nk, selfec
    close(fid)
     
    return
end subroutine
!EOC          
            
