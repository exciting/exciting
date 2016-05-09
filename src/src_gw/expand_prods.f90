!BOP
!
!!ROUTINE: expand_prods
!
!!INTERFACE:
!
subroutine expand_prods(ik,iq,iflag)
!
! !DESCRIPTION:
!
! This subroutine calculates $M^i_{cm}(k,q)$, $M^i_{nc}(k,q)$, 
! and $M^i_{nm}(k,q)$ matrix elements. 
! If flag>0, M^i_{nm} will be transformed to the v-diagonal basis set.
!
!!USES:
    use modmain
    use modgw
      
!!INPUT PARAMETERS:
    implicit none
    integer, intent(in) :: ik, iq
    integer(4), intent(in) :: iflag   ! < 0 -- Mcm
                                      ! = 0 -- Mnc + Mcm 
                                      ! > 0 -- Mnc

!!LOCAL VARIABLES:
    integer :: jk
    real(8) :: tstart, tend

!!REVISION HISTORY:
!
!EOP
!BOC
    call timesec(tstart)

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))

    jk = kqset%kqid(ik,iq) ! index of k-q vector

    ! read eigenvectors from file
    call getevecfvgw(jk,eveck)
    eveckp = conjg(eveck)
    call getevecfvgw(ik,eveck)

    ! calculate products A_{lm}*C_n
    call expand_evec(ik,'t')
    call expand_evec(jk,'c')
    
    !=================================================
    ! Loop over m-blocks in M^i_{nm}
    !=================================================

    ! full matrix
    allocate(minmmat(matsiz,ibgw:nbgw,1:nstsv))
    call calcminm2(ik,iq,ibgw,nbgw,1,nstfv,minmmat)
    
    !if (input%gw%coreflag=='all') then
    !  ! calculate M^i_{nc}
    !  allocate(minc(locmatsiz,nstart:nend,ncg))
    !  call calcminc(ik,iq,nstart,nend,minc)
    !endif

    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    
    call timesec(tend)
    time_eprod = time_eprod+tend-tstart

    return
end subroutine
!EOC

