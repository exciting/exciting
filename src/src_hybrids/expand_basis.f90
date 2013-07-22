
subroutine expand_basis(ik,iq)

    use modmain
    use modgw
    use modmpi,          only: rank
    use mod_hartreefock, only: vcmbsiz, vcmb

    implicit none
    integer(4), intent(in) :: ik, iq

    integer(4) :: jk, i, isym
    real(8)    :: tstart, tend

    call cpu_time(tstart)

!------------------------------------------------------------
! Coulomb potential is precalculated in init_mixedbasis.f90
!------------------------------------------------------------
! global array used for MB transformations
    mbsiz = vcmbsiz(iq)
    allocate(barcvm(matsiz,mbsiz))
    barcvm(:,:) = vcmb(1:matsiz,1:mbsiz,iq)
    
    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(eveck(nmatmax,nstfv,nspnfv))
    allocate(eveckp(nmatmax,nstfv,nspnfv))

    jk=kqid(ik,iq) ! index of k-q vector

    call getevecfv(vklnr(:,jk),vgklnr(:,:,:,jk),eveck)
    eveckp=conjg(eveck)
    call getevecfv(vklnr(:,ik),vgklnr(:,:,:,ik),eveck)

    call expand_evec(ik,'t')
    call expand_evec(jk,'c')

! Calculate the matrix elements $M^i_{nm}(\vec{k},\vec{q})$:
    call calcminm(ik,iq,1)

! Calculate the matrix elements M^i_{nm} where n is a core state
    if(iopcore == 0)then
        call calcminc(ik,iq,1)
    endif
    
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(barcvm)

    call cpu_time(tend)
    if (rank==0) call write_cputime(fgw,tend-tstart, 'EXPAND_BASIS')

    return
    
end subroutine
