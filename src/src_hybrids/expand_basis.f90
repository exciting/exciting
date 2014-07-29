
subroutine expand_basis(ik,iq)

    use modmain
    use modgw
    use modmpi, only: rank
    
    implicit none
    integer(4), intent(in) :: ik, iq

    integer(4) :: jk, i, isym
    real(8)    :: tstart, tend

    call cpu_time(tstart)

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(eveck(nmatmax,nstfv,nspnfv))
    allocate(eveckp(nmatmax,nstfv,nspnfv))

    jk = kqid(ik,iq) ! index of k-q vector

    call getevecfv(vklnr(:,jk),vgklnr(:,:,:,jk),eveck)
    eveckp = conjg(eveck)
    call getevecfv(vklnr(:,ik),vgklnr(:,:,:,ik),eveck)

    call expand_evec(ik,'t')
    call expand_evec(jk,'c')

! Calculate the matrix elements $M^i_{nm}(\vec{k},\vec{q})$:
    call calcminm(ik,iq)

! Calculate the matrix elements M^i_{nm} where n is a core state
    if (iopcore == 0) call calcminc(ik,iq)
    
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    call cpu_time(tend)
    if (rank==0) call write_cputime(fgw,tend-tstart, 'EXPAND_BASIS')

    return
end subroutine
