
subroutine task_second_variation()

    use modmain
    use modgw, only: evalqp, evalks, evalfv, ibgw, nbgw, kset, kqset, Gkqset, &
                     eferks, eferqp
    use modmpi, only: rank

    implicit none
    integer(4) :: ikp, ik, ib
    real(8) :: egap
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:)

    print*, 'nstfv, nstsv=', nstfv, nstsv
    print*, 'ibgw, nbgw=', ibgw, nbgw

    ! ! readjust global variables to fit the number of computed fv-QP states
    ! nstfv = nbgw ! <-- this should be consistent with nstfv from GS
    ! nstsv = nstfv*nspinor

    !---------------------
    ! Read GW FV results
    !---------------------
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstfv,kset%nkpt))
    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(nstfv,kset%nkpt))
    call readevalqp('EVALQP.OUT', kset, ibgw, nbgw, evalks, eferks, evalqp, eferqp)
    deallocate(evalks)

    !---------------------
    ! Read KS SV energies
    !---------------------
    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(nstsv,kset%nkpt))
    filext = "_GW.OUT"
    do ikp = 1, kset%nkpt
        ik = kset%ikp2ik(ikp)
        call getevalsv(kqset%vkl(:,ik), evalks(:,ikp))
    end do
    ! find Fermi energy
    call fermi_exciting(.true., &
                        chgval, &
                        nstsv, kset%nkpt, evalks, &
                        kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
                        eferks, egap, fermidos)
    call bandstructure_analysis('KS SV', 1, nstsv, kset%nkpt, evalks, eferks)

    !------------------------------------
    ! Solve second-variational problem
    !------------------------------------
    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,kset%nkpt))
    if (allocated(evecfv)) deallocate(evecfv)
    allocate(evecfv(nmatmax,nstfv))
    if (allocated(evecsv)) deallocate(evecsv)
    allocate(evecsv(nstsv,nstsv))
    if (allocated(apwalm)) deallocate(apwalm)
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    do ikp = 1, kset%nkpt
        ik = kset%ikp2ik(ikp)
        call match(Gkqset%ngk(1,ik), &
                   Gkqset%gkc(:,1,ik), &
                   Gkqset%tpgkc(:,:,1,ik), &
                   Gkqset%sfacgk(:,:,1,ik),&
                   apwalm)
        call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
        call seceqnsv(ikp, apwalm, evalqp(1:nstfv,ikp), evecfv, evecsv)
    end do
    deallocate(evecfv)
    deallocate(evecsv)
    deallocate(apwalm)

    ! QP energies
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstsv,kset%nkpt))
    evalqp(:,:) = evalsv(:,:)
    deallocate(evalsv)
    call fermi_exciting(.true., &
                        chgval, &
                        nstsv, kset%nkpt, evalqp, &
                        kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
                        eferqp, egap, fermidos)
    call bandstructure_analysis('QP SV', 1, nstsv, kset%nkpt, evalqp, eferqp)

    ! store the eigenvalues to a file for bandstructure plot
    call putevalqp('EVALQPSV.OUT', kset, 1, nstsv, evalks, eferks, evalqp, eferqp)

    deallocate(evalks)
    deallocate(evalqp)

end subroutine