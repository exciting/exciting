
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

    write(*,*) 'nstfv, nstsv=', nstfv, nstsv
    write(*,*) 'ibgw, nbgw=', ibgw, nbgw

    !---------------------
    ! Read GW FV results
    !---------------------
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstfv,kset%nkpt))
    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(nstfv,kset%nkpt))
    call readevalqp('EVALQP.OUT', kset, 1, nstfv, evalks, eferks, evalqp, eferqp)
    deallocate(evalks)
    call bandstructure_analysis('QP FV', 1, nstfv, kset%nkpt, evalqp, eferqp)
    ! do ikp = 1, kset%nkpt
    !     do ib = 1, nstfv
    !         write(81,*) ikp, ib, evalqp(ib,ikp)
    !     end do
    !     write(81,*); write(81,*)
    ! end do

    !---------------------
    ! Read KS SV energies
    !---------------------
    filext = "_GW.OUT"

    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(nstsv,kset%nkpt))
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

    if (allocated(evalsv)) deallocate(evalsv) ! init_gw is done for the full k-grid => vkl, etc. corresponds to vklnr
    allocate(evalsv(nstsv,kqset%nkpt))
    evalsv(:,:) = 0.d0

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
        call seceqnsv(ik, apwalm, evalqp(1:nstfv,ikp), evecfv, evecsv) ! a lot of problems with global GS variables!
        ! write(*,*) 'ik, ngk=', ik, ngk(1,ik), vkl(:,ik)
        ! do ib = 1, nstsv
        !     write(82,*) ib, evalsv(ib,ik)
        ! end do
        ! write(82,*); write(82,*)
    end do
    deallocate(evecfv)
    deallocate(evecsv)
    deallocate(apwalm)

    ! QP energies
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstsv,kset%nkpt))
    do ikp = 1, kset%nkpt
        ik = kset%ikp2ik(ikp)
        evalqp(:,ikp) = evalsv(:,ik)
    end do
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