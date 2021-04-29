
subroutine task_second_variation()

    use modmain
    use modgw,       only: evalqp, evalks, evalfv, ibgw, nbgw, kset, kqset, Gkset, Gkqset, &
                           eferks, eferqp
    use mod_hybrids, only: hybridhf
    use modmpi,      only: rank
    use modxs,       only: isreadstate0

    implicit none
    integer(4) :: ikp, ik, ib, nst
    real(8) :: egap
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:)
    real(8), allocatable :: occqp(:,:)

    ! write(*,*) 'nstfv, nstsv=', nstfv, nstsv
    ! write(*,*) 'ibgw, nbgw=', ibgw, nbgw

    ! Readjust the range of the FV states
    nstfv = nbgw
    nstsv = nstfv*nspinor

    !---------------------
    ! Read GW FV results
    !---------------------
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstfv,kset%nkpt))
    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(nstfv,kset%nkpt))
    call readevalqp('EVALQP.OUT', kset, 1, nstfv, evalks, eferks, evalqp, eferqp)
    deallocate(evalks)
    ! call bandstructure_analysis('QP FV', 1, nstfv, kset%nkpt, evalqp, eferqp)

    !---------------------
    ! Read KS SV energies
    !---------------------
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
    call bandstructure_analysis('KS+SO band structure', 1, nstsv, kset%nkpt, evalks, eferks)

    !------------------------------------
    ! Solve second-variational problem
    !------------------------------------
    if (.not. input%groundstate%spin%realspace) then
        ! Effective Hamiltonian Setup: Radial and Angular integrals (new)
        call MTInitAll(mt_hscf)
        call hmlint(mt_hscf)
    end if

    if (allocated(evalsv)) deallocate(evalsv)
    if (hybridhf) then
        allocate(evalsv(nstsv,kset%nkpt))
    else
        allocate(evalsv(nstsv,kqset%nkpt))
    end if
    evalsv(:,:) = 0.d0

    if (allocated(evecfv)) deallocate(evecfv)
    allocate(evecfv(nmatmax,nstfv))
    if (allocated(evecsv)) deallocate(evecsv)
    allocate(evecsv(nstsv,nstsv))
    if (allocated(apwalm)) deallocate(apwalm)
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    do ikp = 1, kset%nkpt
        if (hybridhf) then
            ! reduced k-grid is assumed by default
            call match(Gkset%ngk(1,ikp), &
                       Gkset%gkc(:,1,ikp), &
                       Gkset%tpgkc(:,:,1,ikp), &
                       Gkset%sfacgk(:,:,1,ikp),&
                       apwalm)
            call get_evec_gw(kset%vkl(:,ikp), Gkset%vgkl(:,:,:,ikp), evecfv)
            call seceqnsv(ikp, apwalm, evalqp(1:nstfv,ikp), evecfv, evecsv)
        else
            ! default GW definition corresponds to reducek=.false.
            ik = kset%ikp2ik(ikp)
            call match(Gkqset%ngk(1,ik), &
                       Gkqset%gkc(:,1,ik), &
                       Gkqset%tpgkc(:,:,1,ik), &
                       Gkqset%sfacgk(:,:,1,ik),&
                       apwalm)
            call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
            call seceqnsv(ik, apwalm, evalqp(1:nstfv,ikp), evecfv, evecsv)
        end if
    end do
    deallocate(evecfv)
    deallocate(evecsv)
    deallocate(apwalm)
    call mt_hscf%release()

    ! QP energies
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstsv,kset%nkpt))
    do ikp = 1, kset%nkpt
        if (hybridhf) then
            evalqp(:,ikp) = evalsv(:,ikp)
        else
            ik = kset%ikp2ik(ikp)
            evalqp(:,ikp) = evalsv(:,ik)
        end if
    end do
    deallocate(evalsv)

    nst = min(int(chgval)+20, nstsv)
    call fermi_exciting(.true., &
                        chgval, &
                        nst, kset%nkpt, evalqp(1:nst,:), &
                        kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
                        eferqp, egap, fermidos)
    call bandstructure_analysis('G0W0+SO band structure', 1, nst, kset%nkpt, evalqp(1:nst,:), eferqp)
    call putevalqp('EVALQPSV.OUT', kset, 1, nstsv, evalks, eferks, evalqp, eferqp)

    ! Calculate state occupation numbers
    if (allocated(occqp)) deallocate(occqp)
    allocate(occqp(nstsv,kset%nkpt))
    call tetiw(kset%nkpt, kset%ntet, nstsv, evalqp, kset%tnodes, &
               kset%wtet, kset%tvol, eferqp, occqp)
    do ikp = 1, kset%nkpt
      do ib = 1, nstsv
        occqp(ib,ikp) = 1.d0/kset%wkpt(ikp)*occqp(ib,ikp)
      end do
    end do

    open(82, file='EVALQPSV.DAT', status='unknown')
    do ikp = 1, kset%nkpt
        write(82,'(a,i4,3f12.4)') '# ikp, vkl=', ikp, kset%vkl(:,ikp)
        do ib = 1, nstsv
            write(82,'(i6,3f18.6)') ib, evalks(ib,ikp)-eferks, evalqp(ib,ikp), occqp(ib,ikp)
        end do
        write(82,*); write(82,*)
    end do
    close(82)

    deallocate(evalks)
    deallocate(evalqp)
    deallocate(occqp)

end subroutine
