
subroutine task_second_variation()
    use modmain
    use modgw, only: evalqp, evalks, evalfv, ibgw, nbgw, kset, Gkset, &
                     eferks, eferqp
    use mod_kpointset

    implicit none

    integer(4) :: ik, ib
    real(8) :: egap
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:)

    ! readjust global variables to fit the number of computed fv-QP states
    nstfv = nbgw
    nstsv = nstfv * nspinor

    !------------------------------------
    ! Read first-variational hamiltonian
    !------------------------------------
    if (allocated(evalqp)) deallocate(evalqp)
    allocate(evalqp(nstfv,kset%nkpt))
    if (allocated(evalks)) deallocate(evalks)
    allocate(evalks(nstfv,kset%nkpt))
    call readevalqp('EVALQP.OUT')
    if (allocated(evalfv)) deallocate(evalfv)
    allocate(evalfv(nstfv,kset%nkpt))
    evalfv(:,:) = evalqp(:,:)
    ! write(*,*) eferks, eferqp
    ! do ik = 1, kset%nkpt
    !     write(*,*) 'ik=', ik
    !     do ib = ibgw, nbgw
    !         write(*,*) ib, evalks(ib,ik), evalqp(ib,ik)
    !     end do
    !     write(*,*)
    ! end do
    deallocate(evalqp)
    deallocate(evalks)

    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,kset%nkpt))
    if (allocated(occsv)) deallocate(occsv)
    allocate(occsv(nstsv,kset%nkpt))
    filext = "_GW.OUT"
    do ik = 1, kset%nkpt
        call getevalsv(kset%vkl(:,ik), evalsv(:,ik))
        call getoccsv(kset%vkl(:,ik), occsv(:,ik))
    end do
    call readfermi()
    evalsv(:,:) = evalsv(:,:)-efermi
    efermi = 0.d0

    filext = "_KS.OUT"
    call writeeval()
    filext = "_GW.OUT"
    allocate(evalks(nstsv,kset%nkpt))
    evalks(:,:) = evalsv(:,:)
    eferks = 0.d0

    !------------------------------------
    ! Solve second-variational problem
    !------------------------------------
    if (allocated(evecfv)) deallocate(evecfv)
    allocate(evecfv(nmatmax,nstfv))
    if (allocated(evecsv)) deallocate(evecsv)
    allocate(evecsv(nstsv,nstsv))
    if (allocated(apwalm)) deallocate(apwalm)
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    do ik = 1, kset%nkpt
        call match(Gkset%ngk(1,ik), &
                   Gkset%gkc(:,1,ik), &
                   Gkset%tpgkc(:,:,1,ik), &
                   Gkset%sfacgk(:,:,1,ik),&
                   apwalm)
        call getevecfv(kset%vkl(:,ik), Gkset%vgkl(:,:,:,ik), evecfv)
        call seceqnsv(ik, apwalm, evalfv(:,ik), evecfv, evecsv)
    end do
    deallocate(evecfv)
    deallocate(apwalm)

    call fermi_exciting(input%groundstate%tevecsv, &
                        chgval, &
                        nstsv, kset%nkpt, evalsv, &
                        kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
                        efermi, egap, fermidos)
    call tetiw(kset%nkpt, kset%ntet, nstsv, evalsv, kset%tnodes, &
               kset%wtet, kset%tvol, efermi, occsv)
    do ik = 1, kset%nkpt
      do ib = 1, nstsv
        occsv(ib,ik) = dble(occmax)/kset%wkpt(ik)*occsv(ib,ik)
      end do
    end do

    ! write out the second-variational eigenvalues and occupation numbers
    filext = "_QP.OUT"
    call writeeval()
    filext = "_GW.OUT"

    nbgw = nstsv
    allocate(evalqp(nstsv,kset%nkpt))
    evalqp(:,:) = evalsv(:,:)
    eferqp = efermi

    ! store the eigenvalues to a file for bandstructure plot
    call putevalqp('EVALQPSV.OUT')

    call bandstructure_analysis('KS', 1, nstsv, kset%nkpt, evalks, eferks)
    call bandstructure_analysis('QP', 1, nstsv, kset%nkpt, evalqp, eferqp)

    deallocate(evalks)
    deallocate(evalqp)

end subroutine