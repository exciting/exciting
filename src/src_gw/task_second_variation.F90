
subroutine task_second_variation()
    use modmain
    use modgw, only: evalqp, evalks, evalfv, ibgw, nbgw, kset, Gkset, &
                     eferks, eferqp
    use modmpi, only: rank

    implicit none

    integer(4) :: ik, ib
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:)

    nstfv = nbgw
    nstsv = nstfv * nspinor
    ! write(*,*) 'nstfv, nstsv = ', nstfv, nstsv

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
    deallocate(evalqp)
    deallocate(evalks)
    ! write(*,*) eferks, eferqp
    ! do ik = 1, kset%nkpt
    !     write(*,*) 'ik=', ik
    !     do ib = ibgw, nbgw
    !         write(*,*) ib, evalks(ib,ik), evalqp(ib,ik)
    !     end do
    !     write(*,*)
    ! end do
    
    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,kset%nkpt))
    filext = "_GW.OUT"
    do ik = 1, kset%nkpt
        call getevalsv(kset%vkl(:,ik), evalsv(:,ik))
    end do
    call occupy
    evalsv(:,:) = evalsv(:,:)-efermi
   
    filext = "_KS.OUT"
    if (rank==0) call writeeval()
    filext = "_GW.OUT"

    allocate(evalks(nstsv,kset%nkpt))
    evalks(:,:) = evalsv(:,:)
    eferks = 0.d0

    !------------------------------------
    ! Solve second-variational problem
    !------------------------------------
    allocate(evecfv(nmatmax,nstfv))
    allocate(evecsv(nstsv,nstsv))
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
    call occupy
    deallocate(evecfv)
    ! write out the second-variational eigenvalues and occupation numbers
    filext = "_QP.OUT"
    if (rank==0) call writeeval()
    filext = "_GW.OUT"

    nbgw = nstsv
    allocate(evalqp(nstsv,kset%nkpt))
    evalqp(:,:) = evalsv(:,:)
    eferqp = efermi

    ! store the eigenvalues to a file for bandstructure plot
    call putevalqp('EVALQPSV.OUT')

    deallocate(evalks)
    deallocate(evalqp)

end subroutine