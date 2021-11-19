
!
!BOP
! !ROUTINE: calc_vnlmat
! !INTERFACE:
!
subroutine calc_vnlmat
! !USES:
    use modmain
    use modgw
    use modfvsystem
    use mod_hybrids
    use modmpi
!
! !DESCRIPTION:
!   Calculates the APW matrix elements of the non-local potential
!EOP
!BOC
    implicit none
    type(evsystem) :: system
    integer :: ik
    integer :: ie1, ie2, nst
    integer :: nmatp
    complex(8), allocatable :: evec(:,:)
    complex(8), allocatable :: temp(:,:), temp1(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:,:)
    complex(8), external :: zdotc

    integer :: ikfirst, iklast
#ifdef MPI
    ikfirst = firstofset(rank,kset%nkpt)
    iklast  = lastofset(rank,kset%nkpt)
#else
    ikfirst = 1
    iklast  = kset%nkpt
#endif

    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vnlmat)) deallocate(vnlmat)
    allocate(vnlmat(nmatmax,nmatmax,ikfirst:iklast))
    vnlmat = zzero
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
    apwalm = zzero
    allocate(evec(nmatmax,nstfv))
    evec = zzero

    do ik = ikfirst, iklast
        if (input.groundstate.hybrid.method.eq."MB") then
            ! matching coefficients
            call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), &
            &          sfacgk(:,:,1,ik), apwalm(:,:,:,:,1))

            ! Hamiltonian and overlap setup
            nmatp = nmat(1,ik)
            call newsystem(system,input%groundstate%solver%packedmatrixstorage,nmatp)
            call MTRedirect(mt_hscf%main,mt_hscf%spinless)
            call hamiltonsetup(system, ngk(1, ik), apwalm, igkig(:, 1, ik), vgkc(:,:,1,ik))
            call overlapsetup(system, ngk(1, ik), apwalm, igkig(:, 1, ik), vgkc(:,:,1,ik))
            !write(*,*) 'overlap=', ik, sum(system%overlap%za)

            ! c
            call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evec)


            nst=min(nmatp, nstfv)
            ! conjg(c)*S
            allocate(temp(nst,nmatp))


           call zgemm('c', 'n', nst, nmatp, nmatp, &
           &          zone, evec(1:nmatp,:), nmatp, &
           &          system%overlap%za, nmatp, &
           &          zzero, temp, nst)

           ! Vnl*conjg(c)*S
           allocate(temp1(nst,nmatp))
           call zgemm('n', 'n', nst, nmatp, nst, &
           &          zone, vxnl(:,:,ik), nstfv, &
           &          temp, nst, zzero, &
           &          temp1, nst)

           ! V^{NL}_{GG'} = conjg[conjg(c)*S]*Vx*conjg(c)*S
           call zgemm('c', 'n', nmatp, nmatp, nst, &
           &          zone, temp, nst, &
           &          temp1, nst, zzero, &
           &          vnlmat(1:nmatp,1:nmatp,ik), nmatp)

            call deletesystem(system)
            deallocate(temp)
            deallocate(temp1)
        else
            ! calculate vnlmat=P+*P
            nmatp = nmat(1,ik)
            call zgemm('c', 'n', nmatp, nmatp, nstsv, &
            &          dcmplx(-1d0,0d0), pace(:,1:nmatp,ik), nstsv, &
            &          pace(:,1:nmatp,ik), nstsv, zzero, &
            &          vnlmat(1:nmatp,1:nmatp,ik), nmatp) ! nmatmax)
        endif

    end do ! ik

    deallocate(apwalm)
    deallocate(evec)

    return
end subroutine
