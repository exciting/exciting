
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
    integer :: ie1, ie2
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
    vnlmat(:,:,:) = zzero
    
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(evec(nmatmax,nstfv))

    do ik = ikfirst, iklast

        ! matching coefficients
        call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik), &
        &          sfacgk(:,:,1,ik),apwalm(:,:,:,:,1))
        !write(*,*) 'apwalm=', ik, sum(apwalm)
            
        ! Hamiltonian and overlap setup 
        nmatp = nmat(1,ik)
        call newsystem(system,input%groundstate%solver%packedmatrixstorage,nmatp)
        call hamiltonandoverlapsetup(system,ngk(1,ik),apwalm, &
        &                            igkig(:,1,ik),vgkc(:,:,1,ik))
        !write(*,*) 'overlap=', ik, sum(system%overlap%za)

        ! S
        if (input%gw%debug.and.(rank==0)) then
            call linmsg(fgw,'-',' Overlap ')
            do ie1 = 1, nmatp, 100
                write(fgw,*) (system%overlap%za(:,:), ie2=1,nmatp,100)
            end do
            call linmsg(fgw,'-','')
        end if
        
        ! c
        call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evec)
        !write(*,*) 'evec=', ik, sum(evec(1:nmatp,:))
        
        if (input%gw%debug.and.(rank==0)) then
            call linmsg(fgw,'-',' EvecFV ')
            do ie1 = 1, nmatp, 100
                write(fgw,*) (evec(ie1,ie2), ie2=1,nstfv,10)
            end do
        end if

        ! conjg(c)*S
        allocate(temp(nstfv,nmatp))
        call zgemm('c','n',nstfv,nmatp,nmatp, &
        &          zone,evec(1:nmatp,:),nmatp, &
        &          system%overlap%za,nmatp, &
        &          zzero,temp,nstfv)
        !write(*,*) 'temp=', ik, sum(temp)

        ! Vnl*conjg(c)*S    
        allocate(temp1(nstfv,nmatp))
        call zgemm('n','n',nstfv,nmatp,nstfv, &
        &          zone,vxnl(:,:,ik),nstfv, &
        &          temp,nstfv,zzero, &
        &          temp1,nstfv)
        !write(*,*) 'temp1=', ik, sum(temp1)

        ! V^{NL}_{GG'} = conjg[conjg(c)*S]*Vx*conjg(c)*S
        call zgemm('c','n',nmatp,nmatp,nstfv, &
        &          zone,temp,nstfv, &
        &          temp1,nstfv,zzero, &
        &          vnlmat(1:nmatp,1:nmatp,ik),nmatp)
        !write(*,*) 'vnlmat=', ik, sum(vnlmat(:,:,ik))
            
        if (input%gw%debug.and.(rank==0)) then
            call linmsg(fgw,'-',' Vx_NL_GG ')
            do ie1 = 1, nmatp, 100
                write(fgw,*) (vnlmat(ie1,ie2,ik), ie2=1,nmatp,100)
            end do
            call linmsg(fgw,'-','')
        end if
        
        call deletesystem(system)
        deallocate(temp)
        deallocate(temp1)
            
    end do ! ik
    
    deallocate(apwalm)
    deallocate(evec)
    
    return
end subroutine
