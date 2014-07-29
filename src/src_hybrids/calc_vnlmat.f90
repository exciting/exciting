
subroutine calc_vnlmat
    
    use modmain
    use modgw
    use modfvsystem
    use mod_hybrids
    use modmpi
            
    implicit none
    type(evsystem) :: system
    integer(4) :: ik
    integer(4) :: ie1, ie2
    complex(8), allocatable :: vnl(:,:)
    complex(8), allocatable :: evec(:,:), overlap(:,:)
    complex(8), allocatable :: temp(:,:), temp1(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:,:)
    
    complex(8), allocatable :: evec_(:,:)
      
    complex(8), external :: zdotc

    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vnlmat)) deallocate(vnlmat)
    allocate(vnlmat(nmatmax,nmatmax,nkpt))

    call newsystem(system,input%groundstate%solver%packedmatrixstorage,nmatmax)
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
    allocate(evec(nmatmax,nstfv))
    allocate(temp(nstfv,nmatmax))
    allocate(temp1(nstfv,nmatmax))

#ifdef MPI
    do ik = firstk(rank), lastk(rank)
#else
    do ik = 1, nkpt
#endif
        call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik), &
        &          sfacgk(:,:,1,ik),apwalm(:,:,:,:,1))
            
! Hamiltonian and overlap set up
        system%overlap%za(:,:) = zzero
        call hamiltonandoverlapsetup(system,ngk(1,ik),apwalm, &
        &                            igkig(:,1,ik),vgkc(:,:,1,ik))

! S
        if ((debug).and.(rank==0)) then
            call linmsg(fgw,'-',' Overlap ')
            do ie1 = 1, nmatmax, 100
                write(fgw,*) (system%overlap%za(:,:), ie2=1,nmatmax,100)
            end do
            call linmsg(fgw,'-','')
        end if

! c
        evec(:,:) = zzero
        call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evec)
        if ((debug).and.(rank==0)) then
            call linmsg(fgw,'-',' EvecFV ')
            do ie1 = 1, nmatmax, 100
                write(fgw,*) (evec(ie1,ie2), ie2=1,nstfv,10)
            end do
        end if

! conjg(c)*S
        call zgemm('c','n',nstfv,nmatmax,nmatmax, &
        &          zone,evec,nmatmax, &
        &          system%overlap%za(:,:),nmatmax, &
        &          zzero,temp,nstfv)

! Vnl*conjg(c)*S
        call zgemm('n','n',nstfv,nmatmax,nstfv, &
        &          zone,vxnl(:,:,ik),nstfv, &
        &          temp,nstfv,zzero, &
        &          temp1,nstfv)

! V^{NL}_{GG'} = conjg[conjg(c)*S]*Vx*conjg(c)*S
        vnlmat(:,:,ik) = zzero
        call zgemm('c','n',nmatmax,nmatmax,nstfv, &
        &          zone,temp,nstfv, &
        &          temp1,nstfv,zzero, &
        &          vnlmat(:,:,ik),nmatmax)
            
        if ((debug).and.(rank==0)) then
            call linmsg(fgw,'-',' Vx_NL_GG ')
            do ie1 = 1, nmatmax, 100
                write(fgw,*) (vnlmat(ie1,ie2,ik), ie2=1,nmatmax,100)
            end do
            call linmsg(fgw,'-','')
        end if
            
    end do ! ik
    
    deallocate(apwalm)
    deallocate(evec)
    deallocate(temp)
    deallocate(temp1)
    call deleteystem(system)

#ifdef MPI
    call mpi_allgatherv_ifc(nkpt,nmatmax*nmatmax,zbuf=vnlmat)
    call barrier
#endif

    return
end subroutine
