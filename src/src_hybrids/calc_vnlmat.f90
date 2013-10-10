
subroutine calc_vnlmat
    
    use modmain
    use modgw
    use modfvsystem
    use mod_hybrids
    use modmpi
            
    implicit none
    type(evsystem) :: system
    integer(4) :: ik, ispn
    integer(4) :: nmatp
    integer(4) :: ie1, ie2
    complex(8), allocatable :: vnl(:,:)
    complex(8), allocatable :: evec(:,:), overlap(:,:)
    complex(8), allocatable :: temp(:,:), temp1(:,:)
    complex(8), allocatable :: apwalm(:,:,:,:,:)
      
    complex(8), external :: zdotc
!_______________________________________________________________________________
!
    allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))

#ifdef MPISEC
    do ik = firstk(rank), lastk(rank)
#else
    do ik = 1, nkpt
#endif
        vnlmat(:,:,ik) = zzero

        do ispn = 1, nspnfv

            nmatp = nmat(ispn,ik)

            call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
           &  sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
            
! Hamiltonian and overlap set up
            call newsystem(system,input%groundstate%solver%packedmatrixstorage,nmatp)
            call hamiltonandoverlapsetup(system,ngk(ispn,ik),apwalm, &
           &  igkig(:,ispn,ik),vgkc(:,:,ispn,ik))

            allocate(vnl(nstfv,nstfv))
            vnl(:,:) = vxnl(:,:,ik)

! S
            allocate(overlap(nmatp,nmatp))
            overlap(:,:) = system%overlap%za(:,:)
            if ((debug).and.(rank==0)) then
                call linmsg(fgw,'-',' Overlap ')
                do ie1 = 1, nmatp, 20
                    write(fgw,*) (overlap(ie1,ie2), ie2=1,nmatp,20)
                end do
                call linmsg(fgw,'-','')
            end if

! c
            allocate(evec(nmatp,nstfv))
            evec(1:nmatp,1:nstfv) = evecfv0(1:nmatp,1:nstfv,1,ik)
            if ((debug).and.(rank==0)) then
                call linmsg(fgw,'-',' EvecFV ')
                do ie1 = 1, nmatp, 20
                    write(fgw,*) (evec(ie1,ie2), ie2=1,nstfv)
                end do
            end if

! conjg(c)*S
            allocate(temp(nstfv,nmatp))
            call zgemm('c','n',nstfv,nmatp,nmatp, &
           &  zone,evec,nmatp,overlap,nmatp,zzero,temp,nstfv)
            deallocate(evec,overlap)

! Vnl*conjg(c)*S
            allocate(temp1(nstfv,nmatp))
            call zgemm('n','n',nstfv,nmatp,nstfv, &
           &  zone,vnl,nstfv,temp,nstfv,zzero,temp1,nstfv)
            deallocate(vnl)

! V^{NL}_{GG'} = conjg[conjg(c)*S]*Vx*conjg(c)*S
            allocate(vnl(nmatp,nmatp))
            call zgemm('c','n',nmatp,nmatp,nstfv, &
           &  zone,temp,nstfv,temp1,nstfv,zzero,vnl,nmatp)
            deallocate(temp,temp1)

            vnlmat(1:nmatp,1:nmatp,ik) = vnl(1:nmatp,1:nmatp)
            deallocate(vnl)

        end do ! ispn
    end do ! ik

#ifdef MPI
    call mpi_allgatherv_ifc(nkpt,nmatmax*nmatmax,zbuf=vnlmat)
    call barrier
#endif
    
    if ((debug).and.(rank==0)) then
        call linmsg(fgw,'-',' Vx_NL_GG ')
        do ie1 = 1, nmatp, 20
            write(fgw,*) (vnlmat(ie1,ie2,ik), ie2=1,nmatp,20)
        end do
        call linmsg(fgw,'-','')
    end if

    return
end subroutine
