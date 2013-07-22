
subroutine add_vxnl(system,ik,nmatp)
    
    use modmain
    use modgw
    use modfvsystem
    use mod_hartreefock
    use modmpi
            
    implicit none
    type(evsystem), intent(InOut) :: system
    integer(4), intent(IN) :: ik
    integer(4), intent(IN) :: nmatp
   
    integer(4) :: ie1, ie2
    
    complex(8), allocatable :: vnl(:,:), vnlmat(:,:)
    complex(8), allocatable :: evec(:,:), overlap(:,:)
    complex(8), allocatable :: temp(:,:), temp1(:,:)
    
    complex(8), external :: zdotc

    allocate(vnl(nstfv,nstfv))
    vnl(:,:) = vxnl(:,:,ik)

!----------------------------------------
! Update the Hamiltonian
!----------------------------------------

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
    allocate(vnlmat(nmatp,nmatp))
    call zgemm('c','n',nmatp,nmatp,nstfv, &
   &  zone,temp,nstfv,temp1,nstfv,zzero,vnlmat,nmatp)
    deallocate(temp,temp1)
    
    if ((debug).and.(rank==0)) then
        call linmsg(fgw,'-',' Vx_NL_GG ')
        do ie1 = 1, nmatp, 20
            write(fgw,*) (vnlmat(ie1,ie2), ie2=1,nmatp,20)
        end do
        call linmsg(fgw,'-','')
    end if
    
! Update Hamiltonian
    system%hamilton%za(:,:) = system%hamilton%za(:,:)+ex_coef*vnlmat(:,:)
    deallocate(vnlmat)

    return
end subroutine
