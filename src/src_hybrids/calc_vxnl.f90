!
!BOP
! !ROUTINE: calc_vxnl
! !INTERFACE:
!
subroutine calc_vxnl
! !USES:    
    use modmain
    use modgw
    use modfvsystem
    use mod_hybrids
    use modmpi
! 
! !DESCRIPTION:
!   Calculates the non-local exchange potential 
!   and the non-local exchange energy for Hartree-Fock based hybrid functionals.
!
!EOP
!BOC
    implicit none

    integer(4) :: ikp, ik, jk, iq
    integer(4) :: ie1, ie2, ie3, icg, icg1
    integer(4) :: ist, l
    integer(4) :: i, j, k, jst, ispn
    integer(4) :: n, nmdim
    real(8)    :: tstart, tend
    complex(8) :: mvm, zt1
    
    real(8) :: evalfv(nstfv)
    complex(8), allocatable :: minm(:,:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8) :: vxfv(nstsv,nstsv), tmat(nstsv,nstsv)
    
    complex(8), external :: zdotc

    call cpu_time(tstart)
    if (rank == 0) call boxmsg(fgw,'-','Calculate Vx_NL')
    
    !------------------------------------------!
    ! Matrix elements of non-local potential   !
    !------------------------------------------!
    if (allocated(vxnl)) deallocate(vxnl)
    allocate(vxnl(nstfv,nstfv,nkpt))
    vxnl(:,:,:) = zzero
    
    if (iopcore<=1) then
      if (allocated(vxnlcc)) deallocate(vxnlcc)
      allocate(vxnlcc(ncg,nkpt))
      vxnlcc(:,:) = zzero
    end if 
    
    if (input%groundstate%tevecsv) then
      if (allocated(bxnl)) deallocate(bxnl)
      allocate(bxnl(nstsv,nstsv,nkpt))
      bxnl(:,:,:) = zzero
    end if
        
    ! Calculate the integration weights using the linearized tetrahedron method
    !if (allocated(evaldft)) deallocate(evaldft)
    !allocate(evaldft(nstsv,nkpt))
    !evaldft(:,:) = evalsv(:,:)
    !call kintw
    !deallocate(evaldft)
    !deallocate(kcw)
    !if (iopcore<2) deallocate(unw)
    
    ! Initialize mixed product basis
    call init_mb
    
!---------------------------------------
! Loop over k-points
!---------------------------------------

#ifdef MPI
    do ikp = firstk(rank), lastk(rank)
#else
    do ikp = 1, nkpt
#endif    
        ik = idikp(ikp)
       
!---------------------------------------
! Integration over BZ
!---------------------------------------
          do iq = 1, nqptnr
          
            Gamma = gammapoint(iq)

! Set the size of the basis for the corresponding q-point
            matsiz = locmatsiz+ngq(iq)
            if (rank==0) write(fgw,101) ikp, iq, locmatsiz, ngq(iq), matsiz

! Calculate the interstitial mixed basis functions
            call diagsgi(iq)

! Calculate the transformation matrix between pw's and the pw mixed basis functions
            call calcmpwipw(iq)

!------------------------------------               
! Calculate the bare Coulomb matrix
!------------------------------------
            call calcbarcmb(iq)
            call setbarcev(0.d0)

!------------------------------------------------------------
! (Re)-Calculate the M^i_{nm}(k,q) matrix elements for given k and q
!------------------------------------------------------------
            call expand_basis(ik,iq)

!------------------------------------------------------------     
! (Re)-Calculate non-local (Hartree-Fock) potential
!------------------------------------------------------------

! Valence electron contribution
            allocate(minm(mbsiz,nstfv,nstfv))
            nmdim = nstfv*nstfv
            call zgemm('c','n',mbsiz,nmdim,matsiz, &
            &          zone,barcvm,matsiz,minmmat,matsiz, &
            &          zzero,minm,mbsiz)
            deallocate(minmmat)

            jk = kqid(ik,iq)
            call getevalfv(vkl(:,indkp(jk)),evalfv)
            
            do ie1 = 1, nstfv
              do ie2 = ie1, nstfv
                zt1 = zzero
                do ie3 = 1, nstfv
                  if (evalfv(ie3)<=efermi) then
                    mvm = zdotc(mbsiz,minm(1:mbsiz,ie1,ie3),1,minm(1:mbsiz,ie2,ie3),1)
                    zt1 = zt1+mvm
                  end if
                end do ! ie3
                vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp)-zt1/dble(nqptnr)
              end do ! ie2
            end do ! ie1
            deallocate(minm)

! Core electron contribution
            if (iopcore<=1) then
                
              allocate(minm(mbsiz,nstfv,ncg))
              nmdim = nstfv*ncg
              call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
              &          zone,barcvm,matsiz,mincmat,locmatsiz, &
              &          zzero,minm,mbsiz)
              deallocate(mincmat)
                
              do ie1 = 1, nstfv
                do ie2 = ie1, nstfv
                  zt1 = zzero
                  do icg = 1, ncg
                    mvm = zdotc(mbsiz,minm(1:mbsiz,ie1,icg),1,minm(1:mbsiz,ie2,icg),1)
                    zt1 = zt1+mvm
                  enddo ! ie3
                  vxnl(ie1,ie2,ikp) = vxnl(ie1,ie2,ikp)-zt1/dble(nqptnr)
                end do ! ie2
              end do ! ie1
              deallocate(minm)
              
              !-----------------------------------------------------
              allocate(minm(mbsiz,ncg,nstfv))
              nmdim = ncg*nstfv
              call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
              &          zone,barcvm,matsiz,micmmat,locmatsiz, &
              &          zzero,minm,mbsiz)
              deallocate(micmmat)
              
              do icg = 1, ncg
                zt1 = zzero
                do ie2 = 1, nstfv
                  if (evalfv(ie2)<=efermi) then
                    mvm = zdotc(mbsiz,minm(1:mbsiz,icg,ie2),1,minm(1:mbsiz,icg,ie2),1)
                    zt1 = zt1+mvm
                  end if
                end do ! ie2
                vxnlcc(icg,ikp) = vxnlcc(icg,ikp)-zt1/dble(nqptnr)
              end do ! icg
              deallocate(minm)
              
              !-----------------------------------------------------
              allocate(minm(mbsiz,ncg,ncg))
              nmdim = ncg*ncg
              call zgemm('c','n',mbsiz,nmdim,locmatsiz, &
              &          zone,barcvm,matsiz,miccmat,locmatsiz, &
              &          zzero,minm,mbsiz)
              deallocate(miccmat)
              
              do icg = 1, ncg
                zt1 = zzero
                do icg1 = 1, ncg
                  mvm = zdotc(mbsiz,minm(1:mbsiz,icg,icg1),1,minm(1:mbsiz,icg,icg1),1)
                  zt1 = zt1+mvm
                end do
                vxnlcc(icg,ikp) = vxnlcc(icg,ikp)-zt1/dble(nqptnr)
              end do ! icg
              deallocate(minm)
                    
            end if ! iopcore

        end do ! iq
        
        do ie1 = 1, nstfv
          do ie2 = ie1+1, nstfv
            vxnl(ie2,ie1,ikp) = conjg(vxnl(ie1,ie2,ikp))
          end do
        end do

        !---------------------------------------------------------      
        ! compute the second-variational \Sigma_x matrix elements
        !---------------------------------------------------------
        if (input%groundstate%tevecsv) then

          allocate(evecsv(nstsv,nstsv))
          call getevecsv(vkl(:,ikp),evecsv)
          
          vxfv(:,:) = zzero
          vxfv(1:nstfv,1:nstfv) = vxnl(:,:,ikp)
          vxfv(nstfv+1:nstsv,nstfv+1:nstsv) = vxnl(:,:,ikp)

          call zgemm( 'n', 'n', &
          &           nstsv, nstsv, nstsv, &
          &           zone, &
          &           vxfv, nstsv, &
          &           evecsv, nstsv, &
          &           zzero, tmat, nstsv)

          call zgemm( 'c', 'n', &
          &           nstsv, nstsv, nstsv, &
          &           zone, &
          &           evecsv, nstsv, &
          &           tmat, nstsv, &
          &           zzero, bxnl(:,:,ikp), nstsv)
          
          bxnl(:,:,ikp) = bxnl(:,:,ikp)-vxfv(:,:)
          
          if (rank==0) write(*,*) 'bxup=', sum(bxnl(1:nstfv,1:nstfv,ikp))
          if (rank==0) write(*,*) 'bxdn=', sum(bxnl(nstfv+1:nstsv,nstfv+1:nstsv,ikp))
          
if (.false.) then
          do i = 1, nstsv
          do j = 1, nstsv
            k = 0
            do ispn = 1, nspinor
              do ist = 1, nstfv
                k = k+1
                l = (ispn-1)*nstfv
                do jst = 1, nstfv
                  l = l + 1
                  bxnl(:,:,ikp) = bxnl(:,:,ikp)+ &
                  &  conjg(evecsv(k,i))*evecsv(l,j)*vxnl(ist,jst,ikp)
                end do
              end do
            end do
          end do
          end do
          deallocate(evecsv)
          
          bxnl(1:nstfv,1:nstfv,ikp) = bxnl(1:nstfv,1:nstfv,ikp)-vxnl(:,:,ikp)
          bxnl(nstfv+1:nstsv,nstfv+1:nstsv,ikp) = &
          &  bxnl(nstfv+1:nstsv,nstfv+1:nstsv,ikp)-vxnl(:,:,ikp)
end if          

        end if ! sv        
        
        !------------------------------------------------------------
        ! Debugging Info
        !------------------------------------------------------------
        if ((debug).and.(rank==0)) then
            call linmsg(fgw,'-','')
            call linmsg(fgw,'-',' Diagonal elements of Vx_NL_nn ')
            write(fgw,*) 'for k-point ', ikp
            i = 0
            do ispn = 1, nspinor
            do ie1 = 1, nstfv
                if (input%groundstate%tevecsv) then
                  i = i+1
                  write(fgw,'(i4,4f12.4)') ie1, vxnl(ie1,ie1,ikp), bxnl(i,i,ikp)
                else
                  write(fgw,'(i4,2f12.4)') ie1, vxnl(ie1,ie1,ikp)
                end if
            end do
            end do
            if (iopcore<=1) then
              call linmsg(fgw,'-',' Diagonal elements of Vx_NL_CC ')
              do icg = 1, ncg
                write(fgw,'(i4,2f12.4)') icg, vxnlcc(icg,ikp)
              end do
            end if
        end if
    
    end do ! ikp

101 format(10x,/, &
    &       'Data for (k,q)-point:',2i4,//,10x,                  &
    &       'Mixed basis:',/,10x,                                &
    &       'Number of atomic basis functions:       ',i4,/,10x, &
    &       'Number of interstitial basis functions: ',i4,/,10x, &
    &       'Total number of basis functions:        ',i4,/)

#ifdef MPI
    call mpi_allgatherv_ifc(nkpt,nstfv*nstfv,zbuf=vxnl)
    call mpi_allgatherv_ifc(nkpt,ncg,zbuf=vxnlcc)
    if (input%groundstate%tevecsv) &
    &  call mpi_allgatherv_ifc(nkpt,nstsv*nstsv,zbuf=bxnl)
    call barrier
#endif

    !-----------------------------------------
    ! Calculate the non-local exchange energy
    !-----------------------------------------
    exnl = 0.d0
    do ikp = 1, nkpt
      call getevalfv(vkl(:,ikp),evalfv)
      do ie1 = 1, nstfv
        if (evalfv(ie1)<=efermi) then
          vxnl(ie1,ie1,ikp) = vxnl(ie1,ie1,ikp)- &
          &                   4.d0*pi*vi*singc2
          exnl = exnl+wkpt(ikp)*vxnl(ie1,ie1,ikp)
        end if
      end do
      if (iopcore<=1) then
        do icg = 1, ncg
          vxnlcc(icg,ikp) = vxnlcc(icg,ikp)- &
          &                 4.d0*pi*vi*singc2
          exnl = exnl+wkpt(ikp)*vxnlcc(icg,ikp)
        end do
      end if
    end do
    
    call cpu_time(tend)
    if (rank==0) call write_cputime(fgw,tend-tstart, 'CALC_VXNL')
    
    return
end subroutine
