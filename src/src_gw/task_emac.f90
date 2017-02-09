
subroutine task_emac()

    use modinput
    use modmain
    use modgw
    use mod_mpi_gw
    use m_getunit
      
    implicit none
    integer :: iq, iom, iop, fid, im
    real(8) :: r0, r1, v(3)
    complex(8) :: e0, e1, e2
    complex(8), allocatable :: emac0_(:,:,:), emac1_(:,:,:), emac2_(:,:,:)
    complex(8), allocatable :: emac0(:,:,:), emac1(:,:,:), emac2(:,:,:)
#ifdef MPI
    integer :: ntype, ncnt, nsize
#endif
    
    if (rank==0) call boxmsg(fgw,'=','Calculate the macroscopic dielectric function')

    !==========================
    ! Perform initialization
    !==========================
    
    ! initialize local GW MPI environment
    call init_mpi_gw
    
    ! prepare GW global data
    call init_gw
    
    ! clear not used anymore global exciting variables
    call clean_gndstate
  
    if (.not.input%gw%rpmat) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    else
      write(*,*)
      write(*,*) 'Info(task_emac): Momentum matrix elements read from files.'
      write(*,*)
    end if
    
    ! occupancy dependent BZ integration weights
    call kintw
    
    !==========
    ! Set q=0
    !==========
    iq = 1
    Gamma = .true.
    
    !========================================
    ! Calculate interstitial basis functions
    !========================================
    matsiz = locmatsiz+Gqset%ngk(1,iq)
    call diagsgi(iq)
    call calcmpwipw(iq)
    
    !======================================
    ! Calculate the bare coulomb potential
    !======================================
    call calcbarcmb(iq)
    
    !========================================
    ! Set v-diagonal MB and reduce its size
    !========================================
    call setbarcev(input%gw%barecoul%barcevtol)
    call delete_coulomb_potential
    
    !===================================
    ! Calculate the dielectric function
    !===================================
    
#ifdef MPI
    call mpi_set_range(nproc_tot, &
    &                  myrank, &
    &                  freq%nomeg, 1, &
    &                  iomstart, iomend, &
    &                  iomcnt, iomdsp)
    !write(*,*) 'iom', myrank, iomstart, iomend
    !write(*,*) 'cnt', myrank, iomcnt(0:nproc_tot-1)
    !write(*,*) 'dsp', myrank, iomdsp(0:nproc_tot-1)
#else
    iomstart = 1
    iomend = freq%nomeg
#endif
    
    ! \epsilon_{00} without LFE
    allocate(emac0_(3,3,iomstart:iomend))
    emac0_ = 0.d0
    ! \epsilon_{00} including LFE
    allocate(emac1_(3,3,iomstart:iomend))
    emac1_ = 0.d0
    ! Averaged eps_{00}^{-1}-1
    allocate(emac2_(3,3,iomstart:iomend))
    emac2_ = 0.d0
    
    call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
    
    select case (trim(input%gw%scrcoul%scrtype))
    
      case('ppm','PPM')
        call calcepsilon_ppm(iq,iomstart,iomend)
        
      case default
        call calcepsilon(iq,iomstart,iomend)
        !---------------------------------
        ! Store \epsilon_{00} without LFE
        !---------------------------------
        do iom = iomstart, iomend
          emac0_(:,:,iom) = epsh(iom,:,:)
        end do ! iom
        !================================
        ! Invert the dielectric function
        !================================
        call calcinveps(iomstart,iomend)
        
    end select
    !---------------------------------
    ! \epsilon_{00} including LFE
    !---------------------------------
    do iom = iomstart, iomend
      emac1_(:,:,iom) = eps00(iom,:,:)
    end do ! iom
    
    !---------------------------------
    ! Averaged \epsilon_{00}^{-1}
    !---------------------------------
    do iom = iomstart, iomend
      emac2_(:,:,iom) = epsh(iom,1,1)
    end do
    
    ! clean unused data
    call delete_dielectric_function(Gamma)
    if (allocated(fnm)) deallocate(fnm)
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
    
    !================================
    ! Collect the results
    !================================
    allocate(emac0(3,3,freq%nomeg))
    allocate(emac1(3,3,freq%nomeg))
    allocate(emac2(3,3,freq%nomeg))
#ifdef MPI
    ncnt = iomend-iomstart+1
    nsize = 3*3
    call MPI_Type_Contiguous(nsize,MPI_DOUBLE_COMPLEX,ntype,ierr)
    call MPI_Type_Commit(ntype,ierr)
    call MPI_GatherV(emac0_,ncnt,ntype,emac0,iomcnt,iomdsp,ntype,0,mycomm,ierr)
    call MPI_GatherV(emac1_,ncnt,ntype,emac1,iomcnt,iomdsp,ntype,0,mycomm,ierr)
    call MPI_GatherV(emac2_,ncnt,ntype,emac2,iomcnt,iomdsp,ntype,0,mycomm,ierr)
    call MPI_Type_Free(ntype,ierr)
#else
    emac0 = emac0_
    emac1 = emac1_
    emac2 = emac2_
#endif
    deallocate(emac0_,emac1_,emac2_)
    
    if (myrank==0) then
      !---------------------------------
      ! Print \epsilon_{00} without LFE
      !---------------------------------
      call getunit(fid)
      open(fid, File='EPS00NLF.OUT', Form='Formatted', &
      &    Action='Write', Status='Replace')
      write(fid,*)
      write(fid,*) 'Macroscopic dieletric tensor (\eps_{00}, head) without local field effect'
      write(fid,*)
      do iom = 1, freq%nomeg
        write(fid, '(" frequency index and value: ",i6,f14.8)') iom, freq%freqs(iom)
        write(fid, '(" real part, imaginary part below")')
        write(fid, '(3f14.8,5x,3f14.8)') &
        &  (dble(emac0(iop,:,iom)), aimag(emac0(iop,:,iom)), iop = 1, 3)
        write(fid,*)
      end do ! iom
      close(fid)
      !-----------------------------------
      ! Print \epsilon_{00} including LFE
      !-----------------------------------
      call getunit(fid)
      open(fid, File='EPS00.OUT', Form='Formatted', &
      &    Action='Write', Status='Replace')
      write(fid,*)
      write(fid,*) 'Macroscopic dieletric tensor (\eps_{00}, head) including local field effect'
      write(fid,*)
      do iom = 1, freq%nomeg
        write(fid, '(" frequency index and value: ",i6,f14.8)') iom, freq%freqs(iom)
        write(fid, '(" real part, imaginary part below")')
        write(fid, '(3f14.8,5x,3f14.8)') &
        &  (dble(emac1(iop,:,iom)), aimag(emac1(iop,:,iom)), iop = 1, 3)
        write(fid,*)
      end do ! iom
      close(fid)
      !=================================
      ! Output the results of averaging
      !=================================
      write(fgw,*)
      write(fgw,*) 'Info(task_emac): '
      write(fgw,*) '  Averaged macroscopic dielectric function with and without local field'
      write(fgw,*)
      write(fgw,*)'# frequency    eps_{00} (diag)    eps_{00}+LFE (diag)    <eps_{00}^{-1}>'
      write(fgw,*)
      call getunit(fid)
      open(fid, File='EPSMACRO.OUT', Form='Formatted', &
      &    Action='Write', Status='Replace')
      write(fid,*)'# frequency    eps_{00} (diag)    eps_{00}+LFE (diag)    <eps_{00}^{-1}>'
      do iom = 1, freq%nomeg
        e0 = (emac0(1,1,iom)+emac0(2,2,iom)+emac0(3,3,iom)) / 3.d0
        e1 = (emac1(1,1,iom)+emac1(2,2,iom)+emac1(3,3,iom)) / 3.d0
        e2 = zone+emac2(1,1,iom)
        write(fgw,10) freq%freqs(iom), e0, e1, e2
        write(fid,10) freq%freqs(iom), e0, e1, e2
      end do
      close(fid)
      10 format(7G12.4)
    end if ! myrank
    deallocate(emac0,emac1,emac2)
    
    return
end subroutine
