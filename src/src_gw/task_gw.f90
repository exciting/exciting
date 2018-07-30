!BOP
!
!!ROUTINE: \verb"task_gw"
!
!!INTERFACE:
!      
subroutine task_gw()
!      
!!DESCRIPTION:
!
! This subroutine performs one GW cycle and calculates the corresponding
! quasiparticle energies.
!
!!USES:
    use modinput
    use modmain,               only : zzero, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
            
!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ikp, iq, fid, ik
    real(8)    :: t0, t1
    integer(4) :: recl
    real(8)    :: ab_plane, ab_norm(3), q0_vol
    integer    :: im
    complex(8) :: vc

!!REVISION HISTORY:
!
! Created Nov 2013 by (DIN)
!
!EOP
!BOC
    !===========================================================================
    ! Initialization
    !===========================================================================
    !----------------------------------------------
    ! Store all important results to the hdf5 file
    !----------------------------------------------
#ifdef _HDF5_
    call hdf5_initialize()
    fgwh5 = "gw_output.h5"
    if (rank==0) then
      call hdf5_create_file(fgwh5)
      call hdf5_create_group(fgwh5,"/","parameters")
      if (rank==0) call write_gw_parameters_hdf5
      call hdf5_create_group(fgwh5,"/","kpoints")
    end if
#endif    
    
    ! prepare GW global data
    call init_gw
    
    !=================================================
    ! Calculate the diagonal matrix elements of the 
    ! DFT exchange-correlation potential
    !=================================================
    ! it is better to do it here to deallocate cfunir and vxcir arrays
    call timesec(t0)
    call calcvxcnn
    call timesec(t1)

    ! clean not used anymore global exciting variables
    call clean_gndstate
    
    if (input%gw%taskname.ne.'g0w0_x') then
      if (.not.input%gw%rpmat) then
        !========================================================
        ! calculate momentum matrix elements and store to a file
        !========================================================
        call calcpmatgw
      end if
    end if
    
    ! occupancy dependent BZ integration weights
    call kintw
    
    !---------------------------------------
    ! treatment of singularities at G+q->0
    !---------------------------------------
    singc1 = 0.d0
    singc2 = 0.d0
    
    if (vccut) then
    
      select case (trim(input%gw%barecoul%cutofftype))
      
        case('0d')
          rccut = 0.5d0*dsqrt(dot_product(avec(:,3),avec(:,3)))
          i_sz = twopi*rccut**2
    
        case('2d')
          !--------------------------------------------------------
          ! Spherically averaged value of the integral around q->0
          !--------------------------------------------------------
          ! cutoff length
          rccut = 0.5d0*dsqrt(dot_product(avec(:,3),avec(:,3)))
          ! ab-plane surface area
          call r3cross(avec(:,1),avec(:,2),ab_norm(:))
          ab_plane = dsqrt(dot_product(ab_norm(:),ab_norm(:)))
          q0_vol   = twopi/dsqrt(pi*ab_plane*kqset%nkpt)
          i_sz = q0_vol*rccut - ((q0_vol*rccut)**2.0d0)/4.0d0
          i_sz = 2.d0*ab_plane*kqset%nkpt*i_sz
          ! if (myrank) then
          !   write(*,*)
          !   write(*,*) 'LIMIT q->0'
          !   write(*,*) ' nqpt = ', kqset%nkpt
          !   write(*,*) ' rccut = ', rccut
          !   write(*,*) ' ab_norm = ', ab_norm
          !   write(*,*) ' ab_plane = ', ab_plane
          !   write(*,*) ' q0_vol = ', q0_vol
          !   write(*,*) ' i_sz = ', i_sz
          !   write(*,*)
          ! end if
      
        case default
          write(*,*) 'ERROR(task_gw): Specified cutoff type is not implemented!'
        
      end select
        
    else
      
        select case (trim(input%gw%selfenergy%singularity))
          case('none')
            singc1 = 0.d0
            singc2 = 0.d0
          case('mpb')
            ! Auxiliary function method
            call setsingc
          case('crg')  
            ! Auxiliary function method
            call calc_q0_singularities
          case default
            write(*,*) 'ERROR(task_gw): Unknown singularity treatment scheme!'
            stop
        end select

    end if
    
    ! initialize self-energy arrays
    call init_selfenergy(ibgw,nbgw,kset%nkpt,freq%nomeg)

    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================    

#ifdef MPI
    call set_mpi_group(kqset%nkpt)
    call mpi_set_range(nproc_row, &
    &                  myrank_row, &
    &                  kqset%nkpt, 1, &
    &                  iqstart, iqend)
    call mpi_set_range(nproc_col, &
    &                  myrank_col, &
    &                  freq%nomeg, 1, &
    &                  iomstart, iomend, &
    &                  iomcnt, iomdsp)
    ! write(*,*) "myrank_row, iqstart, iqend =", myrank_row, iqstart, iqend
    ! write(*,*) "myrank_col, iomstart, iomend =", myrank_col, iomstart, iomend
    ! write(*,*) 'iomcnt: ', iomcnt(0:nproc_col-1)
    ! write(*,*) 'iomdsp: ', iomdsp(0:nproc_col-1)
#else
    iqstart = 1
    iqend = kqset%nkpt
    iomstart = 1
    iomend = freq%nomeg
#endif

    if (myrank==0) then
      call boxmsg(fgw,'=','GW cycle')
      call flushifc(fgw)
    end if
    
    ! each process does a subset
    do iq = iqstart, iqend
  
      if ((input%gw%debug).and.(rank==0)) then
        write(fdebug,*) '(task_gw): q-point cycle, iq = ', iq
      end if
    
      Gamma = gammapoint(kqset%vqc(:,iq))
          
      !========================================
      ! Calculate interstitial basis functions
      !========================================
      matsiz = locmatsiz+Gqset%ngk(1,iq)
      call diagsgi(iq)
      call calcmpwipw(iq)
    
      !======================================
      ! Calculate the bare Coulomb potential
      !======================================
      call calcbarcmb(iq)

      !===============================
      ! Calculate \Sigma^{x}_{kn}(q)
      !===============================
      call calcselfx(iq)
      
      if (input%gw%selfenergy%secordw) then
        call calcselfcSR(iq,iomstart,iomend)
      end if
    
      if (input%gw%taskname.ne.'g0w0_x') then
        !========================================
        ! Set v-diagonal MB and reduce its size
        !========================================
        if (vccut) then
          mbsiz = matsiz
          if (allocated(barc)) deallocate(barc)
          allocate(barc(matsiz,mbsiz))
          do im = 1, matsiz
            vc = cmplx(barcev(im),0.d0,8)
            barc(:,im) = vmat(:,im)*sqrt(vc)
          end do
        else
          call setbarcev(input%gw%barecoul%barcevtol)
        end if
        call delete_coulomb_potential
        !===================================
        ! Calculate the dielectric function
        !===================================
        call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
        select case (trim(input%gw%scrcoul%scrtype))
          case('ppm','PPM')
            call calcepsilon_ppm(iq,iomstart,iomend)
          case default
            call calcepsilon(iq,iomstart,iomend)
            !==========================================
            ! Calculate the screened Coulomb potential
            !==========================================
            call calcinveps(iomstart,iomend)
        end select
        !========================================
        ! Calculate the q-dependent self-energy
        !========================================
        call calcselfc(iq)
        call delete_dielectric_function(Gamma)
        if (allocated(kcw)) deallocate(kcw)
        if (allocated(unw)) deallocate(unw)
      end if
      
      ! clean unused data
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)
      
    end do ! iq
    
    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)
    
#ifdef MPI
    if ((nproc_row>1).and.(myrank_col==0)) then
      write(*,*) "sum self-energy from different q-points"
      call mpi_sum_array(0,selfex,nbandsgw,kset%nkpt,mycomm_row)
      write(*,*) "sum selfex done"  
      if (input%gw%taskname.ne.'g0w0_x') then
        ! G0W0 and GW0 approximations
        call mpi_sum_array(0,selfec,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
        if (input%gw%selfenergy%secordw) then
          call mpi_sum_array(0,selfecw2,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
          call mpi_sum_array(0,selfecSR,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
        end if
        if (input%gw%taskname=='cohsex') then
          ! COHSEX
          call mpi_sum_array(0,sigsx,nbandsgw,kset%nkpt,mycomm_row)
          call mpi_sum_array(0,sigch,nbandsgw,kset%nkpt,mycomm_row)
        end if ! cohsex
        write(*,*) "sum selfec done"  
      end if ! selfec
    endif
#endif
    
    if (myrank==0) then
    
      !===============================
      ! Write self-energies to files
      !===============================
      call timesec(t0)
#ifdef _HDF5_
      do ik = 1, kset%nkpt
        write(cik,'(I4.4)') ik
        path = "/kpoints/"//trim(adjustl(cik))
        if (.not.hdf5_exist_group(fgwh5,"/kpoints",cik)) &
        &  call hdf5_create_group(fgwh5,"/kpoints",cik)
        ! k-point
        call hdf5_write(fgwh5,path,"vkl",kset%vkl(1,ik),(/3/))
        ! exchange
        call hdf5_write(fgwh5,path,"selfex",selfex(1,ik),(/nbandsgw/))
        ! correlation
        if (input%gw%taskname.ne.'g0w0_x') &
        &  call hdf5_write(fgwh5,path,"selfec",selfec(1,1,ik),(/nbandsgw,freq%nomeg/))
      end do
#else
      call write_selfenergy(ibgw,nbgw,kset%nkpt,freq%nomeg)
#endif
      call timesec(t1)
      time_io = time_io+t1-t0
      
      ! KS band structure
      evalks(ibgw:nbgw,:) = evalsv(ibgw:nbgw,:) 
      !call bandanalysis('KS',ibgw,nbgw,evalks(ibgw:nbgw,:),efermi)
      call bandstructure_analysis('KS', &
      &  ibgw,nbgw,kset%nkpt,evalks(ibgw:nbgw,:),efermi)
      
      !=======================================
      ! Calculate the quasiparticle energies
      !=======================================
      call calcevalqp
      
      !------------------------------------------------------
      ! Write quasi-particle energies to file
      !------------------------------------------------------
      call timesec(t0)
      call write_qp_energies('EVALQP.DAT')
      call timesec(t1)
      time_io = time_io+t1-t0
      
      ! G0W0 QP band structure
      select case (input%gw%taskname)
      
        case('g0w0_x')
          !call bandanalysis('G0W0_X',ibgw,nbgw,evalqp(ibgw:nbgw,:),eferqp)
          call bandstructure_analysis('G0W0_X', &
          &  ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)
      
        case('cohsex')
          !call bandanalysis('COHSEX',ibgw,nbgw,evalqp(ibgw:nbgw,:),eferqp)
          call bandstructure_analysis('COHSEX', &
          &  ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)
          
        case('g0w0','gw0')
          !call bandanalysis('G0W0',ibgw,nbgw,evalqp(ibgw:nbgw,:),eferqp)
          call bandstructure_analysis('G0W0', &
          &  ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)
          
      end select
      
      if (input%gw%selfenergy%secordw) then 
      !---------------------------------
      ! Second order screened exchange
      !---------------------------------
        call selfcGW_vs_selfcW2
      end if
 
    end if ! myrank
    
    !--------------------------------------------------------
    ! Calculate quasiparticle energies in GW0 approximation 
    !--------------------------------------------------------
    if (input%gw%taskname=='gw0') then
      
      ! self-consistent cycle
      call calcscgw0
      
      ! print GW0 QP band structure
      if (myrank==0) then
        call timesec(t0)
        !----------------------------------------
        ! Write quasi-particle energies to file
        !----------------------------------------
        call write_qp_energies('EVALQP-GW0.DAT')
        !call bandanalysis('GW0',ibgw,nbgw,evalqp(ibgw:nbgw,:),eferqp)
        call bandstructure_analysis('GW0', &
        &  ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)
        call timesec(t1)
        time_io = time_io+t1-t0
      end if
    end if
    
    if (myrank==0) then
      !----------------------------------------
      ! Save QP energies into binary file
      !----------------------------------------
      call timesec(t0)
      call putevalqp()
#ifdef _HDF5_
      call hdf5_write(fgwh5,"/","efermi",efermi)
      call hdf5_write(fgwh5,"/","eferqp",eferqp)
      do ik = 1, kset%nkpt
        write(cik,'(I4.4)') ik
        path = "/kpoints/"//trim(adjustl(cik))
        ! KS energies
        call hdf5_write(fgwh5,path,"evalks",evalks(1,ik),(/nbandsgw/))
        ! QP energies
        call hdf5_write(fgwh5,path,"evalqp",evalqp(1,ik),(/nbandsgw/))
      end do
#endif
    end if ! myrank
    
    if (allocated(evalsv)) deallocate(evalsv)
    call delete_selfenergy
    
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine
!EOC      
