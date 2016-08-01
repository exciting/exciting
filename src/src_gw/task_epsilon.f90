      
subroutine task_epsilon

    use modinput
    use modmain,               only : zzero, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
            
    implicit none
    integer(4) :: ikp, iq, fid, ik
    real(8)    :: t0, t1
    integer(4) :: recl
    integer :: im
    
    ! mapping array 
    integer, allocatable :: iq2rank(:)
    
    complex(8), allocatable :: epsilon_(:,:,:,:)

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
    
    ! clean not used anymore global exciting variables
    call clean_gndstate
    
    if (.not.input%gw%rpmat) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    end if
    
    ! occupancy dependent BZ integration weights
    call kintw
    
    !===========================================================================
    ! Main loop: BZ integration
    !===========================================================================    

#ifdef MPI
    call set_mpi_group(kqset%nkpt)
    call mpi_set_range(nproc_row, &
    &                  myrank_row, &
    &                  kset%nkpt, 1, &
    &                  iqstart, iqend)
#else
    iqstart = 1
    iqend = kqset%nkpt
#endif
    iomstart = 1
    iomend = freq%nomeg
    
    ! iq <--> MPI rank mapping
    allocate(iq2rank(kset%nkpt))
    iq2rank = -1
    do iq = iqstart, iqend
      iq2rank(iq) = rank
    end do
    
    if (allocated(epsilon_)) deallocate(epsilon_)
    allocate(epsilon_(matsizmax,matsizmax,iomstart:iomend,iqstart:iqend))
    epsilon_(:,:,:,:) = 0.d0

    if (myrank==0) then
      call boxmsg(fgw,'=','q-point cycle')
      call flushifc(fgw)
    end if
    
    ! each process does a subset
    do iq = iqstart, iqend
    
      write(*,*)
      write(*,*) '(task_epsilon): q-point cycle, iq = ', iq
    
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
      
      !========================================
      ! Set v-diagonal MB and reduce its size
      !========================================
      call setbarcev(input%gw%barecoul%barcevtol)
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
      end select
      
      ! save to local array
      epsilon_(1:mbsiz,1:mbsiz,iomstart:iomend,iq) = epsilon(:,:,:)
      
      call delete_dielectric_function(Gamma)
      if (allocated(kcw)) deallocate(kcw)
      if (allocated(unw)) deallocate(unw)
      
      ! clean unused data
      if (allocated(mpwipw)) deallocate(mpwipw)
      if (allocated(barc)) deallocate(barc)
      
    end do ! iq
    
    if (allocated(kiw)) deallocate(kiw)
    if (allocated(ciw)) deallocate(ciw)
    
    !===================================
    ! Write dielectric function to file
    !===================================
    
    ! overwrite existing files
    if (rank==0) then
      call getunit(fid)
      open(fid,File='EPSILON.OUT',form='UNFORMATTED',status='REPLACE')
      close(fid)
    endif
    call barrier
    
    do iq = 1, kset%nkpt
      if (rank==iq2rank(iq)) then
#ifdef _HDF5_
        write(cik,'(I4.4)') ik
        path = "/qpoints/"//trim(adjustl(cik))
        if (.not.hdf5_exist_group(fgwh5,"/qpoints",cik)) &
        &  call hdf5_create_group(fgwh5,"/qpoints",cik)
        call hdf5_write(fgwh5,path,"epsilon", &
        &               epsilon_(1,1,1,iq),(/matsizmax,matsizmax,1:freq%nomeg/))
#endif
        call getunit(fid)
        inquire(iolength=recl) epsilon_(:,:,:,iq)
        open(fid,File="EPSILON.OUT",Action='WRITE',Form='UNFORMATTED',&
        &    Access='DIRECT',Status='OLD',Recl=recl)
        write(fid,rec=iq) epsilon_(:,:,:,iq)
        close(fid)
      end if ! rank
      call barrier
    end do ! iq      
    
    deallocate(epsilon_)    
    if (allocated(evalsv)) deallocate(evalsv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine

