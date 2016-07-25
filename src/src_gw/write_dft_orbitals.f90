
subroutine write_dft_orbitals
    use modinput
    use modmain
    use modgw
    use mod_mpi_gw, only : myrank
!    use mod_hdf5
    use m_getunit
    implicit none
    ! local variables
    integer :: ik, ib, ist, jst, ispn, i, j
    integer :: fid1, fid2, recl
    integer,    allocatable :: idx(:)
    real(8),    allocatable :: evalsv_(:)
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsv(:,:,:), evecsvt(:,:), evecsv_(:,:)
  
    if (allocated(evalsv)) deallocate(evalsv)
    allocate(evalsv(nstsv,nkpt))
    allocate(evecsv(nmatmax,nstsv,nspinor))
    if (input%groundstate%tevecsv) then
      allocate(evecfvt(nmatmax,nstfv))
      allocate(evecsvt(nstsv,nstsv))
      if (.not.ncmag) then
        allocate(evalsv_(nstsv))
        allocate(evecsv_(nmatmax,nstsv))
        allocate(idx(nstsv))
      end if
    end if
    if (allocated(spindex)) deallocate(spindex)
    allocate(spindex(nstsv,nkpt))
    do ib = 1, nstsv
      spindex(ib,:) = ib
    end do
    
!#ifndef _HDF5_
    ! Prepare Fortran binary files
    call getunit(fid1)
    inquire(IoLength=recl) nstsv, vkl(:,1), evalsv(:,1)
    open(fid1,File='GW_EVALSV.OUT',Form='UNFORMATTED',Access='DIRECT',Recl=recl)
    call getunit(fid2)
    inquire(IoLength=recl) nmatmax, nstsv, nspinor, vkl(:,1), evecsv
    open(fid2,File='GW_EVECSV.OUT',Form='UNFORMATTED',Access='DIRECT',Recl=recl)
!#endif

    ! write(*,*) 'nspinor=', nspinor
    ! write(*,*) 'ncmag=', ncmag

    do ik = 1, nkpt

      !-----------------
      ! eigenvalues 
      !-----------------
      evalsv = 0.d0
      call getevalsv(vkl(:,ik),evalsv(:,ik))
      
      !-----------------
      ! eigenvectors
      !-----------------
      evecsv = zzero
      if (input%groundstate%tevecsv) then
        
        evecfvt = zzero
        call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfvt)
        evecsvt = zzero
        call getevecsv(vkl(:,ik),evecsvt)
        
        do j = 1, nstsv
          i = 0
          do ispn = 1, nspinor
          do ist = 1, nstfv
            i = i+1
            call zaxpy(nmatmax, &
            &          evecsvt(i,j), &
            &          evecfvt(:,ist),1, &
            &          evecsv(:,j,ispn),1)
          end do
          end do
        end do
        
        if (.not.ncmag) then
          write(*,*) 'sort evalsv'
        !------------------------------------------------------
        ! Spin-collinear case: order in energy increasing order
        !------------------------------------------------------
          evalsv_(:) = evalsv(:,ik)
          call sortidx(nstsv,evalsv_,idx)
          do ispn = 1, nspinor
            evecsv_(:,:) = evecsv(:,:,ispn)
            do ib = 1, nstsv
              evalsv(ib,ik) = evalsv_(idx(ib))
              evecsv(:,ib,ispn) = evecsv_(:,idx(ib))
              spindex(ib,ik) = idx(ib)
            end do
          end do
        end if ! ncmag
        
      else
      
        call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecsv)
        
      end if
      
!#ifdef _HDF5_
!      write(cik,'(I4.4)') ik
!      if (.not.hdf5_exist_group(fgwh5,"/kpoints",cik)) &
!      &  call hdf5_create_group(fgwh5,"/kpoints",cik)
!      path = "/kpoints/"//trim(adjustl(cik))
!      call hdf5_write(fgwh5,path,"evalsv",evalsv(1,ik),(/nstsv/))
!      call hdf5_write(fgwh5,path,"evecsv",evecsv(1,1),(/nmatmax,nstsv/))
!#else
      ! store in Fortran binary
      write(fid1,Rec=ik) nstsv, vkl(:,ik), evalsv(:,ik)
      write(fid2,Rec=ik) nmatmax, nstsv, nspinor, vkl(:,ik), evecsv
!#endif

    end do ! ik
    
    deallocate(evalsv)
    deallocate(evecsv)
    if (input%groundstate%tevecsv) then
      deallocate(evecfvt)  
      deallocate(evecsvt)
      if (.not.ncmag) then
        deallocate(evalsv_)
        deallocate(evecsv_)
        deallocate(idx)
      end if
    end if
    
!#ifndef _HDF5_
    close(fid1)
    close(fid2)
!#endif

    return
end subroutine
