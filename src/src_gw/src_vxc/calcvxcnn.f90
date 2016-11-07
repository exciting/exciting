!BOP
!
!!ROUTINE: calcvxcnn
!
!!INTERFACE:
!
subroutine calcvxcnn
!
!!DESCRIPTION:
!
! This subroutine calculates the diagonal matrix elements of 
! the exchange correlation potential (only for valence states).

!!USES:
    use modinput
    use modmain, only : apwordmax, lmmaxapw, lmmaxvr, natmtot, nlomax, &
    &                   nstfv, nspinor, nstsv, nmatmax, nspecies, zzero, &
    &                   nmat, natoms, vxcmt, vxcir, zone
    use modgw
    use modmpi
    use m_getunit
    use mod_hdf5

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ik, ik0
    integer(4) :: ist, jst, i, j, k, l, ispn
    integer(4) :: ia, is
    integer(4) :: ngp, fid
    real(8)    :: tstart, tend
    complex(8) :: zsum, zt1
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:), evecsv(:,:), evec(:,:,:)
    complex(8), allocatable :: h(:)
    real(8)    :: t0, t1

!!EXTERNAL ROUTINES: 
    complex(8), external :: zdotc

!!REVISION HISTORY:
! 
! Created    8th. Aug. 2006 by RGA
! Revisited       May  2011 by DIN
!
!EOP
!BOC
    call timesec(tstart)
    
    ! Global array to store <n|Vxc|n>
    if (allocated(vxcnn)) deallocate(vxcnn)
    allocate(vxcnn(nstsv,kset%nkpt))
    vxcnn = zzero
    
    ! allocate exchange-correlation integral arrays
    if (allocated(vxcraa)) deallocate(vxcraa)
    allocate(vxcraa(apwordmax, &
    &               0:input%groundstate%lmaxmat, &
    &               apwordmax, &
    &               0:input%groundstate%lmaxapw, &
    &               0:lmmaxvr, &
    &               natmtot))
    if (allocated(vxcrloa)) deallocate(vxcrloa)
    allocate(vxcrloa(nlomax, &
    &                apwordmax, &
    &                0:input%groundstate%lmaxmat, &
    &                0:lmmaxvr, natmtot))
    if (allocated(vxcrlolo)) deallocate(vxcrlolo)
    allocate(vxcrlolo(nlomax, &
    &                 nlomax, &
    &                 0:lmmaxvr, &
    &                 natmtot))

    ! Calculate radial integrals
    call vxcrad
    
    ! Fourier transform the interstitial part of Vxc
    call genvxcig
    
    allocate(apwalm(Gkset%ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evec(nmatmax,nstsv,nspinor))
    allocate(h(nmatmax))
    
    do ik = firstofset(rank,kset%nkpt), lastofset(rank,kset%nkpt)
        
      ik0 = kset%ikp2ik(ik)
      ngp = Gkset%ngk(1,ik0)

!#ifdef _HDF5_      
!      ! get the eigenvectors from file
!      write(cik,'(I4.4)') ik0
!      path = "/kpoints/"//trim(adjustl(cik))
!      call hdf5_read(fgwh5,path,"evecsv",evecsv(1,1),(/nmatmax,nstsv/))
!#else
      call getevecsvgw_new('GW_EVECSV.OUT',ik0,kqset%vkl(:,ik0),nmatmax,nstsv,nspinor,evec)
!#endif

      ! find the matching coefficients
      call match(ngp, &
      &          Gkset%gkc(:,1,ik0), &
      &          Gkset%tpgkc(:,:,1,ik0), &
      &          Gkset%sfacgk(:,:,1,ik0), &
      &          apwalm)

      do ispn = 1, nspinor
        do i = ibgw, nbgw
          h(:) = zzero
          ! muffin-tin contributions
          do is = 1, nspecies
          do ia = 1, natoms(is)
            call vxcaa(is,ia,ngp,apwalm,evec(:,i,ispn),h)
            call vxcalo(is,ia,ngp,apwalm,evec(:,i,ispn),h)
            call vxclolo(is,ia,ngp,evec(:,i,ispn),h)
          end do
          end do
          ! interstitial contribution
          call vxcistl(ngp,Gkset%igkig(:,1,ik0),evec(:,i,ispn),h)
          vxcnn(i,ik) = vxcnn(i,ik)+zdotc(nmat(1,ik),evec(:,i,ispn),1,h,1)
        end do ! i
      end do ! ispn
      
    enddo ! ik
    
    deallocate(h)
    deallocate(apwalm)
    deallocate(evec)
    deallocate(vxcraa)
    deallocate(vxcrloa)
    deallocate(vxcrlolo)
    
#ifdef MPI
    call mpi_allgatherv_ifc(kset%nkpt,rlen=nstsv,zbuf=vxcnn)
    call barrier
#endif

    if (rank==0) then
#ifdef _HDF5_
      do ik = 1, kset%nkpt
        write(cik,'(I4.4)') ik
        path = "/kpoints/"//trim(adjustl(cik))
        if (.not.hdf5_exist_group(fgwh5,"/kpoints",cik)) &
        &  call hdf5_create_group(fgwh5,"/kpoints",cik)
        call hdf5_write(fgwh5,path,"vxcnn",vxcnn(1,ik),(/nbandsgw/))
      end do
#else    
    ! print results into file VXCNN.OUT
      call getunit(fid)
      open(fid,file='VXCNN.OUT',form='FORMATTED',status='UNKNOWN')
      do ik = 1, kset%nkpt
        write(fid,'("ik=",i4,"    vkl=",3f8.4)') ik, kset%vkl(:,ik)
        do i = ibgw, nbgw
          write(fid,'(i4,2f12.4)') i, vxcnn(i,ik)
        end do
        write(fid,*)
      end do
      close(fid)
#endif
    end if
      
    call timesec(tend)
    time_vxc = time_vxc + (tend-tstart)
    
end subroutine
!EOC


