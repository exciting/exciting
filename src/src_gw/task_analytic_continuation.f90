
subroutine task_analytic_continuation()

    use modinput
    use modmain
    use modgw
    use mod_frequency
    use mod_hdf5
    use mod_mpi_gw
    use m_getunit
    implicit none
    ! local variables
    integer :: ik, ik_, ie, ie_, fid, recl
    real(8) :: egap
    character(20) :: s1, s2
       
    call init0
    call init1
    
    nvelgw = chgval-occmax*dble(ibgw-1)
    nbandsgw = nbgw-ibgw+1
    call init_kqpoint_set
    call generate_freqgrid(freq, &
    &                      input%gw%freqgrid%fgrid, &
    &                      input%gw%freqgrid%fconv, &
    &                      input%gw%freqgrid%nomeg, &
    &                      input%gw%freqgrid%freqmax)
    
    if (myrank==0) then
    
      ! allocate the arrays
      if (allocated(evalsv)) deallocate(evalsv)
      allocate(evalsv(nstsv,kset%nkpt))
      allocate(vxcnn(ibgw:nbgw,kset%nkpt))
      call init_selfenergy(ibgw,nbgw,kset%nkpt,freq%nomeg)
    
      ! read data from files
#ifdef _HDF5_
      !call hdf5_read(fgwh5,"/","efermi",efermi)
      do ik = 1, kset%nkpt
        write(cik,'(I4.4)') ik
        path = "/kpoints/"//trim(adjustl(cik))
        call hdf5_read(fgwh5,path,"evalks",evalks(ibgw,ik),(/nbandsgw/))
        evalsv(:,ik) = evalks(:,ik)
        call hdf5_read(fgwh5,path,"vxcnn",vxcnn(ibgw,ik),(/nbandsgw/))
        call hdf5_read(fgwh5,path,"selfex",selfex(ibgw,ik),(/nbandsgw/))
        call hdf5_read(fgwh5,path,"selfec",selfec(ibgw,1,ik),(/nbandsgw,freq%nomeg/))
      end do ! ik
#else
      do ik = 1, kset%nkpt
        ik_ = kset%ikp2ik(ik)
        call getevalsvgw_new('GW_EVALSV.OUT',ik_,kqset%vkl(:,ik_), &
        &                     nstsv,evalsv(1,ik))
        evalks(ibgw:nbgw,ik) = evalsv(ibgw:nbgw,ik)
      end do
      call getunit(fid)
      open(fid,file='VXCNN.OUT',form='FORMATTED',status='OLD',action='READ')
      do ik = 1, kset%nkpt
        read(fid,*) s1, ik_, s2, kset%vkl(:,ik)
        do ie = ibgw, nbgw
          read(fid,*) ie_, vxcnn(ie,ik)
        end do
        read(fid,*)
      end do
      close(fid)
      call readselfx
      call readselfc
#endif

      ! KS states analysis
      call fermi_exciting(input%groundstate%tevecsv, &
      &                   nvelgw, &
      &                   nbandsgw,kset%nkpt,evalks(ibgw:nbgw,:), &
      &                   kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
      &                   efermi,egap)
      call bandstructure_analysis('KS', &
      &  ibgw,nbgw,kset%nkpt,evalks(ibgw:nbgw,:),efermi)
    
      !======================================
      ! Calculate the quasiparticle energies
      !======================================
      call calcevalqp
      
      !------------------------------------------------------
      ! Write quasi-particle energies to file
      !------------------------------------------------------
      call write_qp_energies('EVALQP-AC.DAT')
      call bandstructure_analysis('G0W0-AC',ibgw,nbgw,kset%nkpt,&
      &                            evalqp(ibgw:nbgw,:),eferqp)
 
      !----------------------------------------
      ! Save QP energies into binary file
      !----------------------------------------
      call putevalqp()
      
      ! clear memory
      deallocate(evalsv)
      deallocate(vxcnn)
      call delete_selfenergy
      
    end if ! myrank
      
    return
end subroutine
