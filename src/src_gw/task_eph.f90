!BOP
!
!!ROUTINE: \verb"task_gw"
!
!!INTERFACE:
!      
subroutine task_eph()
!      
!!DESCRIPTION:
!
! This subroutine performs one GW cycle and calculates the corresponding
! quasiparticle energies.
!
!!USES:
  Use modinput
  Use modmain
  use modgw
  use mod_mpi_gw
  use mod_wannier
  Use FoX_wxml
  use m_wsweight
  use m_getunit
  use m_plotmat
  use m_wannier_interpolate 
  use mod_symmetry
  use mod_dynmat 
  use mod_constants, only: zzero
            
!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ikp, iq, fid, ik, iw, ib, ik1
    real(8)    :: t0, t1
    !integer(4) :: recl
    real(8)    :: ab_plane, ab_norm(3), q0_vol, vtmp(3)
    integer    :: im, iknr, imode
    complex(8) :: vc
    real   (8) :: efold, evtmp, cbmold, evtmp1

    real(8), allocatable :: eval1     (:,:) ! original eigenvalues on the KS mesh 
    real(8), allocatable :: eval2     (:,:) ! Wannier interpolated eigenvalues 
    real(8), allocatable :: evalpath  (:,:) ! Wannier interpolated eigenvalues on a path 
    real(8), allocatable :: evalfv    (:,:) ! auxiliary eigenvalue storage 
    real(8), allocatable :: phfreq    (:,:) ! original eigenvalues on the KS mesh 

    logical    :: l_speceph, l_sigmaeph
    integer    :: nqgrid
    integer    :: ngridkq(3)
!    integer    :: ngridkqtot
    
  !--------------------------------------------------!      
  ! Calculate bandstructure by Wannier interpolation !
  !--------------------------------------------------!
  if( .not. associated( input%properties%wannier)) then
    write(*,*) " Error (bandstr): Wannier functions have not been calculated."
    call terminate
  end if
  call readfermi                !saves fermi energy in variable 'efermi'
  

  !--------------------------------------------
  ! Get the KS eigenvalues on the coarse k grid 
  !--------------------------------------------
  allocate( eval1( nstfv, nkptnr), evalfv( nstfv, nspinor))
  do iknr = 1, wf_kset%nkpt
    call getevalfv( wf_kset%vkl( :, iknr), evalfv)
    eval1( :, iknr) = evalfv( :, 1)
  enddo

  !----------------------------------------------
  ! Assign the variables read from the input file 
  !----------------------------------------------
  ibeph=input%eph%ibeph
  nbeph=input%eph%nbeph
  ibsumeph=input%eph%ibsumeph
  nbsumeph=input%eph%nbsumeph
  nomegeph=input%eph%freqgrideph%nomegeph
  ngridkq=input%eph%ngridqeph
  ngridkqtot=ngridkq(1)*ngridkq(2)*ngridkq(3)
  !print*, "input%eph%ibeph",                input%eph%ibeph 
  !print*, "input%eph%nbeph",                input%eph%nbeph 
  !print*, "input%eph%freqgrideqh%nomegeph", input%eph%freqgrideph%nomegeph

  !----------------------------
  ! Generate the frequency grid 
  !----------------------------
  call generate_freqgrid(freq, &
  &                      'eqdist', &
  &                      'refreq', &
  &                      input%eph%freqgrideph%nomegeph, &
  &                      input%eph%freqgrideph%freqmaxeph)
  do iw = 1, nomegeph 
    freq%freqs(iw)= freq%freqs(iw) - input%eph%freqgrideph%freqmaxeph / 2 
  enddo 
  if (.false. .and. rank==0) call print_freqgrid(freq,fgw) 
  
  !-----------------------
  ! Generate the G vectors
  !-----------------------
  call generate_G_vectors(gset,  &
  &                       bvec,  &
  &                       intgv, &
  &                       input%groundstate%gmaxvr)
  if ((.false.).and.(rank==0)) call print_G_vectors(gset,fid)

  !------------------------------------------
  ! Generate the dense q mesh for the phonons
  !------------------------------------------
  call generate_k_vectors(qsetd,bvec,ngridkq, &
                         input%eph%vqloffeph,.false.,.false.) ! generate a dense reducible k-point mesh
  if (.true.) then
    call getunit(fid)
    open(fid,file='PH_QPOINTS.OUT',action='Write',status='Unknown')
    call print_k_vectors(qsetd,fid)
  endif
  !
  !test the transformation from cartesian to lattice coordinates 
  !do iq = 1, ngridkqtot 
  !  call r3mv  (binv, qsetd%vkcnr (:,iq), vtmp )  
  !  write(6,104), iq, qsetd%vkcnr(:,iq), qsetd%vklnr(:,iq), vtmp
  !enddo
  !104 format(i4,3x,3f8.4,4x,3f8.4,4x,3f8.4)
  
  !-----------------------------------------------
  ! Generate the k-point path for final quantities 
  !-----------------------------------------------
  input%properties%bandstructure%wannier = .FALSE.
  task=20
  call init1        ! initialize k-points along path
  !
  if (.false.)then
    print*, "nkpt path = " , nkpt 
    do ik = 1 , nkpt
      print*, vkl(:, ik)
    enddo
    print*,  eval1 
  endif

  !------------------------------------
  ! Read the dynamical matrix from file 
  !------------------------------------
  call read_dyn()
  !
  !------------------------------------------------
  ! convert the reducible coarse grid from cartesian 
  ! to lattice coord (read from the dynamical matrix files) 
  !------------------------------------------------
  if(allocated(xqlred)) deallocate(xqlred) 
  allocate(xqlred(3,nqredtot)) 
  do iq =1, nqredtot
    call r3mv (binv, xqcred(:,iq), xqlred(:,iq)) 
    write(6,103), iq, xqcred (:,iq), xqlred(:,iq)
  enddo
  103 format(i4,3x,3f8.4,4x,3f8.4)
   

  if (.false.)then 
    ! interpolation on the dense q-mesh on the whole BZ
    if(allocated(dynmatD)) deallocate(dynmatD) 
    allocate(dynmatD(ndynmat,ndynmat,ngridkqtot)) 
    dynmatD=zzero
    CALL wann_ph (dynmat, ndynmat, nqirr, nqredtot, xqcirr, xqirrdeg, dynmatD, ngridkqtot, qsetd%vkcnr, nqred)
    print*, ' Wannier interpolation of the phonons: DONE! '
  else 
    ! interpolation on a path 
    if(allocated(dynmatD)) deallocate(dynmatD) 
    allocate(dynmatD(ndynmat,ndynmat,nkpt)) 
    dynmatD=zzero
    CALL wann_ph (dynmat, ndynmat, nqirr, nqredtot, xqcirr, & 
                  xqirrdeg, dynmatD, nkpt, vkcnr, nqred)
    print*, ' Wannier interpolation of the phonons: DONE! '

    allocate(phfreq(ndynmat,nkpt)) 
    do iq = 1 , nkpt
      print*, 'cdiagh2 for iq = ', iq
      call cdiagh2(ndynmat,dynmatD(:,:,iq),ndynmat,phfreq(:,iq))
      do imode = 1, ndynmat 
        print*, iq, imode, phfreq(imode,iq)
      enddo
    enddo
  endif 

  stop 
  ! interpolate KS eigenvalues on the path 

  ! This vector will store the Wannier-interpolated eigenvalues 
  if( allocated(evalpath)) deallocate( evalpath)
  allocate( evalpath( nstfv, nkpt))
  evalpath = 0.d0
!  call wannier_interpolate_eval( eval1( wf_fst:wf_lst,:), wf_kset%nkpt, wf_kset%vkl, &
!                              evalpath( wf_fst:wf_lst,:), nkpt, vkl, wf_fst, wf_lst)
  call wannier_interpolate_eval ( eval1( wf_fst:wf_lst,:), nkpt, vkl, evalpath( wf_fst:wf_lst,:), lmax=-1)
  !
  ! Print out the band on the path as a test  
  if (.true.) then
    call getunit(fid)
    open(fid,file='band_on_path.OUT',action='Write',status='Unknown')
    do ib = wf_fst, wf_lst
      do ik = 1, nkpt
         write(fid,*) ik, evalpath( ib, ik)    
      enddo
    enddo
    close(fid)
  endif  
  !
  if(allocated(selfeph0))deallocate(selfeph0)
  allocate(selfeph0(ibeph:nbeph,nkpt)) 
  selfeph0(:,:)=0.d0

  if(allocated(selfeph))deallocate(selfeph)
  allocate(selfeph(ibeph:nbeph,nomegeph,nkpt)) 
  selfeph(:,:,:)=0.d0

  if(allocated(speceph))deallocate(speceph)
  allocate(speceph(ibeph:nbeph,nomegeph,nkpt)) 
  speceph(:,:,:)=0.d0

  ! LOOP over all the k-points on the path 

  if (myrank==0) then
    call boxmsg(fgw,'=','EPH self-energy calculations starts ')
    call flushifc(fgw)
  end if

  !----------------------------------------------------
  ! generate the unshifted dense q mesh for BZ integral 
  !----------------------------------------------------
  !
  !print*, " Considering following grids for Wannier interpolation:"
  !print*, " Coarse grid: ", input%groundstate%ngridk
  !print*, " Dense  grid: ", ngridkq
  !print*, " Offset     : ", input%eph%vqloffeph
  !call generate_k_vectors(qsetd,bvec,ngridkqtot, &
  !                        input%eph%vqloffeph,.false.,.false.) ! generate a dense reducible k-point mesh

  do ik = 1 , nkpt 

   print*, " ik = ", ik , " / ", nkpt

  !--------------------------------------------
  ! for each k generate the dense k+q mesh for BZ integral 
  !--------------------------------------------
  !
  !print*, " Considering following grids for Wannier interpolation:"
  !print*, " Coarse grid: ", input%groundstate%ngridk
  !print*, " Dense  grid: ", ngridkqtot
  !print*, " Offset     : ", input%eph%vqloffeph
  ! 
  !call generate_k_vectors(kqsetd,bvec,ngridkqtot, &
  !                       input%eph%vqloffeph,.false.,.false.) ! generate a dense reducible k-point mesh
 

  call generate_kkqmt_vectors (kqsetd,gset,bvec,ngridkq, &
                          input%eph%vqloffeph, .false., vkl(:,ik),.false.)  

  !call getunit(fid)
  !call print_kkqmt_vectors (kqsetd,gset,fid)

  !print to file the k- and q- mesh used in the interpolation 
  if (.false.) then
    call getunit(fid)
    open(fid,file='EPH_QPOINTS-ik'//ik//'.OUT',action='Write',status='Unknown')
    call print_k_vectors(kqsetd%kqmtset,fid)
    call getunit(fid)
    open(fid,file='EPH_KPOINTS-ik'//ik//'.OUT',action='Write',status='Unknown')
    call print_k_vectors(wf_kset,fid)
  endif 
  !print*,  "Fitting from ", wf_kset%nkpt, " points to ", ngridkqtot

  !-----------------------------------------------------------------------
  ! Wannier interpolate the eigenvalues on the new mesh for BZ integration
  !-----------------------------------------------------------------------
  !
  ! This vector will store the Wannier-interpolated eigenvalues 
  if( allocated(eval2)) deallocate(eval2)
  allocate( eval2( nstfv, ngridkqtot))
  eval2 = 0.d0
  print*, 'ngridkqtot = ' , ngridkqtot
  print*, 'kqsetd%kqmtset%vkl', kqsetd%kqmtset%vkl

!  call wannier_interpolate_eval( eval1( wf_fst:wf_lst,:), wf_kset%nkpt, wf_kset%vkl, &
!                                 eval2( wf_fst:wf_lst,:), ngridkqtot, kqsetd%kqmtset%vkl, wf_fst, wf_lst)
  call wannier_interpolate_eval( eval1( wf_fst:wf_lst,:), ngridkqtot, kqsetd%kqmtset%vkl,eval2( wf_fst:wf_lst,:),lmax=-1)
                                 !eval2( wf_fst:wf_lst,:), kqsetd%kset%nkpt, kqsetd%kqmtset%vkl, wf_fst, wf_lst)
 
  ! find the Fermi energy on the new mesh
  if (ik == 1) then
    efnew = -10000.d0
    cbm   =  10000.d0
    do ib = wf_fst, wf_lst
      do ik1 = 1, ngridkqtot 
  
        efold  = efnew
        cbmold = cbm  
        evtmp  = eval2(ib,ik1)
        if (evtmp < efermi .and. evtmp > efold ) efnew = evtmp 
        if (evtmp > efermi .and. evtmp < cbmold) cbm   = evtmp
      enddo
    enddo
    print*, "Calculated Fermi energy : ", efnew * 27.21139 , ' eV' 
    print*, "Conduction band minimum : ", cbm   * 27.21139 , ' eV' 
    !efnew = cbm + 0.100 / 27.21139
  endif

  ! Print interpolated (and non-interpolated) eigenvalues on a aplane in the BZ
  if (.true. .and. ik .eq.1 ) then
    call getunit(fid)
    open(fid,file='band_plane-KS.OUT',action='Write',status='Unknown')
      ib=24 
      do ik1 = 1, wf_kset%nkpt
        if ( abs(wf_kset%vkl(3,ik1)) < 1.d-4) then 
          write(fid,*) wf_kset%vkl(1,ik1), wf_kset%vkl(2,ik1), eval1 ( ib, ik1)-efnew
        endif
      enddo
    close(fid)

    call getunit(fid)
    open(fid,file='band_plane-W.OUT',action='Write',status='Unknown')
      ib=24 
      do ik1 = 1, ngridkqtot
        if ( abs(kqsetd%kset%vkl(3,ik1)) < 1.d-4) then 
          write(fid,*) kqsetd%kset%vkl(1,ik1), kqsetd%kset%vkl(2,ik1), eval2 ( ib, ik1)-efnew
        endif
      enddo
    close(fid)
  endif  

  if (.true.) then
   if (ik.eq.1)then
    call getunit(fid)
    open(fid,file='EPH_KS-ev.OUT',action='Write',status='Unknown')
    do ik1 = 1 , wf_kset%nkpt 
      write(fid,100) ik1, wf_kset%vkl(1:3,ik1), (eval1(24, ik1) -efnew)* 27.21139

    enddo
    close(fid)
  
    call getunit(fid)
    open(fid,file='EPH_Wann-ev.OUT',action='Write',status='Unknown')
    do ik1 = 1 , ngridkqtot 
      write(fid,100) ik1, kqsetd%kqmtset%vkl(1:3,ik1), (eval2(24, ik1) -efnew)* 27.21139
    100 format(i4,3x,3f8.4,1x,3f10.6)
    enddo
    close(fid)
   endif 
  endif 


  !===========================================================================
  ! Main loop: BZ integration
  !===========================================================================    

#ifdef MPI
  call set_mpi_group(kqset%kqmtset%nkpt)
  call mpi_set_range(nproc_row, &
  &                  myrank_row, &
  &                  kqset%kqmtset%nkpt, 1, &
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
  !iqend = kqset%kqmtset%nkpt
  iomstart = 1
  !iomend = freq%nomeg
#endif

  l_speceph=.false.
  if (l_speceph)then 
    call calcspeceph(eval2, evalpath, ik)
  endif
  
  l_sigmaeph=.true.
  if (l_sigmaeph)then 
    call calcselfeph(eval2, evalpath, ik)
  endif 

  enddo 
 
  if(l_sigmaeph)then
    call getunit(fid)
    open(fid,file='ev-eph.OUT',action='Write',status='Unknown')
    do ib = ibeph, nbeph
      do ik = 1, nkpt
          write(fid,102) ik, ib, (evalpath(ib,ik)-efnew)* 27.21139, &
                    real(selfeph0(ib,ik))* 27.21139, aimag(selfeph0(ib,ik))* 27.21139  
        !write(fid,*) " " 
      enddo
    enddo
    close(fid)
  endif
  102 format(i4,2x,i4,3x,3f15.6)


  if(l_speceph)then
    call getunit(fid)
    open(fid,file='SIGMA.OUT',action='Write',status='Unknown')
    do ib = ibeph, nbeph
      do ik = 1, nkpt
        do iw = 1, nomegeph ! loop over frequency
          write(fid,101) ik, ib, (freq%freqs(iw))* 27.21139, (evalpath(ib,ik)-efnew)* 27.21139, &
                         real(selfeph(ib,iw,ik))* 27.21139, aimag(selfeph(ib,iw,ik))* 27.21139, & 
                         speceph(ib,iw,ik)
        enddo
        write(fid,*) " " 
      enddo
    enddo
    close(fid)
  endif
  101 format(i4,2x,i4,3x,7f15.6)

  ! compute the spectral function from the self energy 
  !call spec (selfeph,) 

  
  print*, "Stopping here " 
  call terminate
    
!#ifdef MPI
!    if ((nproc_row>1).and.(myrank_col==0)) then
!      write(*,*) "sum self-energy from different q-points"
!      call mpi_sum_array(0,selfex,nbandsgw,kset%nkpt,mycomm_row)
!      write(*,*) "sum selfex done"  
!      if (input%gw%taskname.ne.'g0w0_x') then
!        ! G0W0 and GW0 approximations
!        call mpi_sum_array(0,selfec,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
!        if (input%gw%selfenergy%secordw) then
!          call mpi_sum_array(0,selfecw2,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
!          call mpi_sum_array(0,selfecSR,nbandsgw,freq%nomeg,kset%nkpt,mycomm_row)
!        end if
!        if (input%gw%taskname=='cohsex') then
!          ! COHSEX
!          call mpi_sum_array(0,sigsx,nbandsgw,kset%nkpt,mycomm_row)
!          call mpi_sum_array(0,sigch,nbandsgw,kset%nkpt,mycomm_row)
!        end if ! cohsex
!        write(*,*) "sum selfec done"  
!      end if ! selfec
!    endif
!#endif
    
    if (myrank==0) then
    
      !===============================
      ! Write self-energies to files
      !===============================
!      call timesec(t0)
!#ifdef _HDF5_
!      do ik = 1, kset%nkpt
!        write(cik,'(I4.4)') ik
!        path = "/kpoints/"//trim(adjustl(cik))
!        if (.not.hdf5_exist_group(fgwh5,"/kpoints",cik)) &
!        &  call hdf5_create_group(fgwh5,"/kpoints",cik)
!        ! k-point
!        call hdf5_write(fgwh5,path,"vkl",kset%vkl(1,ik),(/3/))
!        ! exchange
!        call hdf5_write(fgwh5,path,"selfex",selfex(1,ik),(/nbandsgw/))
!        ! correlation
!        if (input%gw%taskname.ne.'g0w0_x') &
!        &  call hdf5_write(fgwh5,path,"selfec",selfec(1,1,ik),(/nbandsgw,freq%nomeg/))
!      end do
!#else
      !call write_selfenergy(ibgw,nbgw,kset%nkpt,freq%nomeg)
!#endif
!      call timesec(t1)
!      time_io = time_io+t1-t0
      
      ! KS band structure
      evalks(ibgw:nbgw,:) = evalsv(ibgw:nbgw,:) 
      !call bandanalysis('KS',ibgw,nbgw,evalks(ibgw:nbgw,:),efermi)
      call bandstructure_analysis('KS', &
      &  ibgw,nbgw,kset%nkpt,evalks(ibgw:nbgw,:),efermi)
      
      !=======================================
      ! Calculate the quasiparticle energies
      !=======================================
      !call calcevalqp
      
      !------------------------------------------------------
      ! Write quasi-particle energies to file
      !------------------------------------------------------
!      call timesec(t0)
      !call write_qp_energies('EVALQP.DAT')
!      call timesec(t1)
!      time_io = time_io+t1-t0
      
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
      
!      if (input%gw%selfenergy%secordw) then 
!      !---------------------------------
!      ! Second order screened exchange
!      !---------------------------------
!        call selfcGW_vs_selfcW2
!      end if
 
    end if ! myrank
    
    !--------------------------------------------------------
    ! Calculate quasiparticle energies in GW0 approximation 
    !--------------------------------------------------------
    if (input%gw%taskname=='gw0') then
      
      ! self-consistent cycle
      call calcscgw0
      
      ! print GW0 QP band structure
      if (myrank==0) then
!        call timesec(t0)
        !----------------------------------------
        ! Write quasi-particle energies to file
        !----------------------------------------
        call write_qp_energies('EVALQP-GW0.DAT')
        !call bandanalysis('GW0',ibgw,nbgw,evalqp(ibgw:nbgw,:),eferqp)
        call bandstructure_analysis('GW0', &
        &  ibgw,nbgw,kset%nkpt,evalqp(ibgw:nbgw,:),eferqp)
!        call timesec(t1)
!        time_io = time_io+t1-t0
      end if
    end if
    
    if (myrank==0) then
      !----------------------------------------
      ! Save QP energies into binary file
      !----------------------------------------
!      call timesec(t0)
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
    call delete_kkqmt_vectors(kqsetd)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine
!EOC
