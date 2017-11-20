
module mod_exciton_wf
  use modinput
  use modmain
  use modxs
  use modxas
  use modmpi
  use m_putgetexcitons
  implicit none

  integer ::  nrnst1, nrnst3, hamsiz
  real(8),    allocatable :: beval(:)
  complex(8), allocatable :: bevec(:,:)

  ! index maps when sorting KS states by their contribution to excitons
  integer :: nka
  integer, allocatable :: kMap(:)
  integer, allocatable :: nva(:), vMap(:,:)
  integer, allocatable :: nca(:,:), cMap(:,:,:)

  ! Note:
  ! ipair is used further to generate a unique output file name
  ! coding all e-h pair descriptors
  integer, private :: ipair

  public  :: read_exccoeff, plot_excitonWavefunction

  private :: calc_eh_zwf, write_zwfeh, calc_eh_zwfcore

contains

  !--------------------------------------------------------------------------------
  subroutine plot_excitonWavefunction()
    use mod_rgrid
    use mod_xsf_format
    use mod_cube_format
    implicit none    
    ! local
    integer :: i, ip, np, lambda
    integer :: igrid(3)
    integer :: ndim_h, ndim_e
    real(8) :: boxl(4,3)
    real(8),    allocatable :: vvl(:,:)
    complex(8), allocatable :: zwfeh(:)
    character(256) :: fname, label
    type(rgrid)   :: r_h, r_e, r0
    character(24) :: fix

    if (.not.associated(input%xs%excitonPlot)) then
      write(*,*)
      write(*,*) 'Error(mod_exciton_wf::plot_excitonWavefunction): Element excitonPlot is not specified!'
      write(*,*)
      stop
    end if

    if (.not.associated(input%xs%excitonPlot%hole)) then
      write(*,*)
      write(*,*) 'Error(mod_exciton_wf::plot_excitonWavefunction): Real space grid for electrons is not specified!'
      write(*,*)
      stop
    end if

    if (.not.associated(input%xs%excitonPlot%electron)) then
      write(*,*)
      write(*,*) 'Error(mod_exciton_wf::plot_excitonWavefunction): Real space grid for electrons is not specified!'
      write(*,*)
      stop
    end if
    !----------------------------
    ! sanity checks for core excitations
    !----------------------------    
    if (input%xs%BSE%xas) then
      do ipair = 1, size(input%xs%excitonplot%excitonarray)
        fix = trim(input%xs%excitonplot%excitonarray(ipair)%exciton%fix)
        if (fix .eq. 'electron') then
          write(*,*)
          write(*,*) 'Error(mod_exciton_wf::plot_excitonWavefunction): Core hole visualization not implemented yet!'
          write(*,*)
          stop
        end if
      end do
    end if
    !----------------------------
    ! global initialization
    !----------------------------
    call init0
    call init1
    call init2
    call xssave0
    call readfermi

    !----------------------------
    ! HOLE
    !----------------------------
    if (associated(input%xs%excitonPlot%hole%plot1d)) then
      ndim_h = 1
      r_h = gen_1d_rgrid(input%xs%excitonPlot%hole%plot1d)
    
    else if (associated(input%xs%excitonPlot%hole%plot2d)) then
      ndim_h = 2
      r_h = gen_2d_rgrid(input%xs%excitonPlot%hole%plot2d, 0)
    
    else if (associated(input%xs%excitonPlot%hole%plot3d)) then
      ndim_h = 3
      r_h = gen_3d_rgrid(input%xs%excitonPlot%hole%plot3d, 0)
         
    else
      write(*,*)
      write(*,*) 'Error(mod_exciton_wf::plot_excitonWavefunction): Plot type is not specified!'
      write(*,*)
      stop
    end if
    ! call print_rgrid(r_h)

    !----------------------------
    ! ELECTRON
    !----------------------------
    if (associated(input%xs%excitonPlot%electron%plot1d)) then
      ndim_e = 1
      r_e = gen_1d_rgrid(input%xs%excitonPlot%electron%plot1d)
    
    else if (associated(input%xs%excitonPlot%electron%plot2d)) then
      ndim_e = 2
      r_e = gen_2d_rgrid(input%xs%excitonPlot%electron%plot2d, 0)
    
    else if (associated(input%xs%excitonPlot%electron%plot3d)) then
      ndim_e = 3
      r_e = gen_3d_rgrid(input%xs%excitonPlot%electron%plot3d, 0)
         
    else
      write(*,*)
      write(*,*) 'Error(mod_exciton_wf::plot_excitonWavefunction): Plot type is not specified!'
      write(*,*)
      stop
    end if
    ! call print_rgrid(r_e)

    !----------------------------------------
    ! read exciton wavefunction coefficients
    !----------------------------------------
    ! Get data form EXCCOEFF*.OUT
    call get_excitons()
    ! Check compatibility (needs TDA and fixed number of transition over all k)
    if(fcoup_ .eqv. .true. .or. fesel_ .eqv. .true.) then 
      write(*,*) "Error(plot_excitonWavefunction): Data retrieved with get_excitions invalid"
      call terminate
    end if
    ! Copy to corresponding module variables
    nkptnr=nk_bse_
    ! lowest v state (absolute index)
    sta1=koulims_(3,1)
    ! number of v states
    nrnst1=koulims_(4,1)-koulims_(3,1)+1
    ! reference c state 
    istl3=iuref_
    ! lowest c state (relative index)
    sta2 = koulims_(1,1)-iuref_+1
    ! number of c states
    nrnst3=koulims_(2,1)-koulims_(1,1)+1
    if(allocated(beval)) deallocate(beval)
    allocate(beval(hamsize_))
    if(allocated(bevec)) deallocate(bevec)
    allocate(bevec(hamsize_,iex1_:iex2_))
    beval=evals_
    bevec=rvec_
    call clear_excitons
   
    !----------------------------------------
    ! Loop over excitons
    !----------------------------------------
    do ipair = 1, size(input%xs%excitonplot%excitonarray)

      lambda = input%xs%excitonPlot%excitonarray(ipair)%exciton%lambda
      write(*,'(" Exciton index: ",i4)') lambda

      fix = trim(input%xs%excitonplot%excitonarray(ipair)%exciton%fix)
      select case (trim(fix))

        case('hole')
          write(*,*) "hole positions are fixed"
          do ip = 1, r_h%npt
            ! fixed position of electron, r_h
            r0 = gen_1d(1, r_h%vpl(:,ip), 1)
            allocate(zwfeh(r_e%npt))
            if (input%xs%BSE%xas) then
              call calc_eh_zwfcore(lambda, r0, r_e, zwfeh)
            else
              call calc_eh_zwf(lambda, r0, r_e, zwfeh)
            end if
            ! output
            if (rank==0) then
              write(fname,'("excitonWavefunction-", I6.6, "-", I6.6, ".xsf")') lambda, ip
              write(label,'("hole position: ", 3F12.6)') r0%vpl(:,1)
              call write_structure_xsf(fname)
              select case (ndim_e)
                case(1)
                  call write_zwfeh(fname, r_e, abs(zwfeh)**2)
                case(2)
                  call write_2d_xsf(fname, label, r_e%boxl, r_e%ngrid, r_e%npt, abs(zwfeh)**2)
                case(3)
                  call write_3d_xsf(fname, label, r_e%boxl, r_e%ngrid, r_e%npt, abs(zwfeh)**2)
                  write(fname,'("excitonWavefunction-", I6.6, "-", I6.6, ".cube")') lambda, ip
                  call write_3d_cube(fname, label, r_e%boxl, r_e%ngrid, r_e%npt, abs(zwfeh)**2)
              end select
            end if
            call delete_rgrid(r0)
            deallocate(zwfeh)
          end do ! ip
          
        case('electron')
          write(*,*) "Electron positions are fixed"
          do ip = 1, r_e%npt
            ! fixed position of electron, r_h
            r0 = gen_1d(1, r_e%vpl(:,ip), 1)
            allocate(zwfeh(r_h%npt))
            call calc_eh_zwf(lambda, r_h, r0, zwfeh)
            ! output
            if (rank==0) then
              write(fname,'("excitonWavefunction-", I6.6, "-", I6.6, ".xsf")') lambda, ip
              write(label,'("electron position: ", 3F12.6)') r0%vpl(:,1)
              call write_structure_xsf(fname)
              select case (ndim_h)
                case(1)
                  call write_zwfeh(fname, r_h, abs(zwfeh)**2)
                case(2)
                  call write_2d_xsf(fname, label, r_h%boxl, r_h%ngrid, r_h%npt, abs(zwfeh)**2)
                case(3)
                  call write_3d_xsf(fname, label, r_h%boxl, r_h%ngrid, r_h%npt, abs(zwfeh)**2)
                  write(fname,'("excitonWavefunction-", I6.6, "-", I6.6, ".cube")') lambda, ip
                  call write_3d_cube(fname, label, r_h%boxl, r_h%ngrid, r_h%npt, abs(zwfeh)**2)
              end select
            end if
            call delete_rgrid(r0)
            deallocate(zwfeh)
          end do ! ip

        case default
          write(*,*) "Error(mod_exciton_wf::plot_excitonWavefunction): Wrong value is specified!"
          write(*,*) "fix = ", trim(fix)
          stop
      end select

    end do ! ipair

    deallocate(beval,bevec)
    call delete_rgrid(r_h)
    call delete_rgrid(r_e)

    return
  end subroutine

!------------------------------------------------------------------------------
  subroutine calc_eh_zwf(lambda,r_h,r_e,zwfeh)
    use mod_rgrid
    implicit none
    integer,     intent(in) :: lambda
    type(rgrid), intent(in) :: r_h, r_e
    complex(8),  intent(out):: zwfeh(:)
    ! local
    integer :: ivck, ik, iv, ic, ist, jst, ip, npt
    real(8) :: prob
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)
    complex(8), allocatable :: wfmt(:,:,:,:,:)
    complex(8), allocatable :: wfir(:,:,:)
    complex(8), allocatable :: zwfrh(:), zwfre(:)

    integer :: ika, iva, ica

    !--------------------------------------------------------------------------
    ! \Psi^{l}(r_h,r_e) = \sum_{vck} A^{l}_{vck} \psi^{*}_{vk}(r_h) \psi_{ck}(r_e)
    !--------------------------------------------------------------------------
  
    filext = '_QMT001.OUT'

    ! allocate local arrays
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfvt(nmatmax,nstfv))
    allocate(evecsvt(nstsv,nstsv))

    ! calculate the wavefunctions for all states
    allocate(wfmt(lmmaxapw,nrmtmax,natmtot,nspinor,nstsv))
    allocate(wfir(ngrtot,nspinor,nstsv))

    allocate(zwfrh(r_h%npt))
    allocate(zwfre(r_e%npt))

    npt = size(zwfeh)
    zwfeh(:) = 0.d0

    !----------------------------------------------------------
    ! Determine which transitions do contribute to the exciton
    !----------------------------------------------------------
    call sort_transitions(lambda, input%xs%excitonplot%epstol)

    !---------------------
    ! Sum over ik, iv, ic
    !---------------------

#ifdef MPI
    do ika = firstofset(rank,nka), lastofset(rank,nka)
#else
    do ika = 1, nka
#endif
      ik = kMap(ika)
      ! write(*,'(a,i,3f12.4)') 'vkl: ', ik, vkl(:,ik)

      ! read eigenvectors
      call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfvt)
      call getevecsv(vkl(:,ik), evecsvt)
      ! find the matching coefficients
      call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), &
      &          sfacgk(:,:,1,ik), apwalm)
      call genwfsv_new(ik, 1, nstsv, apwalm, evecfvt, evecsvt, wfmt, wfir)

      do iva = 1, nva(ik)
        iv = vMap(iva,ik)
        ist = iv+sta1-1
        ! electron WF 
        call calc_zdata_rgrid(r_h, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zwfrh)

        do ica = 1, nca(iv,ik)
          ic = cMap(ica,iv,ik)
          jst = ic+istl3-1
          ! electron WF 
          call calc_zdata_rgrid(r_e, ik, wfmt(:,:,:,1,jst), wfir(:,1,jst), zwfre)

          ! combined ivck index (bse.f90 function hamidx)
          ivck = ic + nrnst3*(iv-1) + nrnst1*nrnst3*(ik-1)
          if (r_h%npt == 1) then
            ! fixed electron
            do ip = 1, npt
              zwfeh(ip) = zwfeh(ip) + bevec(ivck,lambda)*conjg(zwfrh(1))*zwfre(ip)
            end do
          else
            ! fixed electron
            do ip = 1, npt
              zwfeh(ip) = zwfeh(ip) + bevec(ivck,lambda)*conjg(zwfrh(ip))*zwfre(1)
            end do
          end if
        end do ! ic
      end do ! iv
    end do ! ik

#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, zwfeh, npt, &
    &                  MPI_DOUBLE_COMPLEX,  MPI_SUM, &
    &                  MPI_COMM_WORLD, ierr)
#endif    

    ! total probability
    ! prob = omega*sum(abs(zwfeh)**2)/dble(npt)
    ! write(*,'(a,e15.4/)') 'prob=', prob

    deallocate(zwfrh,zwfre)
    deallocate(apwalm,evecfvt,evecsvt)
    deallocate(wfmt,wfir)

    call clean_transitions

    return
  end subroutine
 !--------------------------------------------------------------------------------
  subroutine calc_eh_zwfcore(lambda,r_h,r_e,zwfeh)
    use mod_rgrid
    use modxas
    implicit none
    integer,     intent(in) :: lambda
    type(rgrid), intent(in) :: r_e, r_h
    complex(8),  intent(out):: zwfeh(:)
    ! local
    integer :: ivck, ik, iv, ic, ist, jst, ip, npt
    real(8) :: prob
    character(80) :: fname
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfvt(:,:)
    complex(8), allocatable :: evecsvt(:,:)
    complex(8), allocatable :: wfmt(:,:,:,:,:)
    complex(8), allocatable :: wfir(:,:,:)
    complex(8), allocatable :: zwfrh(:), zwfre(:)

    integer :: ika, iva, ica

    !--------------------------------------------------------------------------
    ! \Psi^{l}(r_h,r_e) = \sum_{vck} A^{l}_{vck} \psi^{*}_{vk}(r_h) \psi_{ck}(r_e)
    !--------------------------------------------------------------------------
  
    filext = '_QMT001.OUT'

    ! allocate local arrays
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfvt(nmatmax,nstfv))
    allocate(evecsvt(nstsv,nstsv))

    ! calculate the wavefunctions for all states
    allocate(wfmt(lmmaxapw,nrmtmax,natmtot,nspinor,nstsv))
    allocate(wfir(ngrtot,nspinor,nstsv))

    allocate(zwfre(r_e%npt))
    allocate(zwfrh(r_h%npt))

    npt = size(zwfeh)
    zwfeh(:) = 0.d0

    !----------------------------------------------------------
    ! Determine which transitions do contribute to the exciton
    !----------------------------------------------------------
    call sort_transitions(lambda, input%xs%excitonplot%epstol)

    !---------------------
    ! Sum over ik, iv, ic
    !---------------------

#ifdef MPI
    do ika = firstofset(rank,nka), lastofset(rank,nka)
#else
    do ika = 1, nka
#endif
      ik = kMap(ika)
      !write(*,'(a,i,3f12.4)') 'vkl: ', ik, vkl(:,ik)

      ! read eigenvectors
      call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfvt)
      call getevecsv(vkl(:,ik), evecsvt)
      ! find the matching coefficients
      call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), &
      &          sfacgk(:,:,1,ik), apwalm)
      call genwfsv_new(ik, 1, nstsv, apwalm, evecfvt, evecsvt, wfmt, wfir)

      do iva = 1, nva(ik)
        iv = vMap(iva,ik)
        ist = iv+sta1-1
        ! hole WF 
        call calc_zdata_rgrid_core(r_h, ist, ucore(:,ist), zwfrh)
        do ica = 1, nca(iv,ik)
          ic = cMap(ica,iv,ik)
          jst = ic+istl3-1+sta2-1
          ! electron WF 
          call calc_zdata_rgrid(r_e, ik, wfmt(:,:,:,1,jst), wfir(:,1,jst), zwfre)

          ! combined ivck index (bse.f90 function hamidx)
          ivck = ic + nrnst3*(iv-1) + nrnst1*nrnst3*(ik-1)
          ! fixed hole
          do ip = 1, npt
              zwfeh(ip) = zwfeh(ip) + bevec(ivck,lambda)*zwfre(ip)*conjg(zwfrh(1))
          end do
        end do ! ic
      end do ! iv
    end do ! ik

#ifdef MPI
    call MPI_AllReduce(MPI_IN_PLACE, zwfeh, npt, &
    &                  MPI_DOUBLE_COMPLEX,  MPI_SUM, &
    &                  MPI_COMM_WORLD, ierr)
#endif    

    ! total probability
    !prob = omega*sum(zwfeh)/dble(npt)
    !write(*,'(a,e15.4/)') 'prob=', prob

    deallocate(zwfrh,zwfre)
    deallocate(apwalm,evecfvt,evecsvt)
    deallocate(wfmt,wfir)

    call clean_transitions

    return
  end subroutine

  !--------------------------------------------------------------------------------
  subroutine read_exccoeff(fname)
    implicit none
    character*(*) :: fname
    logical :: exist
    integer :: iostat
    integer :: MinNumberExcitons, MaxNumberExcitons
    integer :: hamsiz, ievec

    inquire(file=trim(fname), exist=exist)
    if (.not.exist) then
      write(*,*)
      write(*,'("Error(plot_exciton_wf::read_exccoeff): File ",a," does not exist!")') trim(fname)
      write(*,*)
      stop
    else
      open(50,File=trim(fname), & 
      &    Action='READ',Form='UNFORMATTED', IOstat=iostat)
      if (iostat/=0) then
        write(*,*) iostat
        write(*,'("Error(plot_exciton_wf::read_exccoeff): error reading ",a)') trim(fname)
        write(*,*)
        stop
      end if
      read(50) MinNumberExcitons, MaxNumberExcitons, &
      &        nkptnr, istl3, sta1, sta2, nrnst1, nrnst3, hamsiz
      if (allocated(beval)) deallocate(beval)
      allocate(beval(hamsiz))
      beval = 0.d0
      if (allocated(bevec)) deallocate(bevec)
      allocate(bevec(hamsiz,MinNumberExcitons:MaxNumberExcitons))
      bevec = 0.d0
      do ievec = MinNumberExcitons, MaxNumberExcitons
        read(50) beval(ievec), bevec(1:hamsiz,ievec)
      end do
      close(50)
    end if

    ! write(*,*) "Info(mod_exciton_wf::read_exccoeff):"
    ! write(*,*) MinNumberExcitons, MaxNumberExcitons
    ! write(*,*) nkptnr, istl3, sta1, sta2
    ! write(*,*) nrnst1, nrnst3, hamsiz
    ! write(*,*)

    return
  end subroutine

  !--------------------------------------------------------------------------------
  subroutine write_zwfeh(fname,r_grid,zdata)
    use mod_rgrid
    implicit none
    character*(*), intent(in) :: fname
    type(rgrid),   intent(in) :: r_grid
    real(8),       intent(in) :: zdata(:)
    ! local
    integer :: ip
    if (rank==0) then
      open(80,file=trim(fname),status='Unknown',action='Write')
      do ip = 1, r_grid%npt
        write(80,'(2f16.6)') r_grid%vpd(ip), zdata(ip)
      end do
      close(80)
    end if
    return
  end subroutine

  !----------------------------------------------------------------------------
  subroutine sort_transitions(lambda,epstol)
    implicit none
    integer, intent(in) :: lambda
    real(8), intent(in) :: epstol
    integer :: ik, iv, ic, ivck
    logical :: active

    nka = 0
    allocate(kMap(nkptnr))
    kMap(:) = 0

    allocate(nva(nkptnr))
    nva(:) = 0
    allocate(vMap(nrnst1,nkptnr))
    vMap(:,:) = 0

    allocate(nca(nrnst1,nkptnr))
    nca(:,:) = 0
    allocate(cMap(nrnst3,nrnst1,nkptnr))
    cMap(:,:,:) = 0
    do ik = 1, nkptnr
      active = .false.
      do iv = 1, nrnst1
        active = .false.
        do ic = 1, nrnst3
          ! combined ivck index (bse.f90 function hamidx)
          ivck = ic + nrnst3*(iv-1) + nrnst1*nrnst3*(ik-1)
          if (abs(bevec(ivck,lambda)) > epstol) then
            active = .true.
            nca(iv,ik) = nca(iv,ik)+1
            cMap(nca(iv,ik),iv,ik) = ic
          end if
        end do ! ic
        if (active) then
          nva(ik) = nva(ik)+1
          vMap(nva(ik),ik) = iv
        end if
      end do ! iv
      if (active) then
        nka = nka+1
        kMap(nka) = ik
      end if
    end do

    if (nka==0) then
      write(*,*)
      write(*,*) 'Error(mod_exciton_wf::sort_transitions) None of the transitions'
      write(*,*) '  satisfies to the specified tolerance factor!', input%xs%excitonplot%epstol
      write(*,*) '  Decrease the value of input%xs%excitonplot%epstol'
      write(*,*)
    end if

    ! write(*,*) 'nka=', nka, nkptnr
    ! write(*,*) 'nva=', nva(1), nrnst1
    ! write(*,*) 'nca=', nca(:,1), nrnst3
    ! stop

    return
  end subroutine

  subroutine clean_transitions()
    implicit none
    deallocate(nva,nca)
    deallocate(kMap,vMap,cMap)
  end subroutine

end module
