
subroutine plot_exciton_wf

  use modinput
  use modmain
  use modxs
  use modmpi
  implicit none

  ! local
  integer :: lambda
  real(8) :: r0(3)
  
  integer :: nsta1, nsta2, nrnst1, nrnst3, hamsiz
  integer :: ipair
  character(22) :: fix
    
  real(8),    allocatable :: beval(:)
  complex(8), allocatable :: bevec(:,:)

  call init0
  call init1
  call init2
  call xssave0
  call readfermi

  if (.not.associated(input%xs%excitonPlot)) then
    write(*,*)
    write(*,*) 'Error(plot_exciton_wf): Element excitonplot is not specified!'
    write(*,*)
    stop
  end if

  do ipair = 1, size(input%xs%excitonplot%excitonarray)

    lambda = input%xs%excitonplot%excitonarray(ipair)%exciton%lambda
    write(*,'(" Exciton index: ",i4)') lambda

    r0(:)  = input%xs%excitonplot%excitonarray(ipair)%exciton%origin(:)
    write(*,'(" Origin coordinates: ",3f8.2)') r0

    fix = trim(input%xs%excitonplot%excitonarray(ipair)%exciton%fix)
    select case (trim(fix))
      case('hole')
        write(*,*) "Hole position is fixed"
        
      case('electron')
        write(*,*) "Electron position is fixed"

      case default
        write(*,*) "Error(plot_exciton_wf): Wrong value is specified!"
        write(*,*) "fix = ", trim(fix)
        stop

    end select

    if (associated(input%xs%excitonPlot%plot1d)) then
      write(*,*) '1d plot'
      call calc_1d_exciton_wf(lambda,r0,fix)

    else if (associated(input%xs%excitonPlot%plot2d)) then
      write(*,*) '2d plot'
      call calc_2d_exciton_wf(lambda,r0,fix)

    else if (associated(input%xs%excitonPlot%plot3d)) then
      write(*,*) '3d plot'
      call calc_3d_exciton_wf(lambda,r0,fix)

    else
      write(*,*)
      write(*,*) 'Error(plot_exciton_wf): Plot type is not specified!'
      write(*,*)
      stop
    end if

  end do ! over excitons

contains

  !--------------------------------------------------------------------------------
  subroutine calc_1d_exciton_wf(lambda,r0,fix)
    use mod_rgrid
    implicit none    
    integer      , intent(in) :: lambda
    real(8)      , intent(in) :: r0(3)
    character(*) , intent(in) :: fix
    ! local
    integer       :: iv, nv, ip, np
    character(80) :: fname
    type(rgrid)   :: r_h, r_e
    real(8),    allocatable :: vvl(:,:)
    complex(8), allocatable :: zwfeh(:)
 
    nv = size(input%xs%excitonPlot%plot1d%path%pointarray)
    write(*,*) nv
    if (nv < 1) then
      write (*,*)
      write (*,*) "Error(plot_exciton_wf::calc_1d_exciton_wf): Wrong plot specification!"
      write (*,*)
      stop
    end if
    np = input%xs%excitonPlot%plot1d%path%steps
    write(*,*) np
    If (np < nv) then
      write (*,*)
      write (*,*) "Error(wfplot_new): Wrong plot specification!"
      write (*,*)
      stop
    end if
    allocate(vvl(nv,3))
    do iv = 1, nv
      vvl(iv,:) =input%xs%excitonPlot%plot1d%path%pointarray(iv)%point%coord(:)
      write(*,*) vvl(iv,:)
    end do

    select case (trim(fix))
      case('hole')
        ! fixed position of hole, r_h
        r_h = gen_1d_rgrid(1, r0, 1)
        ! position of electrons, r_e
        r_e = gen_1d_rgrid(nv, vvl, np)
        ! e-h wavefunction
        allocate(zwfeh(r_e%npt))
        
      case('electron')
        ! position of hole, r_h
        r_h = gen_1d_rgrid(nv, vvl, np)
        ! fixed position of electron, r_e
        r_e = gen_1d_rgrid(1, r0, 1)
        ! e-h wavefunction
        allocate(zwfeh(r_h%npt))

    end select

    !call print_rgrid(r_h)
    !call print_rgrid(r_e)

    call read_exccoeff("EXCCOEFF.bin")
    call calc_eh_zwf(lambda,r_h,r_e,zwfeh)

    ! output
    if (rank==0) then
      select case (trim(fix))
        case('hole')
          write(fname,'("wf-eh-pair-1d-",i,"-fix-h.dat")') lambda
          call str_strip(fname)
          open(77,file=trim(fname),status='Unknown',action='Write')
          do ip = 1, r_e%npt
            !write(77,'(3f16.6)') r_e%vpd(ip), zwfeh(ip)
            write(77,'(2f16.6)') r_e%vpd(ip), abs(zwfeh(ip))
          end do
          close(77)

        case('electron')
          write(fname,'("wf-eh-pair-1d-",i,"-fix-e.dat")') lambda
          call str_strip(fname)
          open(77,file=trim(fname),status='Unknown',action='Write')
          do ip = 1, r_h%npt
            !write(77,'(3f16.6)') r_h%vpd(ip), zwfeh(ip)
            write(77,'(2f16.6)') r_h%vpd(ip), abs(zwfeh(ip))
          end do
          close(77)

        end select
    end if

    deallocate(zwfeh)
    call delete_rgrid(r_h)
    call delete_rgrid(r_e)
    deallocate(vvl)

    return
  end subroutine

  !--------------------------------------------------------------------------------
  subroutine calc_2d_exciton_wf(lambda,r0,fix)
    use mod_rgrid
    use mod_xsf_format
    implicit none    
    integer      , intent(in) :: lambda
    real(8)      , intent(in) :: r0(3)
    character(*) , intent(in) :: fix
    ! local
    integer :: ip, np
    integer :: igrid(3)
    real(8) :: boxl(3,3)
    complex(8), allocatable :: zwfeh(:)
    character(80) :: fname
    type(rgrid)   :: r_h, r_e

    igrid(:)  = input%xs%excitonPlot%plot2d%parallelogram%grid(1:3)
    boxl(1,:) = input%xs%excitonPlot%plot2d%parallelogram%origin%coord(1:3)
    boxl(2,:) = input%xs%excitonPlot%plot2d%parallelogram%pointarray(1)%point%coord(1:3)-boxl(1,:)
    boxl(3,:) = input%xs%excitonPlot%plot2d%parallelogram%pointarray(2)%point%coord(1:3)-boxl(1,:)

    select case (trim(fix))
      case('hole')
        ! fixed position of hole, r_h
        r_h = gen_1d_rgrid(1, r0, 1)
        ! position of electrons, r_e
        r_e = gen_2d_rgrid(igrid, boxl)
        ! e-h wavefunction
        allocate(zwfeh(r_e%npt))
        
      case('electron')
        ! position of hole, r_h
        r_h = gen_2d_rgrid(igrid, boxl)
        ! fixed position of electron, r_e
        r_e = gen_1d_rgrid(1, r0, 1)
        ! e-h wavefunction
        allocate(zwfeh(r_h%npt))

    end select

    !call print_rgrid(r_h)
    !call print_rgrid(r_e)

    call read_exccoeff("EXCCOEFF.bin")
    call calc_eh_zwf(lambda,r_h,r_e,zwfeh)

    ! output
    if (rank==0) then
      select case (trim(fix))
        case('hole')
          write(fname,'("wf-eh-pair-2d-",i,"-fix-h.xsf")') lambda
          call str_strip(fname)
          call write_structure_xsf(fname)
          !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))
          call write_2d_xsf(fname, 'module of electron-hole wavefunction', &
          &                 boxl, igrid, r_e%npt, abs(zwfeh))

        case('electron')
          write(fname,'("wf-eh-pair-2d-",i,"-fix-e.xsf")') lambda
          call str_strip(fname)
          call write_structure_xsf(fname)
          !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))
          call write_2d_xsf(fname, 'module of electron-hole wavefunction', &
          &                 boxl, igrid, r_h%npt, abs(zwfeh))

      end select
    end if

    deallocate(zwfeh)
    call delete_rgrid(r_h)
    call delete_rgrid(r_e)

    return
  end subroutine

  !--------------------------------------------------------------------------------
  subroutine calc_3d_exciton_wf(lambda,r0,fix)
    use mod_rgrid
    use mod_xsf_format
    implicit none    
    integer      , intent(in) :: lambda
    real(8)      , intent(in) :: r0(3)
    character(*) , intent(in) :: fix
    ! local
    integer :: ip, np
    integer :: igrid(3)
    real(8) :: boxl(4,3)
    complex(8), allocatable :: zwfeh(:)
    character(80) :: fname
    type(rgrid)   :: r_h, r_e

    igrid(:)  = input%xs%excitonPlot%plot3d%box%grid(1:3)
    boxl(1,:) = input%xs%excitonPlot%plot3d%box%origin%coord(1:3)
    boxl(2,:) = input%xs%excitonPlot%plot3d%box%pointarray(1)%point%coord(1:3)-boxl(1,:)
    boxl(3,:) = input%xs%excitonPlot%plot3d%box%pointarray(2)%point%coord(1:3)-boxl(1,:)
    boxl(4,:) = input%xs%excitonPlot%plot3d%box%pointarray(3)%point%coord(1:3)-boxl(1,:)

    select case (trim(fix))
      case('hole')
        ! fixed position of hole, r_h
        r_h = gen_1d_rgrid(1, r0, 1)
        ! position of electrons, r_e
        r_e = gen_3d_rgrid(igrid, boxl)
        ! e-h wavefunction
        allocate(zwfeh(r_e%npt))
        
      case('electron')
        ! position of hole, r_h
        r_h = gen_3d_rgrid(igrid, boxl)
        ! fixed position of electron, r_e
        r_e = gen_1d_rgrid(1, r0, 1)
        ! e-h wavefunction
        allocate(zwfeh(r_h%npt))

    end select

    !call print_rgrid(r_h)
    !call print_rgrid(r_e)

    call read_exccoeff("EXCCOEFF.bin")
    call calc_eh_zwf(lambda, r_h, r_e, zwfeh)

    ! output
    if (rank==0) then
      select case (trim(fix))
        case('hole')
          write(fname,'("wf-eh-pair-3d-",i,"-fix-h.xsf")') lambda
          call str_strip(fname)
          call write_structure_xsf(fname)
          !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))
          call write_3d_xsf(fname, 'module of the electron-hole wavefunction', &
          &                 boxl, igrid, r_e%npt, abs(zwfeh))

        case('electron')
          write(fname,'("wf-eh-pair-3d-",i,"-fix-e.xsf")') lambda
          call str_strip(fname)
          call write_structure_xsf(fname)
          !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))
          call write_3d_xsf(fname, 'module of the electron-hole wavefunction', &
          &                 boxl, igrid, r_h%npt, abs(zwfeh))

      end select
    end if

    deallocate(zwfeh)
    call delete_rgrid(r_h)
    call delete_rgrid(r_e)

    return
  end subroutine


  subroutine calc_eh_zwf(lambda,r_h,r_e,zwfeh)
    use mod_rgrid
    implicit none
    integer,     intent(in) :: lambda
    type(rgrid), intent(in) :: r_h, r_e
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

    logical :: active
    integer :: nka
    integer, allocatable :: ikMap(:)
    integer, allocatable :: nva(:), ivMap(:,:)
    integer, allocatable :: nca(:,:), icMap(:,:,:)
    integer :: ika, iva, ica

    !write(*,*) istl3, nsta1, nsta2, nrnst1, nrnst3

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
    
    nka = 0
    allocate(ikMap(nkptnr))
    ikMap(:) = 0

    allocate(nva(nkptnr))
    nva(:) = 0
    allocate(ivMap(nrnst1,nkptnr))
    ivMap(:,:) = 0

    allocate(nca(nrnst1,nkptnr))
    nca(:,:) = 0
    allocate(icMap(nrnst3,nrnst1,nkptnr))
    icMap(:,:,:) = 0

    do ik = 1, nkptnr
      do iv = 1, nrnst1
        do ic = 1, nrnst3
          ! combined ivck index (bse.f90 function hamidx)
          ivck = ic + nrnst3*(iv-1) + nrnst1*nrnst3*(ik-1)
          if (abs(bevec(ivck,lambda)) > input%xs%excitonplot%epstol) then
            active = .true.
            nca(iv,ik) = nca(iv,ik)+1
            icMap(nca(iv,ik),iv,ik) = ic
          else
            active = .false.
          end if
        end do
        if (active) then
          nva(ik) = nva(ik)+1
          ivMap(nva(ik),ik) = iv
        end if
      end do
      if (active) then
        nka = nka+1
        ikMap(nka) = ik
      end if
    end do

    write(*,*) 'nka=', nka, nkptnr
    write(*,*) 'nva=', nva(1), nrnst1
    write(*,*) 'nca=', nca(:,1), nrnst3
    ! stop

    !---------------------
    ! Sum over ik, iv, ic
    !---------------------

#ifdef MPI
    do ika = firstofset(rank,nka), lastofset(rank,nka)
#else
    do ika = 1, nka
#endif
      ik = ikMap(ika)
      !write(*,'(a,i,3f12.4)') 'vkl: ', ik, vkl(:,ik)

      ! read eigenvectors
      call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfvt)
      call getevecsv(vkl(:,ik), evecsvt)
      ! find the matching coefficients
      call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), &
      &          sfacgk(:,:,1,ik), apwalm)
      ! call genwfsv_new(ik, nsta1, istl3+nsta2-1, apwalm, evecfvt, evecsvt, wfmt, wfir)
      call genwfsv_new(ik, 1, nstsv, apwalm, evecfvt, evecsvt, wfmt, wfir)

      do iva = 1, nva(ik)
        iv = ivMap(iva,ik)
        ist = iv+nsta1-1
        ! hole WF 
        call calc_zdata_rgrid(r_h, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zwfrh)
        ! if (ik==1) then
        !   write(fname,'("wf-hole-1d-",i,"-",i,".dat")') ik, ist
        !   call str_strip(fname)
        !   open(77,file=trim(fname),status='Unknown',action='Write')
        !   do ip = 1, r_h%npt
        !     write(77,'(3f16.6)') r_h%vpd(ip), zwfrh(ip)
        !   end do
        !   close(77)
        ! end if

        do ica = 1, nca(iv,ik)
          ic = icMap(ica,iv,ik)
          jst = ic+istl3-1+nsta2-1
          ! electron WF 
          call calc_zdata_rgrid(r_e, ik, wfmt(:,:,:,1,jst), wfir(:,1,jst), zwfre)
          ! if (ik==2) then
          !   write(fname,'("wf-elec-1d-",i,"-",i,".dat")') ik, jst
          !   call str_strip(fname)
          !   open(77,file=trim(fname),status='Unknown',action='Write')
          !   do ip = 1, r_e%npt
          !     write(77,'(3f16.6)') r_e%vpd(ip), zwfre(ip)
          !   end do
          !   close(77)
          ! end if

          ! combined ivck index (bse.f90 function hamidx)
          ivck = ic + nrnst3*(iv-1) + nrnst1*nrnst3*(ik-1)
          if (r_h%npt < r_e%npt) then
            ! fixed hole
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
    prob = omega*sum(abs(zwfeh)**2)/dble(npt)
    write(*,'(a,e15.4/)') 'prob=', prob

    deallocate(zwfrh,zwfre)
    deallocate(apwalm,evecfvt,evecsvt)
    deallocate(wfmt,wfir)
    deallocate(beval,bevec)

    deallocate(nva,nca)
    deallocate(ikMap,ivMap,icMap)

    return
  end subroutine

  !--------------------------------------------------------------------------------
  subroutine read_exccoeff(fname)
    implicit none
    character(*) :: fname
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
      &        nkptnr, istl3, nsta1, nsta2, nrnst1, nrnst3, hamsiz
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

    write(*,*) MinNumberExcitons, MaxNumberExcitons
    write(*,*) nkptnr, istl3, nsta1, nsta2
    write(*,*) nrnst1, nrnst3, hamsiz

    return
  end subroutine

end subroutine ! plot_exciton_wf
