
subroutine wfplot_new(ik,ist)
    use modmain
    use modinput
    use modplotlabels
    use mod_rgrid
    use mod_xsf_format
    implicit none
    ! input/output
    integer, intent(in) :: ik, ist
    ! local variables
    integer :: ip
    character(80) :: fname
    ! allocatable arrays      
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: wfmt(:,:,:,:,:)
    complex(8), allocatable :: wfir(:,:,:)
    complex(8), allocatable :: zdata(:)
    !
    type(rgrid) :: grid
    type(plot1d_type), pointer :: plotdef
    type(plotlabels), pointer :: labels

    ! initialise universal variables
    input%groundstate%lradstep = 1
    call init0
    call init1    
    ! read the density and potentials from file
    call readstate
    ! find the new linearisation energies
    call linengy
    ! generate the APW radial functions
    call genapwfr
    ! generate the local-orbital radial functions
    call genlofr

    if ((ik<1) .or. (ik>nkpt)) then
        write (*,*)
        write (*, '("Error(wfplot): k-point out of range : ", I8)') ik
        write (*,*)
        stop
    end if
    if ((ist<1) .or. (ist>nstsv)) then
        write (*,*)
        write (*, '("Error(wfplot): state out of range : ", I8)') ist
        write (*,*)
        stop
    end if

    ! allocate local arrays
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    allocate(evecfv(nmatmax,nstfv))
    allocate(evecsv(nstsv,nstsv))    

    ! get the eigenvectors and values from file
    call getevalsv(vkl(:,ik), evalsv)
    call getevecfv(vkl(:,ik), vgkl(:,:,:,ik), evecfv)
    call getevecsv(vkl(:,ik), evecsv)
    ! find the matching coefficients
    call match(ngk(1,ik), gkc(:,1,ik), tpgkc(:,:,1,ik), &
    &          sfacgk(:,:,1,ik), apwalm)

    ! calculate the wavefunctions for all states
    allocate(wfmt(lmmaxvr,nrmtmax,natmtot,nspinor,nstsv))
    allocate(wfir(ngrtot,nspinor,nstsv))
    !call genwfsv(.false., ngk(1, ik), igkig(:, 1, ik), evalsv, &
    !&            apwalm, evecfv, evecsv, wfmt, wfir)
    call genwfsv_new(ik, ist, ist, apwalm, evecfv, evecsv, wfmt, wfir)

    write(*,*) 'MT', sum(wfmt)
    write(*,*) 'IS', sum(wfir)

    !----------------
    ! 1D case
    !----------------
    If (associated(input%properties%wfplot%plot1d)) Then
      call gen_1d_rgrid(grid)
      !call print_rgrid(grid)
      ! Generate WF on the grid
      allocate(zdata(grid%npt))
      call calc_zdata_rgrid(grid, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
      write(fname,'("data1d-",i,"-",i,".dat")') ik, ist
      call str_strip(fname)
      open(77,file=trim(fname),status='Unknown',action='Write')
      do ip = 1, grid%npt
        write(77,'(3f16.6)') grid%vpd(ip), zdata(ip)
      end do
      close(77)
      call delete_rgrid(grid)
      write(*,*)
      write(*, '("Info(wfplot):")')
      write(*, '(" 1D Wavefunction written to data1d.dat")')
    end if

    !----------------
    ! 2D case
    !----------------
    if (associated(input%properties%wfplot%plot2d)) then
      call gen_2d_rgrid(grid)
      call print_rgrid(grid)
      ! Generate WF on the grid
      allocate(zdata(grid%npt))
      call calc_zdata_rgrid(grid, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
      call delete_rgrid(grid)
      write(fname,'("data2d-",i,"-",i,".xsf")') ik, ist
      call str_strip(fname)
      call write_2d_xsf(fname,grid%npt,dble(zdata))
      write (*,*)
      write (*, '("Info(wfplot):")')
      write (*, '(" 2D wavefunction  written to data2d.xsf")')
      write (*,*)
    end if

    !----------------
    ! 3D case
    !----------------
    if (associated(input%properties%wfplot%plot3d)) then
      call gen_3d_rgrid(grid)
      !call print_rgrid(grid)
      ! Generate WF on the grid
      allocate(zdata(grid%npt))
      call calc_zdata_rgrid(grid, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
      call delete_rgrid(grid)
      write(fname,'("data3d-",i,"-",i,".xsf")') ik, ist
      call str_strip(fname)
      call write_3d_xsf(fname,grid%npt,dble(zdata))
      Write(*,*)
      Write(*, '("Info(wfplot):")')
      Write(*, '(" 3D wavefunction written to data3d.xsf")')
      Write (*,*)
    end if

    write(*, '(" for k-point ", I6, " and state ", I6)') &
    &   ik, ist
    write(*,*)
    
    deallocate(apwalm, evecfv, evecsv, wfmt, wfir)
    deallocate(zdata)

    return
end subroutine

