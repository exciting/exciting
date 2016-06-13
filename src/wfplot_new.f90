
subroutine wfplot_new(ik,ist)
    use modmain
    use modinput
    use modplotlabels
    use mod_rgrid
    use mod_xsf_format
    use mod_cube_format
    use modmpi, only : rank
    implicit none
    ! input/output
    integer, intent(in) :: ik, ist
    ! local variables
    integer :: ip, np, iv, nv
    character(80) :: fname
    integer :: igrid(3)
    real(8) :: boxl(4,3)
    ! allocatable arrays
    real(8),    allocatable :: vvl(:,:)
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
      if (rank==0) then
        write (*,*)
        write (*, '("Error(wfplot): k-point out of range : ", I8)') ik
        write (*,*)
      end if
      stop
    end if
    if ((ist<1) .or. (ist>nstsv)) then
      if (rank==0) then
        write (*,*)
        write (*, '("Error(wfplot): state out of range : ", I8)') ist
        write (*,*)
      end if
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
    allocate(wfmt(lmmaxapw,nrmtmax,natmtot,nspinor,nstsv))
    allocate(wfir(ngrtot,nspinor,nstsv))
    !call genwfsv(.false., ngk(1, ik), igkig(:, 1, ik), evalsv, &
    !&            apwalm, evecfv, evecsv, wfmt, wfir)
    call genwfsv_new(ik, ist, ist, apwalm, evecfv, evecsv, wfmt, wfir)

    !----------------
    ! 1D case
    !----------------
    If (associated(input%properties%wfplot%plot1d)) then
      nv = size(input%properties%wfplot%plot1d%path%pointarray)
      !write(*,*) nv
      if (nv < 1) then
        if (rank==0) then
          write (*,*)
          write (*,*) "Error(wfplot_new): Wrong plot specification!"
          write (*,*)
        end if
        stop
      end if
      np = input%properties%wfplot%plot1d%path%steps
      !write(*,*) np
      If (np < nv) then
        if (rank==0) then
          write (*,*)
          write (*,*) "Error(wfplot_new): Wrong plot specification!"
          write (*,*)
        end if
        stop
      end if
      allocate(vvl(nv,3))
      do iv = 1, nv
        vvl(iv,:) = input%properties%wfplot%plot1d%path%pointarray(iv)%point%coord
        !write(*,*) vvl(iv,:)
      end do

      ! rgrid constructor
      grid = gen_1d_rgrid(nv, vvl, np)
      !call print_rgrid(grid)
      deallocate(vvl)

      ! Generate WF on the grid
      allocate(zdata(grid%npt))
      call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
      call delete_rgrid(grid)

      ! Output
      if (rank==0) then
        write(fname,'("wf1d-",i,"-",i,".dat")') ik, ist
        call str_strip(fname)
        open(77,file=trim(fname),status='Unknown',action='Write')
        do ip = 1, grid%npt
          write(77,'(i,3f16.6)') ip, zdata(ip), abs(zdata(ip))**2
          !write(77,'(3f16.6)') grid%vpd(ip), zdata(ip)
          !write(77,'(2f16.6)') grid%vpd(ip), wkpt(ik)*nkptnr*abs(zdata(ip))**2
        end do
        close(77)
        write(*,*)
        write(*,'("Info(wfplot):")')
        write(*,'(" 1D Wavefunction written to wf1d-ik-ist.dat")')
        write(*,*)
        write(*,'(" for k-point ", I6, " and state ", I6)') ik, ist
        write(*,*)
      end if

      deallocate(zdata)
      
    end if

    !----------------
    ! 2D case
    !----------------
    if (associated(input%properties%wfplot%plot2d)) then

      igrid(1:2) = input%properties%wfplot%plot2d%parallelogram%grid(1:2)
      boxl(1,:)  = input%properties%wfplot%plot2d%parallelogram%origin%coord(1:3)
      boxl(2,:)  = input%properties%wfplot%plot2d%parallelogram%pointarray(1)%point%coord(1:3)-boxl(1,:)
      boxl(3,:)  = input%properties%wfplot%plot2d%parallelogram%pointarray(2)%point%coord(1:3)-boxl(1,:)
      ! test whether box is reasonable ?

      ! rgrid constructor
      grid = gen_2d_rgrid(igrid(1:2), boxl(1:3,:), 0)
      !call print_rgrid(grid)

      ! Generate WF on the grid
      allocate(zdata(grid%npt))
      call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
      call delete_rgrid(grid)

      if (rank==0) then
        write(fname,'("wf2d-",i,"-",i,".xsf")') ik, ist
        call str_strip(fname)
        call write_structure_xsf(fname)
        call write_2d_xsf(fname, 'module squared', boxl(1:3,:), igrid, grid%npt, abs(zdata)**2)
        call write_2d_xsf(fname, 'real',             boxl(1:3,:), igrid, grid%npt, dble(zdata))
        call write_2d_xsf(fname, 'imaginary',        boxl(1:3,:), igrid, grid%npt, aimag(zdata))
        write(*,*)
        write(*,'("Info(wfplot):")')
        write(*,'(" 2D wavefunction  written to wf2d-ik-ist.xsf")')
        write(*,*)
        write(*,'(" for k-point ", I6, " and state ", I6)') ik, ist
        write(*,*)
      end if

      deallocate(zdata)

    end if

    !----------------
    ! 3D case
    !----------------
    if (associated(input%properties%wfplot%plot3d)) then

      igrid(:)  = input%properties%wfplot%plot3d%box%grid
      boxl(1,:) = input%properties%wfplot%plot3d%box%origin%coord
      boxl(2,:) = input%properties%wfplot%plot3d%box%pointarray(1)%point%coord-boxl(1,:)
      boxl(3,:) = input%properties%wfplot%plot3d%box%pointarray(2)%point%coord-boxl(1,:)
      boxl(4,:) = input%properties%wfplot%plot3d%box%pointarray(3)%point%coord-boxl(1,:)

      ! rgrid constructor
      grid = gen_3d_rgrid(igrid, boxl(1:4,:), 0)
      !call print_rgrid(grid)

      ! Generate WF on the grid
      allocate(zdata(grid%npt))
      call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
      call delete_rgrid(grid)

      if (rank==0) then
        write(fname,'("wf3d-",i,"-",i,".xsf")') ik, ist
        call str_strip(fname)
        call write_structure_xsf(fname)
        call write_3d_xsf(fname, 'squared modulus', boxl(1:4,:), igrid, grid%npt, abs(zdata)**2)
        call write_3d_xsf(fname, 'real',            boxl(1:4,:), igrid, grid%npt, dble(zdata))
        call write_3d_xsf(fname, 'imaginary',       boxl(1:4,:), igrid, grid%npt, aimag(zdata))
        write(*,*)
        write(*,'("Info(wfplot):")')
        write(*,'(" 3D wavefunction written to wf3d-ik-ist.xsf")')
        write(*,*)
        write(*,'(" for k-point ", I6, " and state ", I6)') ik, ist
        write(*,*)
        !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))

        ! Gaussian cube-format
        write(fname,'("wf3d-",i,"-",i,".cube")') ik, ist
        call str_strip(fname)
        call write_3d_cube(fname, 'squared modulus', boxl(1:4,:), igrid, grid%npt, abs(zdata)**2)
      end if

      deallocate(zdata)

    end if

    deallocate(apwalm, evecfv, evecsv, wfmt, wfir)
   
    return
end subroutine

