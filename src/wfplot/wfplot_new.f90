
subroutine wfplot_new(ik,ist)
  use modmain
  use modinput
  use modplotlabels
  use mod_rgrid
  use mod_xsf_format
  use mod_cube_format
  use modmpi, only: terminate_if_false, rank
  use wfplot_nice, only: setup_wfplot_gloabls, setup_wfplot_k, calculate_wfplot
  use mod_kpoint, only: vkl 

  implicit none
  ! input/output
  integer, intent(in) :: ik, ist
  ! local variables
  integer :: ip, np, iv, nv
  character(80) :: fname
  integer :: igrid(3)
  real(8) :: boxl(4,3)
  complex(8), allocatable :: evecfv(:, :), evecsv(:, :), apwalm(:, :, :, :)
  complex(8), allocatable :: zdata(:)
  !
  type(rgrid) :: grid
  type(plot1d_type), pointer :: plotdef
  type(plotlabels), pointer :: labels
  
  call setup_wfplot_gloabls(xs_calculation=.false.)
  allocate(evecfv(nmatmax, nstfv))
  allocate(evecsv(nstsv, nstsv))
  allocate(apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))

  !----------------
  ! 1D case
  !----------------
  If (associated(input%properties%wfplot%plot1d)) then
    nv = size(input%properties%wfplot%plot1d%path%pointarray)
    np = input%properties%wfplot%plot1d%path%steps
    call terminate_if_false(1 <= nv, "Error(wfplot_new): Wrong plot specification!")
    call terminate_if_false(nv <= np, "Error(wfplot_new): Wrong plot specification!")

    grid = gen_1d_rgrid(input%properties%wfplot%plot1d)
    call setup_wfplot_k(ik, evecfv, evecsv, apwalm)
    zdata = calculate_wfplot(ik, ist, grid, evecfv, evecsv, apwalm)
    print*, vkl(:, ik)

    ! Output
    if (rank==0) then
      write(fname,'("wf1d-",i4.4,"-",i4.4,".dat")') ik, ist
      open(77,file=trim(fname),status='Unknown',action='Write')
      do ip = 1, grid%npt
        ! path, |psi|^2, Re(psi), Im(psi) 
        write(77,'(4f16.6)') grid%vpd(ip), abs(zdata(ip))**2, zdata(ip)
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

    call delete_rgrid(grid)
    deallocate(zdata)
  end if

  !----------------
  ! 2D case
  !----------------
  if (associated(input%properties%wfplot%plot2d)) then
    grid = gen_2d_rgrid(input%properties%wfplot%plot2d, 0)
    call setup_wfplot_k(ik, evecfv, evecsv, apwalm)
    zdata = calculate_wfplot(ik, ist, grid, evecfv, evecsv, apwalm)
    print*, vkl(:, ik)

    if (rank==0) then
      write(fname,'("wf2d-",i4.4,"-",i4.4,".xsf")') ik, ist
      call str_strip(fname)
      call write_structure_xsf(fname)
      call write_2d_xsf(fname, 'module squared',   grid%boxl(1:3,:), grid%ngrid, grid%npt, abs(zdata)**2)
      call write_2d_xsf(fname, 'real',             grid%boxl(1:3,:), grid%ngrid, grid%npt, dble(zdata))
      call write_2d_xsf(fname, 'imaginary',        grid%boxl(1:3,:), grid%ngrid, grid%npt, aimag(zdata))
      write(*,*)
      write(*,'("Info(wfplot):")')
      write(*,'(" 2D wavefunction  written to wf2d-ik-ist.xsf")')
      write(*,*)
      write(*,'(" for k-point ", I6, " and state ", I6)') ik, ist
      write(*,*)
    end if

    call delete_rgrid(grid)
    deallocate(zdata)

  end if

  !----------------
  ! 3D case
  !----------------
  if (associated(input%properties%wfplot%plot3d)) then

    grid = gen_3d_rgrid(input%properties%wfplot%plot3d, 0)
    call setup_wfplot_k(ik, evecfv, evecsv, apwalm)
    zdata = calculate_wfplot(ik, ist, grid, evecfv, evecsv, apwalm)
    print*, vkl(:, ik)

    if (rank==0) then
      write(fname,'("wf3d-",i4.4,"-",i4.4,".xsf")') ik, ist
      call str_strip(fname)
      call write_structure_xsf(fname)
      call write_3d_xsf(fname, 'squared modulus', grid%boxl(1:4,:), grid%ngrid, grid%npt, abs(zdata)**2)
      call write_3d_xsf(fname, 'real',            grid%boxl(1:4,:), grid%ngrid, grid%npt, dble(zdata))
      call write_3d_xsf(fname, 'imaginary',       grid%boxl(1:4,:), grid%ngrid, grid%npt, aimag(zdata))
      write(*,*)
      write(*,'("Info(wfplot):")')
      write(*,'(" 3D wavefunction written to wf3d-ik-ist.xsf")')
      write(*,*)
      write(*,'(" for k-point ", I6, " and state ", I6)') ik, ist
      write(*,*)
      !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))

      ! Gaussian cube-format
      write(fname,'("wf3d-",i4.4,"-",i4.4,".cube")') ik, ist
      call str_strip(fname)
      call write_3d_cube(fname, 'squared modulus', grid%boxl(1:4,:), grid%ngrid, grid%npt, abs(zdata)**2)
    end if

    call delete_rgrid(grid)
    deallocate(zdata)
  end if
end subroutine
