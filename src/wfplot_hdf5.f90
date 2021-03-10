module wfplot_hdf5
    use modmain
    use modinput
    use modplotlabels
    use mod_rgrid
    use mod_xsf_format
    use mod_cube_format
    use modmpi, only : rank
    use mod_hdf5
    use modxs, only: isreadstate0
  contains
    subroutine write_wfplot_hdf5(ik,ist, plot1d_, plot2d_, plot3d_)
        implicit none
        ! input/output
        integer, intent(in) :: ik, ist
        type(plot1d_type), optional :: plot1d_
        type(plot2d_type), optional :: plot2d_
        type(plot3d_type), optional :: plot3d_
        ! local variables
        integer :: ip, np, iv, nv
        character(80) :: fname
        integer :: igrid(3)
        real(8) :: boxl(4,3)
        character(256) :: gname, cik
        character(256) :: basename="wfplot"
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
#ifdef XS
        call init2
#else
        ! read the density and potentials from file
        call readstate
        ! find the new linearisation energies
        call linengy
        ! generate the APW radial functions
        call genapwfr
        ! generate the local-orbital radial functions
        call genlofr
#endif
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
            write (*, '("Error(wfplot): nstsv : ", I8)') nstsv
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
        If (present(plot1d_)) then
          nv = size(plot1d_%path%pointarray)
          !write(*,*) nv
          if (nv < 1) then
            if (rank==0) then
              write (*,*)
              write (*,*) "Error(wfplot_new): Wrong plot specification!"
              write (*,*)
            end if
            stop
          end if
          np = plot1d_%path%steps
          !write(*,*) np
          If (np < nv) then
            if (rank==0) then
              write (*,*)
              write (*,*) "Error(wfplot_new): Wrong plot specification!"
              write (*,*)
            end if
            stop
          end if

          ! rgrid constructor
          grid = gen_1d_rgrid(plot1d_)
          !call print_rgrid(grid)

          ! Generate WF on the grid
          allocate(zdata(grid%npt))
          call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)

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
        if (present(plot2d_)) then

          ! rgrid constructor
          grid = gen_2d_rgrid(plot2d_, 0)
          !call print_rgrid(grid)

          ! Generate WF on the grid
          allocate(zdata(grid%npt))
          call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)

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
        if (present(plot3d_)) then

          ! rgrid constructor
          grid = gen_3d_rgrid(plot3d_, 0)
          !call print_rgrid(grid)

          ! Generate WF on the grid
          allocate(zdata(grid%npt))
          call calc_zdata_rgrid(grid, ik, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdata)
#ifdef _HDF5_
          if (rank==0) then
            if (.not. hdf5_exist_group(fhdf5, '/', basename)) then
               call hdf5_create_group(fhdf5, '/', basename)
            end if
            write(gname, '(I4.4)') ist
            gname="/"//trim(adjustl(basename))//"/"//trim(adjustl(gname))//"/"
            ! create group for state if necessary
            if (.not. hdf5_exist_group(fhdf5,'/',gname)) then
              call hdf5_create_group(fhdf5,'/',gname)
            end if
            gname=trim(adjustl(gname))//"/"
            ! create subgroup for k-point
            write(cik, '(I4.4)') ik
            if (.not. hdf5_exist_group(fhdf5,gname,cik)) then
              call hdf5_create_group(fhdf5,gname,cik)
            end if
            gname=trim(adjustl(gname))//"/"//trim(adjustl(cik))//"/"
            ! write data
            call hdf5_write(fhdf5,gname,'data',zdata(1),shape(zdata))
          end if
#endif
          call delete_rgrid(grid)
          deallocate(zdata)
          
        end if

        deallocate(apwalm, evecfv, evecsv, wfmt, wfir)
       
        return
    end subroutine
    
    subroutine write_box_hdf5(plot3d_)
      use modinput, only: plot3d_type
      use modmpi, only: rank
      implicit none
      !local variables 
      type(plot3d_type), intent(in) :: plot3d_
      character(256) :: gname
      integer:: i
      type(rgrid) :: grid
      real(8) :: boxc(4,3)
      real(8), parameter :: bohr2ang= 0.529177d0
#ifdef _HDF5_
      if (rank == 0) then
        gname=trim(adjustl("box"))
        if (.not. hdf5_exist_group(fhdf5,'/',gname)) then
          call hdf5_create_group(fhdf5,'/',gname)
        end if
        
        ! create cartesian coordinates for origin and axes
        grid=gen_3d_rgrid(plot3d_, 0)
        do i=1,4
          print *, 'grid%boxl(',i,'=', grid%boxl(i,:)
        end do
        do i=1,4
          call r3mv(input%structure%crystal%basevect, grid%boxl(i,:), boxc(i,:))
        end do
        boxc(:,:)=boxc(:,:)*bohr2ang
        gname=trim(adjustl(gname))//"/"
        call hdf5_write(fhdf5,gname,'origin',boxc(1,1),shape(boxc(1,:)))
        call hdf5_write(fhdf5,gname,'ax1',boxc(2,1),shape(boxc(2,:)))
        call hdf5_write(fhdf5,gname,'ax2',boxc(3,1),shape(boxc(3,:)))
        call hdf5_write(fhdf5,gname,'ax3',boxc(4,1),shape(boxc(4,:)))
        call hdf5_write(fhdf5,gname,'ngrid',grid%ngrid(1),shape(grid%ngrid))
      end if
#endif
    end subroutine

end module
