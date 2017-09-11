
subroutine wannier_plot( ist, cell)
    use modmain
    use modinput
    use modplotlabels
    use mod_rgrid
    use mod_xsf_format
    use mod_cube_format
    use mod_wannier
    use m_wsweight
    use mod_lattice
    use modmpi, only : rank
    implicit none
    ! input/output
    integer, intent(in) :: ist, cell(3)
    ! local variables
    integer :: ip, np, iv, nv, ik, iknr, jst
    character(80) :: fname
    integer :: igrid(3)
    real(8) :: boxl(4,3), cellc(3), s, phi, v1(3), v2(3), v3(3)
    complex(8) :: phase
    ! allocatable arrays
    real(8), allocatable :: vvl(:,:), dist(:)
    real(8), allocatable :: vgkltmp(:,:,:,:)
    complex(8), allocatable :: apwalm(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:)
    complex(8), allocatable :: evecsv(:,:)
    complex(8), allocatable :: wfmt(:,:,:,:,:)
    complex(8), allocatable :: wfir(:,:,:)
    complex(8), allocatable :: zdata(:), zdatatot(:)
    !
    type(rgrid) :: grid
    type(plot1d_type), pointer :: plotdef
    type(plotlabels), pointer :: labels

    ! initialise universal variables
    input%groundstate%lradstep = 1
    call init0
    call init1    
    call readstate
    call linengy
    call genapwfr
    call genlofr

    if ((ist<wf_fst) .or. (ist>wf_lst)) then
      if (rank==0) then
        write (*,*)
        write (*, '("Error (wannierplot): state out of range : ", I8)') ist
        write (*,*)
      end if
      stop
    end if

    call r3mv( input%structure%crystal%basevect, dble( cell), cellc)
    !call wannier_gen_max( input%properties%wannier%state, input%properties%wannier%nst, (/1, 3, 4, 5/))

    ! allocate local arrays

    ! generate plotting grid
    if (associated(input%properties%wannierplot%plot1d)) then
      nv = size(input%properties%wannierplot%plot1d%path%pointarray)
      if (nv < 1) then
        if (rank==0) then
          write (*,*)
          write (*,*) "Error (wannierplot): Wrong plot specification!"
          write (*,*)
        end if
        stop
      end if
      np = input%properties%wannierplot%plot1d%path%steps
      If (np < nv) then
        if (rank==0) then
          write (*,*)
          write (*,*) "Error (wannierplot): Wrong plot specification!"
          write (*,*)
        end if
        stop
      end if

      grid = gen_1d_rgrid(input%properties%wannierplot%plot1d)
    end if

    if (associated(input%properties%wannierplot%plot2d)) then
      grid = gen_2d_rgrid(input%properties%wannierplot%plot2d, 0)
    end if

    if( associated( input%properties%wannierplot%plot3d)) then
      v1 = input%properties%wannierplot%plot3d%box%pointarray(1)%point%coord
      v2 = input%properties%wannierplot%plot3d%box%pointarray(2)%point%coord
      v3 = input%properties%wannierplot%plot3d%box%pointarray(3)%point%coord
      s = input%structure%crystal%scale
      call r3mv( ainv, wf_centers( :, ist) + cellc - s*(v1+v2+v3), input%properties%wannierplot%plot3d%box%origin%coord)
      call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v1-v2-v3), input%properties%wannierplot%plot3d%box%pointarray(1)%point%coord)
      call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v2-v3-v1), input%properties%wannierplot%plot3d%box%pointarray(2)%point%coord)
      call r3mv( ainv, wf_centers( :, ist) + cellc + s*(v3-v1-v2), input%properties%wannierplot%plot3d%box%pointarray(3)%point%coord)
      grid = gen_3d_rgrid( input%properties%wannierplot%plot3d, 0)
    end if

    allocate( zdatatot( grid%npt))
    zdatatot = zzero
    ! calculate the Wannier function on the grid
    allocate( vgkltmp( 3, ngkmax, nspnfv, nkpt))

#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iknr, evecfv, evecsv, apwalm, wfmt, wfir, jst, zdata, phase) reduction(+:zdatatot)
#endif
    allocate( zdata( grid%npt))
    allocate( wfmt( lmmaxapw, nrmtmax, natmtot, nspinor, nstsv))
    allocate( wfir( ngrtot, nspinor, nstsv))
    allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
    allocate( evecfv( nmatmax, nstfv))
    allocate( evecsv( nstsv, nstsv))
#ifdef USEOMP
!!$OMP DO
#endif
    do iknr = 1, wf_kset%nkpt
      call findkpt( wf_kset%vkl( :, iknr), ip, ik)
      ! get the eigenvectors and values from file
      if( input%properties%wannier%input .eq. "groundstate") then
        call getevecfv( wf_kset%vkl( :, iknr), wf_Gkset%vgkl( :, :, :, iknr), evecfv)
        call getevecsv( wf_kset%vkl( :, iknr), evecsv)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, iknr), wf_Gkset%vgkl( :, :, :, iknr), evecfv)
        call getevecsv( wf_kset%vkl( :, iknr), evecsv)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", iknr, wf_kset%vkl( :, iknr), nmatmax, nstfv, nspinor, evecfv)
        call getevecsv( wf_kset%vkl( :, iknr), evecsv)
      else
        call terminate
      end if
      
      ! find the matching coefficients
      call match( wf_Gkset%ngk( 1, iknr), wf_Gkset%gkc( :, 1, iknr), wf_Gkset%tpgkc( :, :, 1, iknr), wf_Gkset%sfacgk( :, :, 1, iknr), apwalm)
      
      ! calculate the wavefunctions for all states
      call genwfsv_new( ik, 1, nstsv, apwalm, evecfv, evecsv, wfmt, wfir)
      
      call ws_weight( dble( cell), dble( cell), wf_kset%vkl( :, iknr), phase, .true.)

      do jst = wf_fst, wf_lst
        call calc_zdata_rgrid( grid, iknr, wfmt(:,:,:,1,jst), wfir(:,1,jst), zdata, nosym=.true.)
        zdatatot = zdatatot + phase*wf_transform( jst, ist, iknr)*zdata
      end do
    end do        
#ifdef USEOMP
!!$OMP END DO 
#endif
    deallocate( zdata, wfmt, wfir, apwalm, evecfv, evecsv)
#ifdef USEOMP
!!$OMP END PARALLEL
#endif
    zdatatot = zdatatot/wf_kset%nkpt
    ! phase correction
    allocate( dist( grid%npt)) 
    phi = 0.d0
    do ip = 1, grid%npt
      dist( ip) = norm2( grid%vpc( :, ip) - wf_centers( :, ist) - cellc(:))
      s = atan2( aimag( zdatatot( ip)), dble( zdatatot( ip)))
      !write(*,'(F23.16)') s
      if( s .gt. 0.5d0*pi) s = s - pi
      if( s .lt. -0.5d0*pi) s = s + pi
      phi = phi + abs( zdatatot( ip))*s
    end do
    phi = phi/sum( abs( zdatatot))
    ip = minloc( dist, 1) 
    s = atan2( aimag( zdatatot( ip)), dble( zdatatot( ip)))
    if( abs( phi - s) .gt. 0.5d0*pi) phi = mod( phi+pi, pi)
    write(*,*) ip
    write(*,'(3F13.6)') wf_centers( :, ist) + cellc
    write(*,'(3F13.6)') grid%vpc( :, ip)
    write(*,'(3F13.6)') minval( grid%vpc( 1, :)), minval( grid%vpc( 2, :)), minval( grid%vpc( 3, :))
    write(*,'(3F13.6)') maxval( grid%vpc( 1, :)), maxval( grid%vpc( 2, :)), maxval( grid%vpc( 3, :))
    phase = exp( -zi*phi)
    zdatatot = phase*zdatatot
    write(*,'(2F13.6)') minval( dble( zdatatot)), maxval( dble( zdatatot))
    write(*,'(2F13.6)') minval( aimag( zdatatot)), maxval( aimag( zdatatot))

    !----------------
    ! 1D case
    !----------------
    if (associated(input%properties%wannierplot%plot1d)) then
      ! Output
      if (rank==0) then
        write(fname,'("wannier1d-",i4.4,".dat")') ist
        open(77,file=trim(fname),status='Unknown',action='Write')
        do ip = 1, grid%npt
          ! path, |psi|^2, Re(psi), Im(psi) 
          write(77,'(4f16.6)') grid%vpd(ip), abs(zdatatot(ip))**2, zdatatot(ip)
          !write(77,'(2f16.6)') grid%vpd(ip), wkpt(ik)*nkptnr*abs(zdata(ip))**2
        end do
        close(77)
        write(*,*)
        write(*,'("Info(wannierplot):")')
        write(*,'(" 1D Wavefunction written to wannier1d-ist.dat")')
        write(*,*)
        write(*,'(" for state ", I6)') ist
        write(*,*)
      end if

      call delete_rgrid( grid)
      deallocate(zdatatot)
      
    end if

    !----------------
    ! 2D case
    !----------------
    if (associated(input%properties%wannierplot%plot2d)) then
      if (rank==0) then
        write(fname,'("wannier2d-",i4.4,".xsf")') ist
        call str_strip(fname)
        call write_structure_xsf(fname)
        call write_2d_xsf(fname, 'module squared',   grid%boxl(1:3,:), grid%ngrid, grid%npt, abs(zdatatot)**2)
        call write_2d_xsf(fname, 'real',             grid%boxl(1:3,:), grid%ngrid, grid%npt, dble(zdatatot))
        call write_2d_xsf(fname, 'imaginary',        grid%boxl(1:3,:), grid%ngrid, grid%npt, aimag(zdatatot))
        write(*,*)
        write(*,'("Info(wannierplot):")')
        write(*,'(" 2D wavefunction  written to wannier2d-ist.xsf")')
        write(*,*)
        write(*,'(" for state ", I6)') ist
        write(*,*)
      end if

      call delete_rgrid( grid)
    end if

    !----------------
    ! 3D case
    !----------------
    if (associated(input%properties%wannierplot%plot3d)) then
      if (rank==0) then
        write(fname,'("wannier3d-",i4.4,".xsf")') ist
        call str_strip(fname)
        call write_structure_xsf(fname)
        call write_3d_xsf(fname, 'squared modulus', grid%boxl(1:4,:), grid%ngrid, grid%npt, abs(zdatatot)**2)
        call write_3d_xsf(fname, 'real',            grid%boxl(1:4,:), grid%ngrid, grid%npt, dble(zdatatot))
        call write_3d_xsf(fname, 'imaginary',       grid%boxl(1:4,:), grid%ngrid, grid%npt, aimag(zdatatot))
        write(*,*)
        write(*,'("Info(wannierplot):")')
        write(*,'(" 3D wavefunction written to wannier3d-ist.xsf")')
        write(*,*)
        write(*,'(" for state ", I6)') ist
        write(*,*)
        !call write_supercell_xsf('supercell.xsf',(/-2,2/),(/-2,2/),(/-2,2/))

        ! Gaussian cube-format
        write(fname,'("wannier3d-",i4.4,".cube")') ist
        call str_strip(fname)
        call write_3d_cube(fname, 'squared modulus', grid%boxl(1:4,:), grid%ngrid, grid%npt, abs(zdatatot)**2)

        call delete_rgrid( grid)
      end if
    end if

    deallocate(zdatatot)
   
    return
end subroutine

