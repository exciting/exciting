
subroutine wannier_plot( fst, lst, cell)
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
    integer, intent(in) :: fst, lst, cell(3)
    ! local variables
    integer :: ip, np, iv, nv, ik, iknr, ist, jst, is, ia, ias, maxdim, l, m, lm, o, ilo1, ngknr, maxnpt
    character(80) :: fname
    integer :: igrid(3)
    real(8) :: boxl(4,3), cellc(3), s, phi, v1(3), v2(3), v3(3)
    integer :: lmcnt( 0:input%groundstate%lmaxapw, nspecies), &
               o2idx( apwordmax, 0:input%groundstate%lmaxapw, nspecies), &
               lo2idx( nlomax, 0:input%groundstate%lmaxapw, nspecies)
    complex(8) :: phase
    ! allocatable arrays
    real(8), allocatable :: vvl(:,:), dist(:)
!    complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), dmatk(:,:,:,:,:), dmatr(:,:,:,:)
    complex(8), allocatable :: evecfv(:,:), apwalm(:,:,:,:), evecsv(:,:), evec(:)
    complex(8), allocatable :: wfmt(:,:,:,:), wfir(:,:)
    complex(8), allocatable :: zdata(:), zdatatot(:,:)
    !
    type(rgrid), allocatable :: grid(:)
    type(plot1d_type), pointer :: plotdef
    type(plotlabels), pointer :: labels

    ! initialise universal variables
    input%groundstate%lradstep = 1
    ! WARNING: make sure that the correct state is read and correct global k-point arrays are definded
    !call init0
    !call init1    
    fname = filext
    if( input%properties%wannier%input .eq. "hybrid") filext = '_PBE.OUT'
    call readstate
    filext = fname
    call readfermi
    call linengy
    call genapwfr
    call genlofr
    call olprad

    if ((fst<wf_fst) .or. (lst>wf_lst)) then
      if (rank==0) then
        write (*,*)
        write (*, '("Error (wannierplot): state out of range : ", I8)') ist
        write (*,*)
      end if
      stop
    end if

    call r3mv( input%structure%crystal%basevect, dble( cell), cellc)
    allocate( grid( fst:lst))

    ! generate plotting grids
    v1 = input%properties%wannierplot%plot3d%box%pointarray(1)%point%coord
    v2 = input%properties%wannierplot%plot3d%box%pointarray(2)%point%coord
    v3 = input%properties%wannierplot%plot3d%box%pointarray(3)%point%coord
    s = input%structure%crystal%scale
    do ist = fst, lst
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

        grid( ist) = gen_1d_rgrid(input%properties%wannierplot%plot1d)
      end if

      if (associated(input%properties%wannierplot%plot2d)) then
        grid( ist) = gen_2d_rgrid(input%properties%wannierplot%plot2d, 0)
      end if

      if( associated( input%properties%wannierplot%plot3d)) then
        call r3mv( ainv, wf_centers( :, ist-wf_fst+1) + cellc - s*(v1+v2+v3), input%properties%wannierplot%plot3d%box%origin%coord)
        call r3mv( ainv, wf_centers( :, ist-wf_fst+1) + cellc + s*(v1-v2-v3), input%properties%wannierplot%plot3d%box%pointarray(1)%point%coord)
        call r3mv( ainv, wf_centers( :, ist-wf_fst+1) + cellc + s*(v2-v3-v1), input%properties%wannierplot%plot3d%box%pointarray(2)%point%coord)
        call r3mv( ainv, wf_centers( :, ist-wf_fst+1) + cellc + s*(v3-v1-v2), input%properties%wannierplot%plot3d%box%pointarray(3)%point%coord)
        grid( ist) = gen_3d_rgrid( input%properties%wannierplot%plot3d, 0)
      end if
    end do

    maxdim = 0
    do is = 1, nspecies
      o2idx( :, :, is) = 0
      lo2idx( :, :, is) = 0
      lmcnt( :, is) = 0
      do l = 0, input%groundstate%lmaxapw
        do o = 1, apword( l, is)
          lmcnt( l, is) = lmcnt( l, is) + 1
          o2idx( o, l, is) = lmcnt( l, is)
        end do
      end do
      do ilo1 = 1, nlorb( is)
        l = lorbl( ilo1, is)
        lmcnt( l, is) = lmcnt( l, is) + 1
        lo2idx( ilo1, l, is) = lmcnt( l, is)
      end do
      maxdim = max( maxdim, maxval( lmcnt( :, is)))
    end do

    maxnpt = 0
    do ist = fst, lst
      maxnpt = max( jst, grid( ist)%npt)
    end do
    allocate( zdatatot( maxnpt, fst:lst))
    zdatatot(:,:) = zzero
    ! calculate the Wannier function on the grid


#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ip, ik, evecfv, evecsv, apwalm, evec, wfmt, wfir, zdata, phase, is, ia, ias)
#endif
    allocate( zdata( maxnpt))
    allocate( wfmt( lmmaxapw, nrmtmax, natmtot, nspinor))
    allocate( wfir( ngrtot, nspinor))
    allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
    allocate( evecfv( nmatmax, nstfv))
    allocate( evec( nmatmax))
#ifdef USEOMP
!$OMP DO
#endif
    do iknr = 1, nkptnr
      call myfindkpt( vklnr( :, iknr), wf_kset, ip, ik)
#ifdef USEOMP
!$OMP CRITICAL (readevec)
#endif
      call wannier_getevec( ik, evecfv)
#ifdef USEOMP
!$OMP END CRITICAL (readevec)
#endif
      
      ! find the matching coefficients
      call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)
      
      call ws_weight( dble( cell), dble( cell), wf_kset%vkl( :, ik), phase, .true.)

      do ist = fst, lst
        call zgemv( 'n', nmatmax, wf_nwf, zone, &
               evecfv( :, wf_fst:wf_lst), nmatmax, &
               wf_transform( :, ist-wf_fst+1, ik), 1, zzero, &
               evec, 1)
        
        ! calculate the wavefunctions for all states
        wfmt(:,:,:,:) = zzero
        wfir(:,:) = zzero
        zdata(:) = zzero
        do is = 1, nspecies
          do ia = 1, natoms( is)
            ias = idxas( ia, is)
            call wavefmt( 1, input%groundstate%lmaxapw, is, ia, wf_Gkset%ngk( 1, ik), apwalm, evec, lmmaxapw, wfmt( :, :, ias, 1))
          end do
        end do
        wfir( 1:wf_Gkset%ngk( 1, ik), 1) = evec( 1:wf_Gkset%ngk( 1, ik))/dsqrt( omega)
        call calc_zdata_rgrid( grid( ist), iknr, wfmt(:,:,:,1), wfir(:,1), zdata( 1:grid( ist)%npt), nosym=.true.)
#ifdef USEOMP
!$OMP CRITICAL (adddata)
#endif
        zdatatot( 1:grid( ist)%npt, ist) = zdatatot( 1:grid( ist)%npt, ist) + conjg( phase)*zdata( 1:grid( ist)%npt)
#ifdef USEOMP
!$OMP END CRITICAL (adddata)
#endif
      end do

    end do        
#ifdef USEOMP
!$OMP END DO 
#endif
    deallocate( zdata, wfmt, wfir, apwalm, evecfv, evec)
#ifdef USEOMP
!$OMP END PARALLEL 
#endif

    zdatatot = zdatatot/wf_kset%nkpt

    ist = fst
    ! phase correction
    allocate( dist( maxnpt))
    do ist = fst, lst
      phi = 0.d0
      do ip = 1, grid( ist)%npt
        dist( ip) = norm2( grid( ist)%vpc( :, ip) - wf_centers( :, ist-wf_fst+1) - cellc(:))
        s = atan2( aimag( zdatatot( ip, ist)), dble( zdatatot( ip, ist)))
        !write(*,'(F23.16)') s
        if( s .gt. 0.5d0*pi) s = s - pi
        if( s .lt. -0.5d0*pi) s = s + pi
        !write(*,'(I,F13.6)') ip, s
        phi = phi + abs( zdatatot( ip, ist))*s
      end do
      phi = phi/sum( abs( zdatatot( :, ist)))
      ip = minloc( dist, 1) 
      s = atan2( aimag( zdatatot( ip, ist)), dble( zdatatot( ip, ist)))
      if( abs( phi - s) .gt. 0.5d0*pi) phi = mod( phi+pi, pi)
      write(*,*) ip
      write(*,'(3F13.6)') wf_centers( :, ist-wf_fst+1) + cellc
      write(*,'(3F13.6)') grid( ist)%vpc( :, ip)
      write(*,'(3F13.6)') minval( grid( ist)%vpc( 1, :)), minval( grid( ist)%vpc( 2, :)), minval( grid( ist)%vpc( 3, :))
      write(*,'(3F13.6)') maxval( grid( ist)%vpc( 1, :)), maxval( grid( ist)%vpc( 2, :)), maxval( grid( ist)%vpc( 3, :))
      phase = exp( -zi*phi)
      zdatatot( :, ist) = phase*zdatatot( :, ist)
      write(*,'(2F13.6)') minval( dble( zdatatot( :, ist))), maxval( dble( zdatatot( :, ist)))
      write(*,'(2F13.6)') minval( aimag( zdatatot( :, ist))), maxval( aimag( zdatatot( :, ist)))
    end do

    !----------------
    ! 1D case
    !----------------
    if (associated(input%properties%wannierplot%plot1d)) then
      ! Output
      if (rank==0) then
        do ist = fst, lst
          write(fname,'("wannier1d-",i4.4,".dat")') ist
          open(77,file=trim(fname),status='Unknown',action='Write')
          do ip = 1, grid( ist)%npt
            ! path, |psi|^2, Re(psi), Im(psi) 
            write(77,'(4f16.6)') grid( ist)%vpd(ip), abs(zdatatot(ip, ist))**2, zdatatot(ip, ist)
            !write(77,'(2f16.6)') grid%vpd(ip), wkpt(ik)*nkptnr*abs(zdata(ip))**2
          end do
          close(77)
          write(*,*)
          write(*,'("Info(wannierplot):")')
          write(*,'(" 1D Wavefunction written to wannier1d-ist.dat")')
          write(*,*)
          write(*,'(" for state ", I6)') ist
          write(*,*)
        end do
      end if

      call delete_rgrid( grid( ist))
    end if

    !----------------
    ! 2D case
    !----------------
    if (associated(input%properties%wannierplot%plot2d)) then
      if (rank==0) then
        do ist = fst, lst
          write(fname,'("wannier2d-",i4.4,".xsf")') ist
          call str_strip(fname)
          call write_structure_xsf(fname)
          call write_2d_xsf(fname, 'module squared',   grid( ist)%boxl(1:3,:), grid( ist)%ngrid, grid( ist)%npt, abs(zdatatot( :, ist))**2)
          call write_2d_xsf(fname, 'real',             grid( ist)%boxl(1:3,:), grid( ist)%ngrid, grid( ist)%npt, dble(zdatatot( :, ist)))
          call write_2d_xsf(fname, 'imaginary',        grid( ist)%boxl(1:3,:), grid( ist)%ngrid, grid( ist)%npt, aimag(zdatatot( :, ist)))
          write(*,*)
          write(*,'("Info(wannierplot):")')
          write(*,'(" 2D wavefunction  written to wannier2d-ist.xsf")')
          write(*,*)
          write(*,'(" for state ", I6)') ist
          write(*,*)
        end do
      end if

      call delete_rgrid( grid( ist))
    end if

    !----------------
    ! 3D case
    !----------------
    if (associated(input%properties%wannierplot%plot3d)) then
      if (rank==0) then
        do ist = fst, lst
          write(fname,'("wannier3d-",i4.4,".xsf")') ist
          call str_strip(fname)
          call write_structure_xsf(fname)
          call write_3d_xsf(fname, 'squared modulus', grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, abs(zdatatot( :, ist))**2)
          call write_3d_xsf(fname, 'real',            grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, dble(zdatatot( :, ist)))
          call write_3d_xsf(fname, 'imaginary',       grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, aimag(zdatatot( :, ist)))
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
          call write_3d_cube(fname, 'squared modulus', grid( ist)%boxl(1:4,:), grid( ist)%ngrid, grid( ist)%npt, abs(zdatatot( :, ist))**2)
          call delete_rgrid( grid( ist))
        end do
      end if
    end if

    deallocate( zdatatot, grid)
   
    return
end subroutine

