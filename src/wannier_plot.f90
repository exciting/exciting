
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
    integer :: ip, np, iv, nv, ik, iknr, jst, is, ia, ias, maxdim, l, m, lm, o, ilo1, ngknr
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
    complex(8), allocatable :: zdata(:), zdatatot(:)
    !
    type(rgrid) :: grid
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

    if ((ist<wf_fst) .or. (ist>wf_lst)) then
      if (rank==0) then
        write (*,*)
        write (*, '("Error (wannierplot): state out of range : ", I8)') ist
        write (*,*)
      end if
      stop
    end if

    call r3mv( input%structure%crystal%basevect, dble( cell), cellc)

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

    allocate( zdatatot( grid%npt))
    zdatatot = zzero
    ! calculate the Wannier function on the grid

!    allocate( zdata( grid%npt))
!    allocate( wfmt( lmmaxapw, nrmtmax, natmtot, nspinor, wf_fst:wf_lst))
!    allocate( wfir( ngrtot, nspinor, wf_fst:wf_lst))
!
!      allocate( dmatk( wf_fst:wf_lst, maxdim, lmmaxapw, natmtot, wf_kset%nkpt))
!
!      dmatk = zzero
!
!      ! build k-point density coefficients
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, ngknr, apwalm, evecfv, is, ia, ias, l, m, lm, o, ilo1)
!#endif
!      allocate( evecfv( nmatmax, nstfv, nspnfv))
!      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
!#ifdef USEOMP                
!!$OMP DO
!#endif
!      do ik = 1, wf_kset%nkpt
!        ngknr = wf_Gkset%ngk( 1, ik)
!
!        ! get matching coefficients
!        call match( ngknr, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
!          
!#ifdef USEOMP                
!!$OMP CRITICAL (readevec)
!#endif
!        ! read eigenvector      
!        if( input%properties%wannier%input .eq. "groundstate") then
!          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
!        else if( input%properties%wannier%input .eq. "hybrid") then
!          call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
!        else if( input%properties%wannier%input .eq. "gw") then
!          call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspnfv, evecfv)
!        else
!          call terminate
!        end if
!#ifdef USEOMP                
!!$OMP END CRITICAL (readevec)
!#endif
!
!        do is = 1, nspecies
!          do ia = 1, natoms( is)
!            ias = idxas( ia, is)
!            ! APW contribution
!            do l = 0, input%groundstate%lmaxapw
!              do m = -l, l
!                lm = idxlm( l, m)
!                do o = 1, apword( l, is)
!                  call zgemv( 'T', ngknr, wf_nst, zone, &
!                       evecfv( 1:ngknr, wf_fst:wf_lst, 1), ngknr, &
!                       apwalm( 1:ngknr, o, lm, ias, 1), 1, zzero, &
!                       dmatk( :, o2idx( o, l, is), lm, ias, ik), 1)
!                end do
!              end do
!            end do
!            ! LO contribution
!            do ilo1 = 1, nlorb( is)
!              l = lorbl( ilo1, is)
!              do m = -l, l
!                lm = idxlm( l, m)
!                dmatk( :, lo2idx( ilo1, l, is), lm, ias, ik) = evecfv( ngknr+idxlo( lm, ilo1, ias), wf_fst:wf_lst, 1)
!              end do
!            end do
!          end do
!        end do
!      end do
!#ifdef USEOMP                
!!$OMP END DO
!#endif
!      deallocate( evecfv, apwalm)
!#ifdef USEOMP                
!!$OMP END PARALLEL
!#endif
!
!      ! build R-point density coefficients
!      allocate( dmatr( wf_fst:wf_lst, maxdim, lmmaxapw, natmtot))
!      dmatr = zzero
!      do ik = 1, wf_kset%nkpt
!        call ws_weight( dble( cell), dble( cell), wf_kset%vkl( :, ik), phase, kgrid=.true.)
!        do is = 1, nspecies
!          do ia = 1, natoms( is)
!            ias = idxas( ia, is)
!            do lm = 1, lmmaxapw
!              call zgemm( 'T', 'N', wf_nst, maxdim, wf_nst, conjg( phase), &
!                   wf_transform( :, :, ik), wf_nst, &
!                   dmatk( :, :, lm, ias, ik), wf_nst, zone, &
!                   dmatr( :, :, lm, ias), wf_nst)
!            end do
!          end do
!        end do
!      end do
!      
!      deallocate( dmatk)
!
!      call plotmat( transpose( dmatr( ist, :, :, 1)))
!      write(*,*)
!      call plotmat( transpose( dmatr( ist, :, :, 2)))
!
!      ! build radial functions
!      wfmt(:,:,:,:,:) = zzero
!      do is = 1, nspecies
!        do ia = 1, natoms( is)
!          ias = idxas( ia, is)
!          do l = 0, input%groundstate%lmaxapw
!            do m = -l, l
!              lm = idxlm( l, m)
!              do o = 1, apword( l, is)
!                wfmt( lm, :, ias, 1, ist) = wfmt( lm, :, ias, 1, ist) + dmatr( ist, o2idx( o, l, is), lm, ias)*apwfr( :, 1, o, l, ias)
!              end do
!            end do
!          end do
!          do ilo1 = 1, nlorb( is)
!            l = lorbl( ilo1, is)
!            do m = -l, l
!              lm = idxlm( l, m)
!              wfmt( lm, :, ias, 1, ist) = wfmt( lm, :, ias, 1, ist) + dmatr( ist, lo2idx( ilo1, l, is), lm, ias)*lofr( :, 1, ilo1, ias)
!            end do
!          end do
!        end do
!      end do
!
!      call calc_zdata_rgrid( grid, 1, wfmt(:,:,:,1,ist), wfir(:,1,ist), zdatatot, nosym=.true.)



    !call print_rgrid( grid)
    !stop


#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ip, ik, evecfv, evecsv, apwalm, evec, wfmt, wfir, zdata, phase, is, ia, ias)
#endif
    allocate( zdata( grid%npt))
    allocate( wfmt( lmmaxapw, nrmtmax, natmtot, nspinor))
    allocate( wfir( ngrtot, nspinor))
    allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
    allocate( evecfv( nmatmax, nstfv))
    allocate( evec( nmatmax))
    allocate( evecsv( nstsv, nstsv))
#ifdef USEOMP
!$OMP DO
#endif
    do iknr = 1, nkptnr
      call myfindkpt( vklnr( :, iknr), wf_kset, ip, ik)
#ifdef USEOMP
!$OMP CRITICAL (readevec)
#endif
      write(*,'(2(3F13.6,3x))') wf_kset%vkl( :, ik), vklnr( :, iknr)
      ! get the eigenvectors and values from file
      if( input%properties%wannier%input .eq. "groundstate") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspinor, evecfv)
      else
        call terminate
      end if
      call getevecsv( wf_kset%vkl( :, ik), evecsv)
#ifdef USEOMP
!$OMP END CRITICAL (readevec)
#endif
      
      ! find the matching coefficients
      call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm)
      
      call ws_weight( dble( cell), dble( cell), wf_kset%vkl( :, ik), phase, .true.)

      call zgemv( 'N', nmatmax, wf_nst, zone, &
           evecfv( :, wf_fst:wf_lst), nmatmax, &
           wf_transform( :, ist, ik), 1, zzero, &
           evec, 1)
      
      ! calculate the wavefunctions for all states
      wfmt(:,:,:,:) = zzero
      wfir(:,:) = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          call wavefmt( 1, input%groundstate%lmaxapw, is, ia, wf_Gkset%ngk( 1, ik), apwalm, evec, lmmaxapw, wfmt( :, :, ias, 1))
        end do
      end do
      wfir( 1:wf_Gkset%ngk( 1, ik), 1) = evec( 1:wf_Gkset%ngk( 1, ik))/dsqrt( omega)
      call calc_zdata_rgrid( grid, iknr, wfmt(:,:,:,1), wfir(:,1), zdata, nosym=.true.)
#ifdef USEOMP
!$OMP CRITICAL (adddata)
#endif
      zdatatot = zdatatot + conjg( phase)*zdata
#ifdef USEOMP
!$OMP END CRITICAL (adddata)
#endif

    end do        
#ifdef USEOMP
!$OMP END DO 
#endif
    deallocate( zdata, wfmt, wfir, apwalm, evecfv, evecsv, evec)
#ifdef USEOMP
!$OMP END PARALLEL 
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
      write(*,'(I,F13.6)') ip, s
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

