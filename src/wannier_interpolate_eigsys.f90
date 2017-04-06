module m_wannier_interpolate_eigsys
    implicit none
    contains

subroutine wannier_interpolate_eigsys( evalin, nqin, vqlin, evalout, evecout)
  use mod_wannier
  use m_wsweight
  use m_plotmat
  use mod_kqpts
  use mod_lattice
  use mod_eigenvalue_occupancy

  implicit none
  integer, intent( in) :: nqin
  real(8), intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
  real(8), intent( in) :: vqlin( 3, nqin)
  real(8), intent( out) :: evalout( wf_fst:wf_lst, nqin)
  complex(8), intent( out) :: evecout( nmatmax, wf_fst:wf_lst, nqin)
  
  integer :: nrpt, ix, iy, iz, ik, iq, ir, ngqwf, ngkwf, ngkmaxtmp
  complex(8) :: ftweight
  real(8) :: vqcin(3), v(3), ffac, t1

  real(8), allocatable :: rptl(:,:), vkctmp(:,:)
  complex(8), allocatable :: auxmat(:,:), ueu(:,:,:), phase(:,:), hamilton(:,:,:), evecint(:,:,:), evecfv(:,:,:), apwolp(:,:), uv(:,:), cuv(:,:), rhs(:,:)
  complex(8), allocatable :: evecq(:,:)
  
  integer :: o, is, ia, ias, l, m, lm, lmo, nlmomax, igq, igk
  integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:), lm2l(:)
  integer, allocatable :: igxigwf(:)
  real(8), allocatable :: vgxlwf(:,:,:), vgqcwf(:,:,:), gxcwf(:), tpgxcwf(:,:), sfacwf(:,:)
  real(8), allocatable :: vgkcwf(:,:,:)
  complex(8), allocatable :: matchq(:,:,:), matchk(:,:,:), apwalm(:,:,:,:,:)

  integer :: lapack_lwork, lapack_info
  integer, allocatable :: lapack_ipiv(:)
  complex(8), allocatable :: lapack_work(:)

  logical :: ongrid

  !**********************************************
  ! interpolated eigenenergies and corresponding 
  ! eigenvectors in Wannier basis
  !**********************************************

  ! generate set of lattice vectors 
  nrpt = wf_kset%nkpt
  allocate( rptl( 3, nrpt))
  ia = 0
  do iz = -wf_kset%ngridk(3)/2, -wf_kset%ngridk(3)/2+wf_kset%ngridk(3)-1
    do iy = -wf_kset%ngridk(2)/2, -wf_kset%ngridk(2)/2+wf_kset%ngridk(2)-1
      do ix = -wf_kset%ngridk(1)/2, -wf_kset%ngridk(1)/2+wf_kset%ngridk(1)-1
        ia = ia + 1
        rptl( :, ia) = (/ dble( ix), dble( iy), dble( iz)/)
      end do
    end do
  end do
  
  ! calculate Hamlitonian matrix elements in Wannier representation 
  allocate( ueu( wf_fst:wf_lst, wf_fst:wf_lst, wf_kset%nkpt))
  allocate( phase( wf_kset%nkpt, nqin))
  allocate( hamilton( wf_fst:wf_lst, wf_fst:wf_lst, nqin))
#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iy, auxmat, iq, ftweight) reduction(+:phase)
#endif
  allocate( auxmat( wf_fst:wf_lst, wf_fst:wf_lst))
#ifdef USEOMP
!!$OMP DO
#endif
  do ik = 1, wf_kset%nkpt
    do iy = wf_fst, wf_lst
      ueu( iy, :, ik) = wf_transform( iy, :, ik)*evalin( iy, ik)
    end do
    call zgemm( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
         wf_transform( :, :, ik), wf_nst, &
         ueu( :, :, ik), wf_nst, zzero, &
         auxmat, wf_nst)
    ueu( :, :, ik) = auxmat
    do iq = 1, nqin
      phase( ik, iq) = zzero
      do ir = 1, nrpt
        call ws_weight( rptl( :, ir), rptl( :, ir), vqlin( :, iq)-wf_kset%vkl( :, ik), ftweight, kgrid=.true.)
        phase( ik, iq) = phase( ik, iq) + conjg( ftweight)
      end do
    end do
  end do
#ifdef USEOMP
!!$OMP END DO
#endif
  deallocate( auxmat)
#ifdef USEOMP
!!$OMP END PARALLEL
#endif
  phase = phase/wf_kset%nkpt
  do iy = wf_fst, wf_lst
    call zgemm( 'N', 'N', wf_nst, nqin, wf_kset%nkpt, zone, &
         ueu( iy, :, :), wf_nst, &
         phase, wf_kset%nkpt, zzero, &
         hamilton( iy, :, :), wf_nst)
  end do
  deallocate( ueu)

  ! interpolation
  allocate( evecint( wf_fst:wf_lst, wf_fst:wf_lst, nqin))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq)
!$OMP DO
#endif
  do iq = 1, nqin 
    call diaghermat( wf_nst, hamilton( :, :, iq), evalout( :, iq), evecint( :, :, iq))
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  deallocate( rptl, hamilton)
  
  !**********************************************
  ! interpolated eigenvectors in LAPW+lo basis
  !**********************************************

  ! determine ngkmax for all k and q
  allocate( vkctmp( 3, nkpt))
  vkctmp = vkc
  ngkmaxtmp = ngkmax
  deallocate( vkc)
  allocate( vkc( 3, nqin))
  do iq = 1, nqin
    call r3mv( bvec, vqlin( :, iq), vkc( :, iq))
  end do
  call getngkmax
  deallocate( vkc)
  allocate( vkc( 3, nkpt))
  vkc = vkctmp
  write(*,*) ngkmax, ngkmaxtmp
  ngkmax = 2*max( ngkmax, ngkmaxtmp)

  ! generate APW radial functions
  call genapwfr

  ! count combined (l,m,o) indices and build index maps
  allocate( nlmo( nspecies))
  allocate( lmo2l( (input%groundstate%lmaxapw + 1)**2*apwordmax, nspecies), &
            lmo2m( (input%groundstate%lmaxapw + 1)**2*apwordmax, nspecies), &
            lmo2o( (input%groundstate%lmaxapw + 1)**2*apwordmax, nspecies))
  nlmomax = 0
  do is = 1, nspecies
    nlmo( is) = 0
    do l = 0, input%groundstate%lmaxapw
      do o = 1, apword( l, is)
        do m = -l, l
          nlmo( is) = nlmo( is) + 1
          lmo2l( nlmo( is), is) = l
          lmo2m( nlmo( is), is) = m
          lmo2o( nlmo( is), is) = o
        end do
      end do
    end do
    nlmomax = max( nlmomax, nlmo( is))
  end do

  allocate( evecfv( nmatmax, nstsv, nspinor))
  allocate( igxigwf( ngkmax))
  allocate( vgxlwf( 3, ngkmax, nspinor), vgqcwf( 3, ngkmax, nspinor), gxcwf( ngkmax), tpgxcwf( 2, ngkmax))
  allocate( vgkcwf( 3, ngkmax, nspinor))
  allocate( sfacwf( ngkmax, natmtot))
  allocate( uv( wf_fst:wf_lst, wf_fst:wf_lst))
  allocate( cuv( nmatmax, wf_fst:wf_lst))
  allocate( rhs( ngkmax, wf_fst:wf_lst))
  allocate( matchq( nlmomax, ngkmax, natmtot))
  allocate( matchk( nlmomax, ngkmax, natmtot))
  allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
  allocate( evecq( nmatmax, wf_fst:wf_lst))

  do iq = 1, nqin
    ! generate G+q vectors
    call r3mv( bvec, vqlin( :, iq), vqcin)
    call gengpvec( vqlin( :, iq), vqcin, ngqwf, igxigwf, vgxlwf(:,:,1), vgqcwf(:,:,1), gxcwf, tpgxcwf)
    if( ngqwf .gt. ngkmax) then
      write( *, '("ngkmax ",I)') iq
      ngqwf = ngkmax
    endif
    ngqwf = min( ngqwf, ngkmax)

    ! get matching coefficients for G+q        
    call gensfacgp( ngqwf, vgqcwf, ngkmax, sfacwf)
    call match( ngqwf, gxcwf, tpgxcwf, sfacwf, apwalm(:, :, :, :, 1))
    matchq = zzero
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        do lmo = 1, nlmo( is)
          l = lmo2l( lmo, is)
          m = lmo2m( lmo, is)
          o = lmo2o( lmo, is)
          lm = idxlm( l, m)
          matchq( lmo, 1:ngqwf, ias) = apwalm( 1:ngqwf, o, lm, ias, 1)
        end do
      end do
    end do

    rhs = zzero
    evecq = zzero
    ongrid = .false.

    write(*, '(I,3F13.6)') iq, vqlin( :, iq) 
    !write(*,*) "-----------------"
    do ik = 1, wf_kset%nkpt
      !write(*, '(I,3F13.6)') ik, wf_kset%vkl( :, ik) 

      ! generate G+k vectors
      call gengpvec( wf_kset%vkl( :, ik), wf_kset%vkc( :, ik), ngkwf, igxigwf, vgxlwf(:,:,1), vgkcwf(:,:,1), gxcwf, tpgxcwf)
      ! read eigenvector      
      if( input%properties%wannier%input .eq. "groundstate") then
        call getevecfv( wf_kset%vkl( :, ik), vgxlwf, evecfv)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, ik), vgxlwf, evecfv)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspinor, evecfv)
      else
        call terminate
      end if
    
      call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
           wf_transform( :, :, ik), wf_nst, &
           !wf_transform( :, :, ik), wf_nst, zzero, &
           !evecint( :, :, iq), wf_nst, &
           evecint( :, :, iq), wf_nst, zzero, &
           uv, wf_nst)
      
      if( norm2( wf_kset%vkl( :, ik) - vqlin( :, iq)) .lt. input%structure%epslat) then
      !  write( *, '(3F13.6)') vqlin( :, iq)
      !  write( *, '(3F13.6)') wf_kset%vkl( :, ik)
      !  call plotmat( uv)
      !  write( *, '(2F13.6)') phase( ik, iq)
        ongrid = .true.
        write(*,*) "original"
        call plotmat( evecfv( :, wf_fst:wf_lst, 1), .true.)
        write(*,*)
      end if

      call zgemm( 'N', 'N', nmatmax, wf_nst, wf_nst, zone, &
           evecfv( :, wf_fst:wf_lst, 1), nmatmax, &
           uv, wf_nst, zzero, &
           cuv, nmatmax)

      !----------!
      ! APW part !
      !----------!
      ! get matching coefficients for G+k        
      call gensfacgp( ngkwf, vgkcwf, ngkmax, sfacwf)
      call match( ngkwf, gxcwf, tpgxcwf, sfacwf, apwalm(:, :, :, :, 1))
      matchk = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do lmo = 1, nlmo( is)
            l = lmo2l( lmo, is)
            m = lmo2m( lmo, is)
            o = lmo2o( lmo, is)
            lm = idxlm( l, m)
            matchk( lmo, 1:ngkwf, ias) = apwalm( 1:ngkwf, o, lm, ias, 1)
          end do
        end do
      end do
      ! build APW q-k overlap
      if( allocated( apwolp)) deallocate( apwolp)
      allocate( apwolp( ngqwf, ngkwf))
        ! muffin tin region
      apwolp = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          call zgemm( 'C', 'N', ngqwf, ngkwf, nlmo( is), zone, &
               matchq( :, 1:ngqwf, ias), nlmo( is), &
               matchk( :, 1:ngkwf, ias), nlmo( is), zone, &
               apwolp, ngqwf)
        end do
      end do
        ! interstitial region
      do igk = 1, ngkwf
        do igq = 1, ngqwf
          v(:) = vgkcwf( :, igk, 1) - vgqcwf( :, igq, 1)
          if( norm2( v) .lt. input%structure%epslat) apwolp( igq, igk) = apwolp( igq, igk) + zone
          do is = 1, nspecies
            t1 = norm2( v)*rmt( is)
            if( norm2( v) .lt. input%structure%epslat) then
              ffac = fourpi/(3.d0*omega)*rmt( is)**3
            else
              ffac = fourpi/omega*(sin( t1) - t1*cos( t1))/norm2( v)**3
            end if
            do ia = 1, natoms( is)
              t1 = -dot_product( v, atposc( :, ia, is))
              apwolp( igq, igk) = apwolp( igq, igk) - cmplx( cos( t1), sin( t1), 8)*ffac
            end do
          end do
        end do
      end do
      ! sum up right hand side of linear system of equations
      call zgemm( 'N', 'N', ngqwf, wf_nst, ngkwf, phase( ik, iq), &
           apwolp, ngqwf, &
           cuv( 1:ngkwf, :), ngkwf, zone, &
           rhs( 1:ngqwf, :), ngqwf)

      !---------!
      ! LO part !
      !---------!
      !evecq( (ngqwf+1):(ngqwf+nlotot), :) = evecq( (ngqwf+1):(ngqwf+nlotot), :) + phase( ik, iq)*cuv( (ngqwf+1):(ngqwf+nlotot), :)
      evecq( (ngqwf+1):nmatmax, :) = evecq( (ngqwf+1):nmatmax, :) + phase( ik, iq)*cuv( (ngqwf+1):nmatmax, :)
      
    end do
    ! build APW q-q overlap
    if( allocated( apwolp)) deallocate( apwolp)
    allocate( apwolp( ngqwf, ngqwf))
      ! muffin tin region
    apwolp = zzero
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        call zgemm( 'C', 'N', ngqwf, ngqwf, nlmo( is), zone, &
             matchq( :, 1:ngqwf, ias), nlmo( is), &
             matchq( :, 1:ngqwf, ias), nlmo( is), zone, &
             apwolp, ngqwf)
      end do
    end do
      ! interstitial region
    do igk = 1, ngqwf
      do igq = 1, ngqwf
        v(:) = vgqcwf( :, igk, 1) - vgqcwf( :, igq, 1)
        if( norm2( v) .lt. input%structure%epslat) apwolp( igq, igk) = apwolp( igq, igk) + zone
        do is = 1, nspecies
          t1 = norm2( v)*rmt( is)
          if( norm2( v) .lt. input%structure%epslat) then
            ffac = fourpi/(3.d0*omega)*rmt( is)**3
          else
            ffac = fourpi/omega*(sin( t1) - t1*cos( t1))/norm2( v)**3
          end if
          do ia = 1, natoms( is)
            t1 = -dot_product( v, atposc( :, ia, is))
            apwolp( igq, igk) = apwolp( igq, igk) - cmplx( cos( t1), sin( t1), 8)*ffac
          end do
        end do
      end do
    end do
    ! solve linear system of equations
    if( allocated( lapack_ipiv)) deallocate( lapack_ipiv)
    allocate( lapack_ipiv( ngqwf))
    if( allocated( lapack_work)) deallocate( lapack_work)
    allocate( lapack_work( 1))
    call zhesv( 'U', ngqwf, wf_nst, apwolp, ngqwf, lapack_ipiv, rhs( 1:ngqwf, :), ngqwf, lapack_work, -1, lapack_info)
    lapack_lwork = lapack_work( 1)
    if( allocated( lapack_work)) deallocate( lapack_work)
    allocate( lapack_work( lapack_lwork))
    call zhesv( 'U', ngqwf, wf_nst, apwolp, ngqwf, lapack_ipiv, rhs( 1:ngqwf, :), ngqwf, lapack_work, lapack_lwork, lapack_info)
    if( lapack_info .ne. 0) then
      write( *, '(" ERROR (wannier_interpolate_eigsys): Lapack-routine ZHESV returned info ",I," for q-point ",I)') lapack_info, iq
    end if
    ! assign solutions to interpolated eigenvector
    evecq( 1:ngqwf, :) = rhs( 1:ngqwf, :)

    !write( *, '(2F13.6)') dot_product( conjg( evecq( :, wf_fst)), evecq( :, wf_fst))
    !write( *, '(2F13.6)') dot_product( conjg( evecq( :, wf_fst)), evecq( :, wf_lst))
    if( ongrid) then
      write(*,*) "interpolated"
      call plotmat( evecq, .true.)
      write(*,*)
    end if
    evecout( :, :, iq) = evecq

  end do

  deallocate( evecint, nlmo, lmo2l, lmo2m, lmo2o, evecfv, igxigwf, vgxlwf, gxcwf, tpgxcwf, vgkcwf, vgqcwf, sfacwf, uv, cuv, rhs, matchq, matchk, apwalm, evecq, lapack_work, lapack_ipiv)
  ngkmax = ngkmaxtmp

  return
end subroutine wannier_interpolate_eigsys

end module m_wannier_interpolate_eigsys
