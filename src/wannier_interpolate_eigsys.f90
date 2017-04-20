module m_wannier_interpolate_eigsys
    implicit none
    contains

subroutine wannier_interpolate_eigsys( evalin, int_kset, int_Gkset, evalout, evecout)
  use mod_wannier
  use m_wsweight
  use m_plotmat
  use mod_kqpts
  use mod_lattice
  use mod_eigenvalue_occupancy

  implicit none
  real(8), intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
  type( k_set), intent( in) :: int_kset
  type( Gk_set), intent( in) :: int_Gkset
  real(8), intent( out) :: evalout( wf_fst:wf_lst, int_kset%nkpt)
  complex(8), intent( out) :: evecout( int_Gkset%ngkmax+nlotot, wf_fst:wf_lst, int_kset%nkpt)
  
  integer :: nrpt, ix, iy, iz, ik, iq, ir
  integer :: ngkmaxwf, nmatmaxwf, ngkmaxint, nmatmaxint, ngkmax0, nmatmax0 
  complex(8) :: ftweight
  real(8) :: vqcin(3), v(3), vi(3), ffac, t1

  real(8), allocatable :: rptl(:,:), vkctmp(:,:)
  complex(8), allocatable :: auxmat(:,:), ueu(:,:,:), phase(:,:), hamilton(:,:,:), evecint(:,:,:), evecfv(:,:,:), apwolp(:,:), uv(:,:), cuv(:,:), rhs(:,:)
  complex(8), allocatable :: evecq(:,:)
  
  integer :: o, is, ia, ias, l, m, lm, lmo, nlmomax, igq, igk, ngktmp, ngqtmp
  integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:), lm2l(:)
  complex(8), allocatable :: matchq(:,:,:), matchk(:,:,:), apwalm(:,:,:,:,:)
  complex(8), allocatable :: p1(:,:), p2(:,:)

  integer :: lapack_lwork, lapack_info
  integer, allocatable :: lapack_ipiv(:)
  complex(8), allocatable :: lapack_work(:)

  logical :: ongrid
  character(256) :: fname

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
  allocate( phase( wf_kset%nkpt, int_kset%nkpt))
  allocate( hamilton( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iy, auxmat)
#endif
  allocate( auxmat( wf_fst:wf_lst, wf_fst:wf_lst))
#ifdef USEOMP
!$OMP DO
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
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( auxmat)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

  allocate( p1( wf_kset%nkpt, nrpt), p2( int_kset%nkpt, nrpt))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iq, ir)
!$OMP DO
#endif
  do ir = 1, nrpt
    do iq = 1, int_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), int_kset%vkl( :, iq), p2( iq, ir), kgrid=.true.)
    end do
    do ik = 1, wf_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), p1( ik, ir), kgrid=.true.)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  call zgemm( 'N', 'C', wf_kset%nkpt, int_kset%nkpt, nrpt, zone, &
       p1, wf_kset%nkpt, &
       p2, int_kset%nkpt, zzero, &
       phase, wf_kset%nkpt)

  deallocate( p1, p2)

  phase = phase/wf_kset%nkpt
  do iy = wf_fst, wf_lst
    call zgemm( 'N', 'N', wf_nst, int_kset%nkpt, wf_kset%nkpt, zone, &
         ueu( iy, :, :), wf_nst, &
         phase, wf_kset%nkpt, zzero, &
         hamilton( iy, :, :), wf_nst)
  end do
  deallocate( ueu)

  ! interpolation
  allocate( evecint( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq)
!$OMP DO
#endif
  do iq = 1, int_kset%nkpt 
    call diaghermat( wf_nst, hamilton( :, :, iq), evalout( :, iq), evecint( :, :, iq))
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  deallocate( rptl, hamilton)
!  allocate( evecint( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))
!  allocate( phase( wf_kset%nkpt, int_kset%nkpt))
!  evecint = zone
!  evalout = zzero
!  phase = zone
  
  !**********************************************
  ! interpolated eigenvectors in LAPW+lo basis
  !**********************************************

  call readfermi
  call linengy
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
  
  ngkmaxwf = wf_Gkset%ngkmax
  nmatmaxwf = ngkmaxwf+nlotot
  ngkmaxint = int_Gkset%ngkmax
  nmatmaxint = ngkmaxint+nlotot
  ngkmax0 = ngkmax
  nmatmax0 = nmatmax
  ngkmax = ngkmaxwf
  nmatmax = nmatmaxwf
  !nmatmax = max( nmatmaxwf, nmatmaxint)

  !allocate( evecfv( nmatmaxwf, nstsv, nspinor))
  !allocate( uv( wf_fst:wf_lst, wf_fst:wf_lst))
  !allocate( cuv( nmatmaxwf, wf_fst:wf_lst))
  allocate( matchq( nlmomax, ngkmaxint, natmtot))
  !allocate( matchk( nlmomax, ngkmaxwf, natmtot))
  !allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
  allocate( evecq( nmatmaxint, wf_fst:wf_lst))
  !allocate( apwolp( ngkmax, ngkmax))
  allocate( rhs( ngkmaxint, wf_fst:wf_lst))

  write(*,*) int_Gkset%ngk( 1, :)
  do iq = 1, int_kset%nkpt
    write(*, '(I)') iq 
    ngkmax = ngkmaxint
    ! get matching coefficients for G+q        
    allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
    write(*,*) int_kset%nkpt, int_Gkset%ngk( 1, iq)
    ngqtmp = int_Gkset%ngk( 1, iq)
    apwalm = zzero
    call match( ngqtmp, int_Gkset%gkc( :, 1, iq), int_Gkset%tpgkc( :, :, 1, iq), int_Gkset%sfacgk( :, :, 1, iq), apwalm( :, :, :, :, 1))
    matchq = zzero
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        do lmo = 1, nlmo( is)
          l = lmo2l( lmo, is)
          m = lmo2m( lmo, is)
          o = lmo2o( lmo, is)
          lm = idxlm( l, m)
          matchq( lmo, 1:ngqtmp, ias) = apwalm( 1:ngqtmp, o, lm, ias, 1)
        end do
      end do
    end do
    deallocate( apwalm)

    rhs = zzero
    evecq = zzero
    ngkmax = ngkmaxwf
    ongrid = .false.

    write(*, '(3F13.6)') int_kset%vkl( :, iq) 
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, evecfv, uv, cuv, ngktmp, apwalm, matchk, is, ia, ias, lmo, l, m, o, lm, apwolp, igk, igq, v, t1, ffac), reduction(+:rhs, evecq)
#endif    
  allocate( evecfv( nmatmaxwf, nstsv, nspinor))
  allocate( uv( wf_fst:wf_lst, wf_fst:wf_lst))
  allocate( cuv( nmatmaxwf, wf_fst:wf_lst))
  allocate( matchk( nlmomax, ngkmaxwf, natmtot))
  allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
  allocate( apwolp( ngkmaxint, ngkmaxwf))
#ifdef USEOMP
!$OMP DO
#endif    
    do ik = 1, wf_kset%nkpt
      !write( *, '(I3.3,2x)') ik

      ! read eigenvector      
      if( input%properties%wannier%input .eq. "groundstate") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmaxwf, nstfv, nspinor, evecfv)
      else
        call terminate
      end if
    
      call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
           wf_transform( :, :, ik), wf_nst, &
           evecint( :, :, iq), wf_nst, zzero, &
           uv, wf_nst)
      
      call zgemm( 'N', 'N', nmatmaxwf, wf_nst, wf_nst, zone, &
           evecfv( :, wf_fst:wf_lst, 1), nmatmaxwf, &
           uv, wf_nst, zzero, &
           cuv, nmatmaxwf)

      !----------!
      ! APW part !
      !----------!
      ! get matching coefficients for G+k        
      ngktmp = wf_Gkset%ngk( 1, ik)
      apwalm = zzero
      call match( ngktmp, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
      matchk = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do lmo = 1, nlmo( is)
            l = lmo2l( lmo, is)
            m = lmo2m( lmo, is)
            o = lmo2o( lmo, is)
            lm = idxlm( l, m)
            matchk( lmo, 1:ngktmp, ias) = apwalm( 1:ngktmp, o, lm, ias, 1)
          end do
        end do
      end do

      !if( norm2( wf_kset%vkl( :, ik) - int_kset%vkl( :, iq)) .lt. input%structure%epslat) then
      !  write( *, '(3F13.6)') int_kset%vkl( :, iq)
      !  write( *, '(3F13.6)') wf_kset%vkl( :, ik)
      !  call plotmat( uv)
      !  write( *, '(2F13.6)') phase( ik, iq)
      !  ongrid = .true.
      !  write(*,*) "original"
      !  call plotmat( evecfv( :, wf_fst:wf_lst, 1), .true.)
      !end if

      ! build APW q-k overlap
        ! muffin tin region
      apwolp = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          call zgemm( 'C', 'N', ngqtmp, ngktmp, nlmo( is), zone, &
               matchq( 1:nlmo( is), 1:ngqtmp, ias), nlmo( is), &
               matchk( 1:nlmo( is), 1:ngktmp, ias), nlmo( is), zone, &
               apwolp( 1:ngqtmp, 1:ngktmp), ngqtmp)
        end do
      end do
        ! interstitial region
        ! contributes just if q = k
      if( norm2( wf_kset%vkl( :, ik) - int_kset%vkl( :, iq)) .lt. input%structure%epslat) then
        do igk = 1, ngktmp
          do igq = 1, ngqtmp
            v(:) = wf_Gkset%vgkc( :, igk, 1, ik) - int_Gkset%vgkc( :, igq, 1, iq)
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
      end if
      !write( fname, '("olpqk_",I3.3,"_",I3.3)') iq, ik
      !call writematlab( apwolp, fname)
      ! sum up right hand side of linear system of equations
      call zgemm( 'N', 'N', ngqtmp, wf_nst, ngktmp, phase( ik, iq), &
           apwolp( 1:ngqtmp, 1:ngktmp), ngqtmp, &
           cuv( 1:ngktmp, :), ngktmp, zone, &
           rhs( 1:ngqtmp, :), ngqtmp)

      !---------!
      ! LO part !
      !---------!
      !evecq( (ngqwf+1):(ngqwf+nlotot), :) = evecq( (ngqwf+1):(ngqwf+nlotot), :) + phase( ik, iq)*cuv( (ngqwf+1):(ngqwf+nlotot), :)
      evecq( (ngqtmp+1):(ngqtmp+nlotot), :) = evecq( (ngqtmp+1):(ngqtmp+nlotot), :) + phase( ik, iq)*cuv( (ngktmp+1):(ngktmp+nlotot), :)
      
    end do
#ifdef USEOMP
!$OMP END DO
#endif    
  deallocate( evecfv, uv, cuv, matchk, apwalm, apwolp)
#ifdef USEOMP
!$OMP END PARALLEL
#endif    
    write(*,*)
    write(*, '("rhs built")')
    ! build APW q-q overlap
    allocate( apwolp( ngkmaxint, ngkmaxint))
      ! muffin tin region
    apwolp = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, ia, ias), reduction(+:apwolp)
!$OMP DO
#endif    
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        call zgemm( 'C', 'N', ngqtmp, ngqtmp, nlmo( is), zone, &
             matchq( 1:nlmo( is), 1:ngqtmp, ias), nlmo( is), &
             matchq( 1:nlmo( is), 1:ngqtmp, ias), nlmo( is), zone, &
             apwolp( 1:ngqtmp, 1:ngqtmp), ngqtmp)
      end do
    end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
      ! interstitial region
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( igk, igq, v, is, t1, ffac, ia), reduction(+:apwolp)
!$OMP DO
#endif    
    do igk = 1, ngqtmp
      do igq = 1, ngqtmp
        v(:) = int_Gkset%vgkc( :, igk, 1, iq) - int_Gkset%vgkc( :, igq, 1, iq)
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
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
    write(*, '("olp built")')

    !write( fname, '("olpqq_",I3.3)') iq
    !call writematlab( apwolp, fname)
    !write( fname, '("rhs_",I3.3)') iq
    !call writematlab( rhs, fname)

    ! solve linear system of equations
    if( allocated( lapack_ipiv)) deallocate( lapack_ipiv)
    allocate( lapack_ipiv( int_Gkset%ngk( 1, iq)))
    if( allocated( lapack_work)) deallocate( lapack_work)
    allocate( lapack_work( 1))
    call zhesv( 'U', ngqtmp, wf_nst, apwolp( 1:ngqtmp, 1:ngqtmp), ngqtmp, lapack_ipiv, rhs( 1:ngqtmp, :), ngqtmp, lapack_work, -1, lapack_info)
    lapack_lwork = lapack_work( 1)
    if( allocated( lapack_work)) deallocate( lapack_work)
    allocate( lapack_work( lapack_lwork))
    call zhesv( 'U', ngqtmp, wf_nst, apwolp( 1:ngqtmp, 1:ngqtmp), ngqtmp, lapack_ipiv, rhs( 1:ngqtmp, :), ngqtmp, lapack_work, lapack_lwork, lapack_info)
    if( lapack_info .ne. 0) then
      write( *, '(" ERROR (wannier_interpolate_eigsys): Lapack-routine ZHESV returned info ",I," for q-point ",I)') lapack_info, iq
    end if
    ! assign solutions to interpolated eigenvector
    evecq( 1:ngqtmp, :) = rhs( 1:ngqtmp, :)
    write(*, '("system solved")')

    !write( *, '(2F13.6)') dot_product( conjg( evecq( :, wf_fst)), evecq( :, wf_fst))
    !write( *, '(2F13.6)') dot_product( conjg( evecq( :, wf_fst)), evecq( :, wf_lst))
    !if( ongrid) then
    !  write(*,*) "interpolated"
    !  call plotmat( evecq, .true.)
    !  write(*,*)
    !end if
    evecout( :, :, iq) = evecq
    !call plotmat( evecout( :, :, iq))
    write(*,*)
    deallocate( apwolp)

    !write(*,*) "-----------------"
  end do

  deallocate( evecint, nlmo, lmo2l, lmo2m, lmo2o, rhs, matchq, evecq, lapack_work, lapack_ipiv)

  ngkmax = ngkmax0
  nmatmax = nmatmax0

  return
end subroutine wannier_interpolate_eigsys

end module m_wannier_interpolate_eigsys
