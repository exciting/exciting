module m_wannier_interpolate_eigsys2
    implicit none
    contains

subroutine wannier_interpolate_eigsys2( evalin, int_kset, int_Gkset, evalout, evecout, bc)
  use modmain
  use mod_wannier
  use m_wsweight
  use m_plotmat
  use mod_kqpts
  use mod_lattice
  use mod_eigenvalue_occupancy
  use mod_eigensystem
  use mod_atoms
  use mod_muffin_tin
  use mod_lattice
  use mod_constants

  implicit none
  real(8), intent( in) :: evalin( wf_fst:wf_lst, wf_kset%nkpt)
  type( k_set), intent( in) :: int_kset
  type( Gk_set), intent( in) :: int_Gkset
  real(8), intent( out) :: evalout( wf_fst:wf_lst, int_kset%nkpt)
  complex(8), intent( out) :: evecout( int_Gkset%ngkmax+nlotot, wf_fst:wf_lst, int_kset%nkpt)
  real(8), intent( out) :: bc( 0:3, natmtot, wf_fst:wf_lst, int_kset%nkpt)
  
  integer :: nrpt, ix, iy, iz, ik, iq, ir, v(3)
  integer :: ngkmaxwf, nmatmaxwf, ngkmaxint, nmatmaxint, ngkmax0, nmatmax0 
  complex(8) :: ftweight, z1
  real(8) :: vqcin(3), vi(3), t
  real(8), allocatable :: ffacg(:)

  real(8), allocatable :: rptl(:,:), vkctmp(:,:), olp_ll(:,:)
  complex(8), allocatable :: auxmat(:,:), ueu(:,:,:), phase(:,:), hamilton(:,:,:), evecint(:,:,:), evecfv(:,:,:), uv(:,:), cuv(:,:), rhs(:,:), rhs_cpy(:,:), olp_alk(:,:), olp_alq(:,:), olp_tot(:,:)
  complex(8), allocatable :: cfuntmp(:)
  
  integer :: o, is, ia, ias, ilo1, ilo2, l, m, lm, lmo, nlmomax, igq, igk, ngktmp, ngqtmp
  integer, allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:), lm2l(:)
  complex(8), allocatable :: matchq(:,:,:), matchk(:,:,:), apwalm(:,:,:,:,:)
  complex(8), allocatable :: p1(:,:), p2(:,:)

  integer :: lapack_lwork, lapack_info
  integer, allocatable :: lapack_ipiv(:)
  complex(8), allocatable :: lapack_work(:)

  integer :: lmax, lmmax
  complex(8), allocatable :: dmat_tmp(:,:,:), dmat(:,:), evecsv(:,:), evecfv2(:,:,:), apwalm2(:,:,:,:,:)
  complex(8), allocatable :: dmatk(:,:,:,:,:,:)

  logical :: ongrid
  character(256) :: fname

  integer, allocatable :: equik(:)
  integer :: nsym

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
    auxmat = zzero
    do iy = wf_fst, wf_lst
      auxmat( iy, :) = wf_transform( iy, :, ik)*evalin( iy, ik)
    end do
    call zgemm( 'C', 'N', wf_nst, wf_nst, wf_nst, zone, &
         wf_transform( :, :, ik), wf_nst, &
         auxmat, wf_nst, zzero, &
         ueu( :, :, ik), wf_nst)
  end do
#ifdef USEOMP
!$OMP END DO
#endif
  deallocate( auxmat)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

  allocate( p1( nrpt, wf_kset%nkpt), p2( nrpt, int_kset%nkpt))
  p2 = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, iq, ir)
!$OMP DO
#endif
  do ir = 1, nrpt
    do iq = 1, int_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), int_kset%vkl( :, iq), p2( ir, iq), kgrid=.true.)
    end do
    do ik = 1, wf_kset%nkpt
      call ws_weight( rptl( :, ir), rptl( :, ir), wf_kset%vkl( :, ik), p1( ir, ik), kgrid=.true.)
    end do
  end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  call zgemm( 'C', 'N', wf_kset%nkpt, int_kset%nkpt, nrpt, zone, &
       p1, nrpt, &
       p2, nrpt, zzero, &
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
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq)
!!$OMP DO
#endif
  do iq = 1, int_kset%nkpt 
    !write(*,'(I,3F13.6)') iq, int_kset%vkl( :, iq)
    !call plotmat( hamilton( :, :, iq), .true.)
    !write(*,*)
    call diaghermat( wf_nst, hamilton( :, :, iq), evalout( :, iq), evecint( :, :, iq))
    !call plotmat( evecint( :, :, iq), .true.)
    !write(*,*)
  end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
  deallocate( rptl, hamilton)
  evecout = zzero
  !write( fname, '("phase")')
  !call writematlab( phase, fname)
  
  !**********************************************
  ! interpolated bandcharacter  
  !**********************************************

  !evecint( :, 3, 2:(int_kset%nkpt-1)) = evecint( :, 3, 2:(int_kset%nkpt-1)) +   evecint( :, 4, 2:(int_kset%nkpt-1))
  !evecint( :, 4, 2:(int_kset%nkpt-1)) = evecint( :, 3, 2:(int_kset%nkpt-1)) - 2*evecint( :, 4, 2:(int_kset%nkpt-1))
  !evecint( :, 5, 2:(int_kset%nkpt-1)) = evecint( :, 5, 2:(int_kset%nkpt-1)) +   evecint( :, 6, 2:(int_kset%nkpt-1))
  !evecint( :, 6, 2:(int_kset%nkpt-1)) = evecint( :, 5, 2:(int_kset%nkpt-1)) - 2*evecint( :, 6, 2:(int_kset%nkpt-1))
  !evecint( :, 3:6, 2:(int_kset%nkpt-1)) = evecint( :, 3:6, 2:(int_kset%nkpt-1))/sqrt( 2.d0)

  allocate( equik( 100))
  do iq = 1, int_kset%nkpt
    call myfindsym( int_kset%vkl( :, iq), nsym, equik)
    write(*,*), iq, nsym!, equik( 1:nsym)
  end do

  do is = 1, nspecies
    do ia = 1, natoms( is)
      ias = idxas( ia, is)
      call wannier_interpolate_bandchar( 3, is, ia, int_kset, evecint, phase, bc( :, ias, :, :))
    end do
  end do

  !do iq = 1, int_kset%nkpt
  !  write(*,'(I,3F13.6)') iq, int_kset%vkl( :, iq)
  !  call plotmat( evecint( :, :, iq), .true.)
  !  write(*,*)
  !  do ik = 1, wf_kset%nkpt
  !    write(*,'(I,3F13.6)') ik, wf_kset%vkl( :, ik)
  !    call plotmat( matmul( wf_transform( :, :, ik), evecint( :, :, iq)), .true.)
  !    write(*,*)
  !  end do
  !  write(*,*)
  !end do
  return
  
  !**********************************************
  ! interpolated eigenvectors in LAPW+lo basis
  !**********************************************

  call readstate
  call readfermi
  call linengy
  call genapwfr
  call genlofr
  call olprad

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

  allocate( matchq( nlmomax, ngkmaxint, natmtot))
  allocate( rhs( nmatmaxint, wf_fst:wf_lst))
  allocate( rhs_cpy( nmatmaxint, wf_fst:wf_lst))
  allocate( olp_alq( ngkmaxint, nlotot))
  allocate( cfuntmp( wf_Gset%ngvec))
  allocate( ffacg( wf_Gset%ngvec))

  ! build characteristic function
  cfuntmp(:) = zzero
  cfuntmp(1) = zone
  do is = 1, nspecies
    do ix = 1, wf_Gset%ngvec
      if( wf_Gset%gc( ix) .gt. input%structure%epslat) then
        t = wf_Gset%gc( ix)*rmt( is)
        ffacg( ix) = fourpi/omega*(sin( t) - t*cos( t))/(wf_Gset%gc( ix)**3)
      else
        ffacg( ix) = fourpi/omega/3.d0*rmt( is)**3
      end if
    end do
    do ia = 1, natoms( is)
      do ix = 1, wf_Gset%ngvec
        t = -dot_product( wf_Gset%vgc( :, ix), atposc( :, ia, is))
        z1 = cmplx( cos( t), sin( t), 8)
        cfuntmp( ix) = cfuntmp( ix) - z1*ffacg( ix)
      end do
    end do
  end do
  deallocate( ffacg)

  ! build k- and q-independent LO-LO overlap and APW-LO radial integrals
  allocate( olp_ll( nlotot, nlotot))
  olp_ll = 0.d0
  do is = 1, nspecies
    do ia = 1, natoms( is)
      ias = idxas( ia, is)
      do ilo1 = 1, nlorb( is)
        l = lorbl( ilo1, is)
        ! LO-LO overlap
        do ilo2 = 1, nlorb( is)
          if( l .eq. lorbl( ilo2, is)) then
            lm = idxlm( l, -l)
            do lmo = idxlm( l, -l), idxlm( l, l)
              ix = idxlo( lm, ilo1, ias)+lmo-lm
              iy = idxlo( lm, ilo2, ias)+lmo-lm
              olp_ll( ix, iy) = olp_ll( ix, iy) + ololo( ilo1, ilo2, ias)
            end do
          end if
        end do
      end do
    end do
  end do
  write( fname, '("olp_ll")')
  !call writematlab( zone*olp_ll( :, :), fname)

  write( fname, '("phase")')
  !call writematlab( phase, fname)

  write(*,*) int_Gkset%ngk( 1, :)
  do iq = 1, int_kset%nkpt
    write(*, '(I)') iq 
    ngkmax = ngkmaxint
    ! get matching coefficients for G+q        
    allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
    !write(*,*) int_kset%nkpt, int_Gkset%ngk( 1, iq)
    ngqtmp = int_Gkset%ngk( 1, iq)
    apwalm = zzero
    call match( ngqtmp, int_Gkset%gkc( :, 1, iq), int_Gkset%tpgkc( :, :, 1, iq), int_Gkset%sfacgk( :, :, 1, iq), apwalm( :, :, :, :, 1))
    matchq = zzero
    olp_alq = zzero
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        ! build q-dependent APW-LO overlap
        do ilo1 = 1, nlorb( is)
          l = lorbl( ilo1, is)
          lm = idxlm( l, -l)
          lmo = idxlm( l, l)
          ix = idxlo( lm, ilo1, ias)
          iy = idxlo( lmo, ilo1, ias)
          do o = 1, apword( l, is)
            olp_alq( 1:ngqtmp, ix:iy) = olp_alq( 1:ngktmp, ix:iy) + conjg( apwalm( 1:ngqtmp, o, lm:lmo, ias, 1)*oalo( o, ilo1,              ias))
          end do
        end do
        do lmo = 1, nlmo( is)
          l = lmo2l( lmo, is)
          m = lmo2m( lmo, is)
          o = lmo2o( lmo, is)
          lm = idxlm( l, m)
          matchq( lmo, 1:ngqtmp, ias) = apwalm( 1:ngqtmp, o, lm, ias, 1)
        end do
        write( fname, '("matchq/matchq",2("_",I3.3))') iq, ias
        !call writematlab( matchq( :, 1:ngqtmp, ias), fname)
      end do
    end do

    deallocate( apwalm)
    write( fname, '("olp_alq/olp_alq_",I3.3)') iq
    !call writematlab( olp_alq( 1:ngqtmp, :), fname)

    write( fname, '("evecint/evecint_",I3.3)') iq
    !call writematlab( evecint( :, :, iq), fname)
    
    rhs = zzero
    ngkmax = ngkmaxwf
    ongrid = .false.

    !write(*, '(3F13.6)') int_kset%vkl( :, iq) 
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ik, evecfv, uv, cuv, ngktmp, apwalm, matchk, ix, iy, is, ia, ias, lmo, l, m, o, lm, olp_tot, igk, igq, v, olp_alk), reduction(+:rhs)
#endif    
    allocate( evecfv( nmatmaxwf, nstsv, nspinor))
    allocate( uv( wf_fst:wf_lst, wf_fst:wf_lst))
    allocate( cuv( nmatmaxwf, wf_fst:wf_lst))
    allocate( matchk( nlmomax, ngkmaxwf, natmtot))
    allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspinor))
    allocate( olp_tot( nmatmaxint, nmatmaxwf))
    allocate( olp_alk( ngkmaxwf, nlotot))
#ifdef USEOMP
!$OMP DO
#endif    
    do ik = 1, wf_kset%nkpt
      ngktmp = wf_Gkset%ngk( 1, ik)
      cuv = zzero
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
    
      write( fname, '("evec0/evec0_",I3.3)') ik
      !call writematlab( evecfv( 1:(ngktmp+nlotot), wf_fst:wf_lst, 1), fname)
    
      write( fname, '("transform/transform_",I3.3)') ik
      !call writematlab( wf_transform( :, :, ik), fname)
    
      call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, zone, &
           wf_transform( :, :, ik), wf_nst, &
           evecint( :, :, iq), wf_nst, zzero, &
           uv, wf_nst)

      !if( norm2( int_kset%vkl( :, iq) - wf_kset%vkl( :, ik)) .lt. input%structure%epslat) then
      !  call plotmat( uv)
      !end if
      
      call zgemm( 'N', 'N', ngktmp+nlotot, wf_nst, wf_nst, zone, &
           evecfv( 1:(ngktmp+nlotot), wf_fst:wf_lst, 1), ngktmp+nlotot, &
           uv, wf_nst, zzero, &
           cuv( 1:(ngktmp+nlotot), :), ngktmp+nlotot)

      !----------!
      ! APW part !
      !----------!
      ! get matching coefficients for G+k        
      apwalm = zzero
      call match( ngktmp, wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
      matchk = zzero
      olp_alk = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          ! build k-dependent APW-LO overlap
          do ilo1 = 1, nlorb( is)
            l = lorbl( ilo1, is)
            lm = idxlm( l, -l)
            lmo = idxlm( l, l)
            ix = idxlo( lm, ilo1, ias)
            iy = idxlo( lmo, ilo1, ias)
            do o = 1, apword( l, is)
              olp_alk( 1:ngktmp, ix:iy) = olp_alk( 1:ngktmp, ix:iy) + conjg( apwalm( 1:ngktmp, o, lm:lmo, ias, 1)*oalo( o, ilo1,              ias))
            end do
          end do
          do lmo = 1, nlmo( is)
            l = lmo2l( lmo, is)
            m = lmo2m( lmo, is)
            o = lmo2o( lmo, is)
            lm = idxlm( l, m)
            matchk( lmo, 1:ngktmp, ias) = apwalm( 1:ngktmp, o, lm, ias, 1)
          end do
          write( fname, '("matchk/matchk",2("_",I3.3))') ik, ias
          !call writematlab( matchk( :, 1:ngktmp, ias), fname)
        end do
      end do

      write( fname, '("olp_alk/olp_alk_",I3.3))') ik
      !call writematlab( olp_alk( 1:ngktmp, :), fname)

      ! build APW-APW q-k overlap
        ! muffin tin region
      olp_tot = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          call zgemm( 'C', 'N', ngqtmp, ngktmp, nlmo( is), zone, &
               matchq( 1:nlmo( is), 1:ngqtmp, ias), nlmo( is), &
               matchk( 1:nlmo( is), 1:ngktmp, ias), nlmo( is), zone, &
               olp_tot( 1:ngqtmp, 1:ngktmp), ngqtmp)
        end do
      end do
        ! interstitial region
        ! contributes just if q = k
      if( norm2( wf_kset%vkl( :, ik) - int_kset%vkl( :, iq)) .lt. input%structure%epslat) then
        do igk = 1, ngktmp
          do igq = 1, ngqtmp
            v(:) = wf_Gset%ivg( :, int_Gkset%igkig( igq, 1, iq)) - ivg( :, wf_Gkset%igkig( igk, 1, ik))
            ix = ivgig( v(1), v(2), v(3))
            if( (ix .gt. 0) .and. (ix .le. ngvec)) then
              olp_tot( igq, igk) = olp_tot( igq, igk) + cfunig( ix)
            end if
          end do
        end do
      end if

      ! assign rest of total overlap
      olp_tot( (ngqtmp+1):(ngqtmp+nlotot), (ngktmp+1):(ngktmp+nlotot)) = zone*olp_ll(:,:)
      olp_tot( 1:ngqtmp, (ngktmp+1):(ngktmp+nlotot)) = olp_alq( 1:ngqtmp, :)
      olp_tot( (ngqtmp+1):(ngqtmp+nlotot), 1:ngktmp) = conjg( transpose( olp_alk( 1:ngktmp, :)))
      write( fname, '("olpqk/olpqk",2("_",I3.3))') iq, ik
      !call writematlab( olp_tot( 1:(ngqtmp+nlotot), 1:(ngktmp+nlotot)), fname)

      ! sum up right hand side of linear system of equations
      call zgemm( 'N', 'N', ngqtmp+nlotot, wf_nst, ngktmp+nlotot, phase( ik, iq), &
           olp_tot( 1:(ngqtmp+nlotot), 1:(ngktmp+nlotot)), ngqtmp+nlotot, &
           cuv( 1:(ngktmp+nlotot), :), ngktmp+nlotot, zone, &
           rhs( 1:(ngqtmp+nlotot), :), ngqtmp+nlotot)

    end do
#ifdef USEOMP
!$OMP END DO
#endif    
    deallocate( evecfv, uv, cuv, matchk, apwalm, olp_tot, olp_alk)
#ifdef USEOMP
!$OMP END PARALLEL
#endif    
    !write(*,*)
    !write(*, '("rhs built")')
    ! build APW-APW q-q overlap
    rhs_cpy = rhs
    allocate( olp_tot( nmatmaxint, nmatmaxint))
      ! muffin tin region
    olp_tot = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, ia, ias), reduction(+:olp_tot)
!$OMP DO
#endif    
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        call zgemm( 'C', 'N', ngqtmp, ngqtmp, nlmo( is), zone, &
             matchq( 1:nlmo( is), 1:ngqtmp, ias), nlmo( is), &
             matchq( 1:nlmo( is), 1:ngqtmp, ias), nlmo( is), zone, &
             olp_tot( 1:ngqtmp, 1:ngqtmp), ngqtmp)
      end do
    end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
      ! interstitial region
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( igk, igq, v, ix, is, ia), reduction(+:olp_tot)
!$OMP DO
#endif    
    do igk = 1, ngqtmp
      do igq = 1, ngqtmp
        v(:) = wf_Gset%ivg( :, int_Gkset%igkig( igq, 1, iq)) - wf_Gset%ivg( :, int_Gkset%igkig( igk, 1, iq))
        ix = ivgig( v(1), v(2), v(3))
        if( (ix .gt. 0) .and. (ix .le. ngvec)) then
          olp_tot( igq, igk) = olp_tot( igq, igk) + cfunig( ix)
        end if
      end do
    end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif    
    
    ! assign rest of total overlap
    olp_tot( (ngqtmp+1):(ngqtmp+nlotot), (ngqtmp+1):(ngqtmp+nlotot)) = zone*olp_ll(:,:)
    olp_tot( 1:ngqtmp, (ngqtmp+1):(ngqtmp+nlotot)) = olp_alq( 1:ngqtmp, :)
    olp_tot( (ngqtmp+1):(ngqtmp+nlotot), 1:ngqtmp) = conjg( transpose( olp_alq( 1:ngqtmp, :)))
    !write(*, '("olp built")')

    write( fname, '("olpqq/olpqq_",I3.3)') iq
    call writematlab( olp_tot( 1:(ngqtmp+nlotot), 1:(ngqtmp+nlotot)), fname)
    write( fname, '("rhs/rhs_",I3.3)') iq
    !call writematlab( rhs( 1:(ngqtmp+nlotot), :), fname)

    ! solve linear system of equations
    if( allocated( lapack_ipiv)) deallocate( lapack_ipiv)
    allocate( lapack_ipiv( ngqtmp+nlotot))
    if( allocated( lapack_work)) deallocate( lapack_work)
    allocate( lapack_work( 1))
    call zhesv( 'U', ngqtmp+nlotot, wf_nst, &
         olp_tot( 1:(ngqtmp+nlotot), 1:(ngqtmp+nlotot)), ngqtmp+nlotot, lapack_ipiv, &
         rhs( 1:(ngqtmp+nlotot), :), ngqtmp+nlotot, lapack_work, -1, lapack_info)
    lapack_lwork = lapack_work( 1)
    if( allocated( lapack_work)) deallocate( lapack_work)
    allocate( lapack_work( lapack_lwork))
    call zhesv( 'U', ngqtmp+nlotot, wf_nst, &
         olp_tot( 1:(ngqtmp+nlotot), 1:(ngqtmp+nlotot)), ngqtmp+nlotot, lapack_ipiv, &
         rhs( 1:(ngqtmp+nlotot), :), ngqtmp+nlotot, lapack_work, lapack_lwork, lapack_info)
    if( lapack_info .ne. 0) then
      write( *, '(" ERROR (wannier_interpolate_eigsys): Lapack-routine ZHESV returned info ",I," for q-point ",I)') lapack_info, iq
    end if
    ! assign solutions to interpolated eigenvector
    !do ix = wf_fst, wf_lst
    !  z1 = dot_product( rhs( 1:(ngqtmp+nlotot), ix), rhs_cpy( 1:(ngqtmp+nlotot), ix))
    !  write(*,'(2(2F13.6,3x))') z1, zone/sqrt( z1)
    !  evecout( 1:(ngqtmp+nlotot), ix, iq) = rhs( 1:(ngqtmp+nlotot), ix)/sqrt( z1)
    !  rhs( 1:(ngqtmp+nlotot), ix) = rhs_cpy( 1:(ngqtmp+nlotot), ix)/sqrt( z1)
    !end do
    !write(*, '("system solved")')
    evecout( 1:(ngqtmp+nlotot), :, iq) = rhs( 1:(ngqtmp+nlotot), :)

    deallocate( olp_tot)

    write( fname, '("evec/evec_",I3.3)') iq
    !call writematlab( evecout( 1:(ngqtmp+nlotot), :, iq), fname)

    !allocate( auxmat( wf_nst, wf_nst))
    !call zgemm( 'C', 'N', wf_nst, wf_nst, ngqtmp+nlotot, zone, &
    !     evecout( 1:(ngqtmp+nlotot), :, iq), ngqtmp+nlotot, &
    !     rhs( 1:(ngqtmp+nlotot), :), ngqtmp+nlotot, zzero, &
    !     auxmat, wf_nst)
    !write( fname, '("check/check_",I3.3)') iq
    !call writematlab( auxmat, fname)
    !deallocate( auxmat)
  end do

  deallocate( evecint, nlmo, lmo2l, lmo2m, lmo2o, rhs, matchq, lapack_work, lapack_ipiv)

  ngkmax = ngkmax0
  nmatmax = nmatmax0

  return
end subroutine wannier_interpolate_eigsys2

end module m_wannier_interpolate_eigsys2
