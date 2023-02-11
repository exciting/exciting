module mod_pwmat
  use modinput
  use modmpi
  use mod_atoms
  use mod_APW_LO
  use mod_muffin_tin
  use gaunt
  use constants, only: zzero, zone, zi, zil, fourpi
  use mod_eigensystem, only : idxlo, nmatmax_ptr
  use mod_lattice, only : bvec
  use mod_Gvector, only : ngvec, gc, cfunig, ivg, ivgig, ngrtot, ngrid, intgv
  use mod_Gkvector, only : gkmax, ngkmax_ptr
  use mod_eigenvalue_occupancy, only : nstfv
  use mod_misc, only : filext
  use mod_kpointset
  use modxs, only : fftmap_type

  implicit none
  private

  integer              :: pwmat_fst1, pwmat_lst1, pwmat_fst2, pwmat_lst2, pwmat_nst1, pwmat_nst2
  integer              :: ng, nlammax, nlmlammax, pwmat_lmaxapw, pwmat_lmaxexp
  real(8)              :: pwmat_gmax
  real(8)              :: vecql(3)
  character(256)       :: pwmat_fname
  logical              :: pwmat_usefft
  type( k_set)         :: pwmat_kset
  type( fftmap_type)   :: pwmat_fftmap

  integer, allocatable    :: nlam(:,:), nlmlam(:), lam2apwlo(:,:,:), idxlmlam(:,:,:), lm2l(:)
  integer, allocatable    :: vecgl(:,:)
  complex(8), allocatable :: rignt(:,:,:,:), pwmat_cfunir(:)

  real(8) :: ts, te, t0, t1
  real(8) :: ttot, tinit, tri, trignt, tgk, tprep, tmt, tgg, tir, tio

  public :: pwmat_init, pwmat_prepare, pwmat_init_qg, pwmat_genpwmat, pwmat_destroy

! variable

! methods
  contains
      subroutine pwmat_init( apwlmax, explmax, kset, fst1, lst1, fst2, lst2, fname, gmax, fft)
          integer, intent( in) :: apwlmax, explmax, fst1, lst1, fst2, lst2
          type( k_set), intent( in) :: kset
          character(*), optional, intent( in) :: fname
          real(8), optional, intent( in) :: gmax
          logical, optional, intent( in) :: fft

          integer :: l1, m1, o1, lm1, is, ilo, lam

          ttot   = 0.d0
          tinit  = 0.d0
          tri    = 0.d0
          trignt = 0.d0
          tgk    = 0.d0
          tprep  = 0.d0
          tmt    = 0.d0
          tgg    = 0.d0
          tir    = 0.d0
          tio    = 0.d0

          call timesec( ts)

          pwmat_fname = 'EVEC'
          if( present( fname)) pwmat_fname = trim( fname)

          pwmat_gmax = input%groundstate%gmaxvr
          if( present( gmax)) pwmat_gmax = gmax

          pwmat_usefft = .false.
          if( present( fft)) pwmat_usefft = fft

          pwmat_lmaxapw = min( apwlmax, input%groundstate%lmaxapw)
          pwmat_lmaxexp = min( explmax, pwmat_lmaxapw)
          pwmat_kset = kset
          pwmat_fst1 = fst1
          pwmat_lst1 = lst1
          pwmat_fst2 = fst2
          pwmat_lst2 = lst2
          pwmat_nst1 = pwmat_lst1-pwmat_fst1+1
          pwmat_nst2 = pwmat_lst2-pwmat_fst2+1

          ! count combined (l,m,o) indices and build index maps
          allocate( nlam( 0:pwmat_lmaxapw, nspecies))
          allocate( lam2apwlo( apwordmax+nlomax, 0:pwmat_lmaxapw, nspecies))
          nlammax = 0
          do is = 1, nspecies
            do l1 = 0, pwmat_lmaxapw
              nlam( l1, is) = 0
              lam2apwlo( :, l1, is) = 0
              do o1 = 1, apword( l1, is)
                nlam( l1, is) = nlam( l1, is)+1
                lam2apwlo( nlam( l1, is), l1, is) = o1
              end do
              do ilo = 1, nlorb( is)
                if( lorbl( ilo, is) .eq. l1) then
                  nlam( l1, is) = nlam( l1, is) + 1
                  lam2apwlo( nlam( l1, is), l1, is) = -ilo
                end if
              end do
              nlammax = max( nlammax, nlam( l1, is))
            end do
          end do

          allocate( lm2l( (pwmat_lmaxapw + 1)**2))
          allocate( nlmlam( nspecies))
          allocate( idxlmlam( (pwmat_lmaxapw + 1)**2, nlammax, nspecies)) 
          idxlmlam(:,:,:) = 0
          nlmlammax = 0
          nlmlam(:) = 0
          do l1 = 0, pwmat_lmaxapw
            do m1 = -l1, l1
              lm1 = idxlm( l1, m1)
              lm2l( lm1) = l1
              do is = 1, nspecies
                do lam = 1, nlam( l1, is)
                  nlmlam( is) = nlmlam( is) + 1
                  idxlmlam( lm1, lam, is) = nlmlam( is)
                end do
              end do
            end do
          end do
          nlmlammax = maxval( nlmlam, 1)
          !write(*,'("pwmat init: nlmlammax = :",I4)') nlmlammax

          ! check if Gaunt coefficients are available and create them if not
          if( .not. gaunt_coeff_yyy%check_bounds( pwmat_lmaxapw, pwmat_lmaxexp, pwmat_lmaxapw)) &
            gaunt_coeff_yyy = non_zero_gaunt_yyy( pwmat_lmaxapw, pwmat_lmaxexp, pwmat_lmaxapw)

          call timesec( t1)
          tinit = tinit + t1 - ts

          return
      end subroutine pwmat_init

      subroutine pwmat_prepare( ik, evec)
          integer, intent( in)      :: ik
          complex(8), intent( in)   :: evec( nmatmax_ptr, nstfv)

          integer :: ngp, is, ia, ias, lm, l, lam, o, ilo, fst, lst, nst
          integer :: ig, igp

          integer, allocatable :: igpig(:) 
          real(8), allocatable :: vgpl(:,:), vgpc(:,:), gpc(:), tpgpc(:,:)
          complex(8), allocatable :: evecmt(:,:,:), evecir(:,:)
          complex(8), allocatable :: sfacgp(:,:), apwalm(:,:,:,:), auxvec(:)

          fst = min( pwmat_fst1, pwmat_fst2)
          lst = max( pwmat_lst1, pwmat_lst2)
          nst = lst-fst+1

          call timesec( t0)
          allocate( evecmt( nlmlammax, fst:lst, natmtot))
          allocate( igpig( ngkmax_ptr), vgpl( 3, ngkmax_ptr), vgpc( 3, ngkmax_ptr), gpc( ngkmax_ptr), tpgpc( 2, ngkmax_ptr))
          allocate( sfacgp( ngkmax_ptr, natmtot))
          allocate( apwalm( ngkmax_ptr, apwordmax, lmmaxapw, natmtot))
          allocate( auxvec( fst:lst))

          ! generate the G+k-vectors
          call gengpvec( pwmat_kset%vkl( :, ik), pwmat_kset%vkc( :, ik), ngp, igpig, vgpl, vgpc, gpc, tpgpc)
          ! generate the structure factors
          call gensfacgp( ngp, vgpc, ngkmax_ptr, sfacgp)
          ! find matching coefficients for k-point k+q
          call match( ngp, gpc, tpgpc, sfacgp, apwalm)

          evecmt = zzero
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              do lm = 1, (pwmat_lmaxapw + 1)**2
                l = lm2l( lm)
                do lam = 1, nlam( l, is)
                  ! lam belongs to apw
                  if( lam2apwlo( lam, l, is) .gt. 0) then
                    o = lam2apwlo( lam, l, is)
                    call zgemv( 't', ngp, nst, zone, &
                           evec( :, fst:lst), nmatmax_ptr, &
                           apwalm( :, o, lm, ias), 1, zzero, &
                           auxvec, 1)
                    evecmt( idxlmlam( lm, lam, is), :, ias) = auxvec
                  end if
                  ! lam belongs to lo
                  if( lam2apwlo( lam, l, is) .lt. 0) then
                    ilo = -lam2apwlo( lam, l, is)
                    evecmt( idxlmlam( lm, lam, is), :, ias) = evec( ngp+idxlo( lm, ilo, ias), fst:lst)
                  end if
                end do
              end do
            end do
          end do

          call timesec( t1)
          tprep = tprep + t1 - t0
          call pwmat_putevecmt( pwmat_kset%vkl( :, ik), fst, lst, evecmt)
          call timesec( t0)
          tio = tio + t0 - t1
          deallocate( evecmt, vgpl, vgpc, gpc, tpgpc, sfacgp, apwalm, auxvec)

          if( pwmat_usefft) then
            if( .not. allocated( pwmat_cfunir)) then
              call genfftmap( pwmat_fftmap, pwmat_gmax)
              write(*,*) pwmat_fftmap%ngrid, pwmat_fftmap%ngrtot
              write(*,*) ngrid, ngrtot
              allocate( pwmat_cfunir( pwmat_fftmap%ngrtot+1))
              pwmat_cfunir = zzero
              do ig = 1, ngvec
                if( gc( ig) .lt. pwmat_gmax) pwmat_cfunir( pwmat_fftmap%igfft( ig)) = cfunig( ig)
              end do
              call zfftifc( 3, pwmat_fftmap%ngrid, 1, pwmat_cfunir)
            end if

            allocate( evecir( pwmat_fftmap%ngrtot+1, fst:lst))
            evecir = zzero
            do is = fst, lst
              do igp = 1, ngp
                ig = igpig( igp)
                evecir( pwmat_fftmap%igfft( ig), is) = evec( igp, is)
              end do
              call zfftifc( 3, pwmat_fftmap%ngrid, 1, evecir( :, is))
            end do
            
            call timesec( t1)
            tprep = tprep + t1 - t0
            call pwmat_putevecir( pwmat_kset%vkl( :, ik), fst, lst, evecir)
            call timesec( t0)
            tio = tio + t0 - t1
            deallocate( evecir)
          end if
          deallocate( igpig)
          call timesec( t1)
          tprep = tprep + t1 - t0
          return
      end subroutine pwmat_prepare
      
      subroutine pwmat_init_qg( vecql_, vecgl_, ng_)
          real(8), intent( in) :: vecql_(3)
          integer, intent( in) :: vecgl_(3,*), ng_

          integer :: l1, l2, l3, m1, m2, lm1, lm2, lm3, is, ia, ias, i, ig, vg(3)
          integer :: lam1, lam2, idx1, idx2
          real(8) :: qgc, tp(2), vr(3)
          complex(8) :: strf, ylm( (pwmat_lmaxexp + 1)**2)

          real(8), allocatable :: radint(:,:,:)
          character(256) :: fname

          ng = ng_
          vecql = vecql_
          call r3frac( input%structure%epslat, vecql, vg)
          if( allocated( vecgl)) deallocate( vecgl)
          allocate( vecgl( 3, ng))
          do ig = 1, ng
            vecgl( :, ig) = vecgl_( :, ig) + vg
          end do
          !write(*,'("vq and vg: ",3f13.6,3i)') vecql, vecgl
          
          ! compute radial integrals times Gaunt and expansion prefactor
          call timesec( t0)
          if( allocated( rignt)) deallocate( rignt)
          allocate( rignt( nlmlammax, nlmlammax, natmtot, ng))
          rignt = zzero
          allocate( radint( nlammax, nlammax, 0:pwmat_lmaxexp))
          ! generate radial functions
          !call readstate
          !call readfermi
          !call linengy
          !call genapwfr
          !call genlofr
          ! multiply radial integral with Gaunt coefficients and expansion prefactor
          ! generate spherical harmonics

          do ig = 1, ng
            call r3mv( bvec, vecql + dble( vecgl( :, ig)), vr)
            call sphcrd( vr, qgc, tp)
            call genylm( pwmat_lmaxexp, tp, ylm)

            do is = 1, nspecies
              do ia = 1, natoms( is)
                ias = idxas( ia, is)
                strf = fourpi*exp( -zi*dot_product( vr, atposc( :, ia, is)))
                ! generate radial integral
                do l1 = 0, pwmat_lmaxapw
                  do l2 = 0, pwmat_lmaxapw
                    call timesec( t1)
                    trignt = trignt + t1 - t0
                    call pwmat_genri( is, ia, l1, l2, ig, radint)
                    call timesec( t0)
                    tri = tri + t0 - t1
                    !write(*,'("init_qg: ",3I)') ias, l1, l2
                    do m1 = -l1, l1
                      lm1 = idxlm( l1, m1)
                      do m2 = -l2, l2
                        lm2 = idxlm( l2, m2)

                        do lam1 = 1, nlam( l1, is)
                          idx1 = idxlmlam( lm1, lam1, is)
                          do lam2 = 1, nlam( l2, is)
                            idx2 = idxlmlam( lm2, lam2, is)
                          
                            do i = 1, gaunt_coeff_yyy%num(lm1, lm2)
                              lm3 = gaunt_coeff_yyy%lm2(i, lm1, lm2)
                              l3 = lm2l(lm3)
                              rignt(idx1, idx2, ias, ig) = rignt(idx1, idx2, ias, ig) + &
                                   conjg( zil(l3) ) * conjg( ylm(lm3) ) * cmplx( gaunt_coeff_yyy%val(i, lm1, lm2) * radint( lam1, lam2, l3), 0, 8 )
                            end do
                            
                          end do
                        end do

                      end do
                    end do
                  end do
                end do

                rignt( :, :, ias, ig) = strf*rignt( :, :, ias, ig)
              end do ! atoms
            end do !species

          end do

          deallocate( radint)
          call timesec( t1)
          trignt = trignt + t1 - t0
          return
      end subroutine pwmat_init_qg
      
      subroutine pwmat_genri( is, ia, l1, l2, ig, radint)
          integer, intent( in) :: is, ia, l1, l2, ig
          real(8), intent( out) :: radint( nlammax, nlammax, 0:pwmat_lmaxexp)

          integer :: ias, ir, nr, l3, o1, o2, ilo1, ilo2, lam1, lam2
          real(8) :: x, qgc, tp(2), vr(3)

          real(8), allocatable :: jlqgr(:,:), fr(:), gf(:), cf(:,:)
          
          radint(:,:,:) = 0.d0
          call r3mv( bvec, vecql + dble( vecgl( :, ig)), vr)
          
          ias = idxas( ia, is)
          call sphcrd( vr, qgc, tp)
          nr = nrmt( is)
          ! generate spherical Bessel functions
          allocate( jlqgr( 0:pwmat_lmaxexp, nr))
          do ir = 1, nr
            x = qgc*spr( ir, is)
            call sbessel( pwmat_lmaxexp, x, jlqgr( :, ir))
          end do

#ifdef USEOMP
!!$omp parallel default( shared) private( lam1, lam2, o1, o2, ilo1, ilo2, l3, ir, fr, gf, cf)
#endif
          allocate( fr( nr), gf( nr), cf( 3, nr))
#ifdef USEOMP
!!$omp do collapse( 2)
#endif
          do lam1 = 1, nlam( l1, is)
            do lam2 = 1, nlam( l2, is)
              
              ! lam1 belongs to apw
              if( lam2apwlo( lam1, l1, is) .gt. 0) then
                o1 = lam2apwlo( lam1, l1, is)
                ! lam2 belongs to apw
                if( lam2apwlo( lam2, l2, is) .gt. 0) then
                  o2 = lam2apwlo( lam2, l2, is)
                  do l3 = 0, pwmat_lmaxexp
                    do ir = 1, nr
                      fr( ir) = apwfr( ir, 1, o1, l1, ias)*jlqgr( l3, ir)*apwfr( ir, 1, o2, l2, ias)*spr( ir, is)*spr( ir, is)
                    end do
                    call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                    radint( lam1, lam2, l3) = gf( nr)
                  end do
                end if
                ! lam2 belongs to lo
                if( lam2apwlo( lam2, l2, is) .lt. 0) then
                  ilo2 = -lam2apwlo( lam2, l2, is)
                  do l3 = 0, pwmat_lmaxexp
                    do ir = 1, nr
                      fr( ir) = apwfr( ir, 1, o1, l1, ias)*jlqgr( l3, ir)*lofr( ir, 1, ilo2, ias)*spr( ir, is)*spr( ir, is)
                    end do
                    call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                    radint( lam1, lam2, l3) = gf( nr)
                  end do
                end if
              end if
              ! lam1 belongs to lo
              if( lam2apwlo( lam1, l1, is) .lt. 0) then
                ilo1 = -lam2apwlo( lam1, l1, is)
                ! lam2 belongs to apw
                if( lam2apwlo( lam2, l2, is) .gt. 0) then
                  o2 = lam2apwlo( lam2, l2, is)
                  do l3 = 0, pwmat_lmaxexp
                    do ir = 1, nr
                      fr( ir) = lofr( ir, 1, ilo1, ias)*jlqgr( l3, ir)*apwfr( ir, 1, o2, l2, ias)*spr( ir, is)*spr( ir, is)
                    end do
                    call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                    radint( lam1, lam2, l3) = gf( nr)
                  end do
                end if
                ! lam2 belongs to lo
                if( lam2apwlo( lam2, l2, is) .lt. 0) then
                  ilo2 = -lam2apwlo( lam2, l2, is)
                  do l3 = 0, pwmat_lmaxexp
                    do ir = 1, nr
                      fr( ir) = lofr( ir, 1, ilo1, ias)*jlqgr( l3, ir)*lofr( ir, 1, ilo2, ias)*spr( ir, is)*spr( ir, is)
                    end do
                    call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                    radint( lam1, lam2, l3) = gf( nr)
                  end do
                end if
              end if

            end do
          end do
#ifdef USEOMP
!!$omp end do
#endif
          deallocate( fr, gf, cf)
#ifdef USEOMP
!!$omp end parallel
#endif
          deallocate( jlqgr)

          return
      end subroutine pwmat_genri
      
      subroutine pwmat_genpwmat( ik, evec1, evec2, pwmat)
          integer, intent( in)     :: ik
          complex(8), intent( in)  :: evec1( nmatmax_ptr, *), evec2( nmatmax_ptr, *)
          complex(8), intent( out) :: pwmat( pwmat_fst1:pwmat_lst1, pwmat_fst2:pwmat_lst2, ng)

          integer :: ngp, ngpq, is, ia, ias 
          real(8) :: veckl(3), veckc(3), veckql(3), veckqc(3)
          integer :: shift(3), g(3), gs(3), igk, igq, ig

          integer, allocatable :: igpig(:), igpqig(:) 
          real(8), allocatable :: vgpql(:,:), vgpqc(:,:), gkqc(:), tpgkqc(:,:)
          complex(8), allocatable :: auxmat(:,:), evecmt1(:,:,:), evecmt2(:,:,:), cfunmat(:,:)

          pwmat = zzero
          
          veckl = pwmat_kset%vkl( :, ik)
          veckc = pwmat_kset%vkc( :, ik)
          ! k+q-vector in lattice coordinates
          veckql = veckl + vecql
          ! map vector components to [0,1) interval
          call r3frac( input%structure%epslat, veckql, shift)
          ! k+q-vector in Cartesian coordinates
          call r3mv( bvec, veckql, veckqc)
          
          !--------------------------------------!
          !      muffin-tin matrix elements      !
          !--------------------------------------!

          allocate( igpig( ngkmax_ptr), igpqig( ngkmax_ptr), vgpql( 3, ngkmax_ptr), vgpqc( 3, ngkmax_ptr), gkqc( ngkmax_ptr), tpgkqc( 2, ngkmax_ptr))

          call timesec( t0)
          ! generate the G+k+q-vectors
          call gengpvec( veckql, veckqc, ngpq, igpqig, vgpql, vgpqc, gkqc, tpgkqc)
          ! generate the G+k-vectors
          call gengpvec( veckl, veckc, ngp, igpig, vgpql, vgpqc, gkqc, tpgkqc)
          call timesec( t1)
          tgk = tgk + t1 - t0
          allocate( evecmt1( nlmlammax, pwmat_fst1:pwmat_lst1, natmtot))
          allocate( evecmt2( nlmlammax, pwmat_fst2:pwmat_lst2, natmtot))
          allocate( auxmat( pwmat_fst1:pwmat_lst1, nlmlammax))
          call timesec( t0)
          call pwmat_getevecmt( veckl, pwmat_fst1, pwmat_lst1, evecmt1)
          call pwmat_getevecmt( veckql, pwmat_fst2, pwmat_lst2, evecmt2)
          call timesec( t1)
          tio = tio + t1 - t0

          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)

              do ig = 1, ng
                call zgemm( 'c', 'n', pwmat_nst1, nlmlam( is), nlmlam( is), zone, &
                       evecmt1( :, :, ias), nlmlammax, &
                       rignt( :, :, ias, ig), nlmlammax, zzero, &
                       auxmat, pwmat_nst1)
                call zgemm( 'n', 'n', pwmat_nst1, pwmat_nst2, nlmlam( is), zone, &
                       auxmat, pwmat_nst1, &
                       evecmt2( :, :, ias), nlmlammax, zone, &
                       pwmat( :, :, ig), pwmat_nst1)
              end do
            end do
          end do
          call timesec( t0)
          tmt = tmt + t0 - t1
          deallocate( evecmt1, evecmt2, auxmat)
          
          !--------------------------------------!
          !     interstitial matrix elements     !
          !--------------------------------------!
 
          allocate( cfunmat( ngp, ngpq))
          allocate( auxmat( pwmat_fst1:pwmat_lst1, ngpq))
          do ig = 1, ng
            gs = shift + vecgl( :, ig)
            g = 0
            cfunmat(:,:) = zzero
            call timesec( t0)
#ifdef USEOMP
!$omp parallel default( shared) private( igk, igq, g)
!$omp do collapse( 2)
#endif
            do igk = 1, ngp
              do igq = 1, ngpq
                g = ivg( :, igpig( igk)) - ivg( :, igpqig( igq)) + gs
                if( (g(1) .ge. intgv(1,1)) .and. (g(1) .le. intgv(1,2)) .and. &
                    (g(2) .ge. intgv(2,1)) .and. (g(2) .le. intgv(2,2)) .and. &
                    (g(3) .ge. intgv(3,1)) .and. (g(3) .le. intgv(3,2))) then
                  cfunmat( igk, igq) = cfunig( ivgig( g(1), g(2), g(3)))
                end if
              end do
            end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
            call timesec( t1)
            tgg = tgg + t1 - t0

            call zgemm( 'c', 'n', pwmat_nst1, ngpq, ngp, zone, &
                   evec1, nmatmax_ptr, &
                   cfunmat, ngp, zzero, &
                   auxmat, pwmat_nst1)
            call zgemm( 'n', 'n', pwmat_nst1, pwmat_nst2, ngpq, zone, &
                   auxmat, pwmat_nst1, &
                   evec2, nmatmax_ptr, zone, &
                   pwmat( :, :, ig), pwmat_nst1)
            call timesec( t0)
            tir = tir + t0 - t1
          end do

          deallocate( auxmat, cfunmat, igpig, igpqig)

          return
      end subroutine pwmat_genpwmat

      subroutine pwmat_putevecmt( vpl, fst, lst, evecmt)
          use m_getunit
          real(8), intent( in)    :: vpl(3)
          integer, intent( in)    :: fst, lst
          complex(8), intent( in) :: evecmt( nlmlammax, fst:lst, natmtot)

          integer :: ik, un, recl

          call findkptinset( vpl, pwmat_kset, un, ik)

          inquire( iolength=recl) pwmat_kset%vkl( :, ik), nlmlammax, fst, lst, evecmt
          call getunit( un)
          open( un, file=trim( pwmat_fname)//'MT'//trim( filext), action='write', form='unformatted', access='direct', recl=recl)
          write( un, rec=ik) pwmat_kset%vkl( :, ik), nlmlammax, fst, lst, evecmt
          close( un)
          !if( ik .eq. 12) then
          !  write(*,'(3f13.6,3i)') pwmat_kset%vkl( :, ik), nlmlammax, fst, lst
          !  call plotmat( evecmt( :, :, 3))
          !  write(*,*)
          !end if

          return
      end subroutine pwmat_putevecmt

      subroutine pwmat_putevecir( vpl, fst, lst, evecir)
          use m_getunit
          real(8), intent( in)    :: vpl(3)
          integer, intent( in)    :: fst, lst
          complex(8), intent( in) :: evecir( pwmat_fftmap%ngrtot+1, fst:lst)

          integer :: ik, un, recl

          call findkptinset( vpl, pwmat_kset, un, ik)

          inquire( iolength=recl) pwmat_kset%vkl( :, ik), pwmat_fftmap%ngrtot+1, fst, lst, evecir
          call getunit( un)
          open( un, file=trim( pwmat_fname)//'IR'//trim( filext), action='write', form='unformatted', access='direct', recl=recl)
          write( un, rec=ik) pwmat_kset%vkl( :, ik), pwmat_fftmap%ngrtot+1, fst, lst, evecir
          close( un)

          return
      end subroutine pwmat_putevecir

      subroutine pwmat_getevecmt( vpl, fst, lst, evecmt)
          use m_getunit
          real(8), intent( in)     :: vpl(3)
          integer, intent( in)     :: fst, lst
          complex(8), intent( out) :: evecmt( nlmlammax, fst:lst, natmtot)

          integer :: ik, un, recl, nlmlammax_, fst_, lst_
          real(8) :: vpl_(3)
          logical :: exist

          complex(8), allocatable :: tmp(:,:,:)

          call findkptinset( vpl, pwmat_kset, un, ik)

          inquire( file=trim( pwmat_fname)//'MT'//trim( filext), exist=exist)
          if( exist) then
            inquire( iolength=recl) vpl_, nlmlammax_, fst_, lst_
            call getunit( un)
            open( un, file=trim( pwmat_fname)//'MT'//trim( filext), action='read', form='unformatted', access='direct', recl=recl)
            read( un, rec=1) vpl_, nlmlammax_, fst_, lst_
            close( un)
            if( nlmlammax_ .ne. nlmlammax) then
              write(*,*)
              write(*,'("Error (pwmat_geteveshort): Different number of basis functions:")')
              write(*,'(" current:",i8)') nlmlammax
              write(*,'(" file   :",i8)') nlmlammax_
              stop
            end if
            if( (fst .lt. fst_) .or. (lst .gt. lst_)) then
              write(*,*)
              write(*,'("Error (pwmat_geteveshort): Band indices out of range:")')
              write(*,'(" current:",2i8)') fst, lst
              write(*,'(" file   :",2i8)') fst_, lst_
              stop
            end if
            allocate( tmp( nlmlammax, fst_:lst_, natmtot))
            inquire( iolength=recl) vpl_, nlmlammax_, fst_, lst_, tmp
            call getunit( un)
            open( un, file=trim( pwmat_fname)//'MT'//trim( filext), action='read', form='unformatted', access='direct', recl=recl)
            read( un, rec=ik) vpl_, nlmlammax_, fst_, lst_, tmp
            close( un)
            evecmt = tmp( :, fst:lst, :)
            deallocate( tmp)
            !if( ik .eq. 12) then
            !  write(*,'(3f13.6,3i)') vpl_, nlmlammax_, fst_, lst_
            !  call plotmat( evecmt( :, :, 3))
            !  write(*,*)
            !end if
          else
            write(*,*)
            write(*,'("Error (pwmat_getevecmt): File does not exist:",a)') trim( pwmat_fname)//'MT'//trim( filext)
            stop
          end if

          return
      end subroutine pwmat_getevecmt

      subroutine pwmat_getevecir( vpl, fst, lst, evecir)
          use m_getunit
          real(8), intent( in)     :: vpl(3)
          integer, intent( in)     :: fst, lst
          complex(8), intent( out) :: evecir( pwmat_fftmap%ngrtot+1, fst:lst)

          integer :: ik, un, recl, ngrtot_, fst_, lst_
          real(8) :: vpl_(3)
          logical :: exist

          complex(8), allocatable :: tmp(:,:)

          call findkptinset( vpl, pwmat_kset, un, ik)

          inquire( file=trim( pwmat_fname)//'IR'//trim( filext), exist=exist)
          if( exist) then
            inquire( iolength=recl) vpl_, ngrtot_, fst_, lst_
            call getunit( un)
            open( un, file=trim( pwmat_fname)//'IR'//trim( filext), action='read', form='unformatted', access='direct', recl=recl)
            read( un, rec=1) vpl_, ngrtot_, fst_, lst_
            close( un)
            if( ngrtot_ .ne. pwmat_fftmap%ngrtot+1) then
              write(*,*)
              write(*,'("Error (pwmat_getevecir): Different number of plane waves:")')
              write(*,'(" current:",i8)') pwmat_fftmap%ngrtot+1
              write(*,'(" file   :",i8)') ngrtot
              stop
            end if
            if( (fst .lt. fst_) .or. (lst .gt. lst_)) then
              write(*,*)
              write(*,'("Error (pwmat_geteveshort): Band indices out of range:")')
              write(*,'(" current:",2i8)') fst, lst
              write(*,'(" file   :",2i8)') fst_, lst_
              stop
            end if
            allocate( tmp( pwmat_fftmap%ngrtot+1, fst_:lst_))
            inquire( iolength=recl) vpl_, ngrtot_, fst_, lst_, tmp
            call getunit( un)
            open( un, file=trim( pwmat_fname)//'IR'//trim( filext), action='read', form='unformatted', access='direct', recl=recl)
            read( un, rec=ik) vpl_, ngrtot_, fst_, lst_, evecir
            close( un)
          else
            write(*,*)
            write(*,'("Error (pwmat_getevecir): File does not exist:",a)') trim( pwmat_fname)//'IR'//trim( filext)
            stop
          end if

          return
      end subroutine pwmat_getevecir

      subroutine pwmat_destroy
          use m_getunit
          integer :: un
          logical :: exist

          if( mpiglobal%rank .eq. 0) then
            inquire( file=trim( pwmat_fname)//'MT'//trim( filext), exist=exist)
            if( exist) then
              call getunit( un)
              open( un, file=trim( pwmat_fname)//'MT'//trim( filext))
              close( un, status='delete')
            end if
            inquire( file=trim( pwmat_fname)//'IR'//trim( filext), exist=exist)
            if( exist) then
              call getunit( un)
              open( un, file=trim( pwmat_fname)//'IR'//trim( filext))
              close( un, status='delete')
            end if
          end if

          call timesec( te)
          ttot = te - ts
          if( mpiglobal%rank .eq. 0 .and. .false.) then
            write(*,'("          PWMAT TIMING          ")')
            write(*,'("================================")')
            write(*,'("initialization   : ",f13.6)') tinit
            write(*,'("preparation      : ",f13.6)') tprep
            write(*,'("radial integrals : ",f13.6)') tri
            write(*,'("gaunt integrals  : ",f13.6)') trignt
            write(*,'("G+k vectors      : ",f13.6)') tgk
            write(*,'("muffin-tin part  : ",f13.6)') tmt
            write(*,'("char. fun. GG    : ",f13.6)') tgg
            write(*,'("interstitial part: ",f13.6)') tir
            write(*,'("I/O              : ",f13.6)') tio
            write(*,'("--------------------------------")')
            write(*,'("sum              : ",f13.6)') tinit + tri + trignt + tmt + tgk + tprep + tgg + tir + tio
            write(*,'("rest             : ",f13.6)') ttot - (tinit + tri + trignt + tmt + tgk + tprep + tgg + tir + tio)
            write(*,'("--------------------------------")')
            write(*,'("total            : ",f13.6)') ttot
            write(*,'("================================")')
          end if

          if( allocated( nlam)) deallocate( nlam)
          if( allocated( nlmlam)) deallocate( nlmlam)
          if( allocated( lam2apwlo)) deallocate( lam2apwlo)
          if( allocated( idxlmlam)) deallocate( idxlmlam)
          if( allocated( lm2l)) deallocate( lm2l)
          if( allocated( rignt)) deallocate( rignt)
          
          return
      end subroutine pwmat_destroy
end module mod_pwmat
