module mod_pwmat
  use modinput
  use mod_atoms
  use mod_APW_LO
  use mod_muffin_tin
  use mod_constants
  use mod_eigensystem, only : idxlo, nmatmax
  use mod_lattice, only : bvec
  use mod_Gvector, only : ngvec, gc, cfunig, ivg, ivgig, ngrtot, ngrid, intgv
  use mod_Gkvector, only : gkmax, ngkmax
  use m_plotmat

  implicit none
  private

  integer :: nlammax, nlmlammax, pw_lmaxapw, pw_lmaxexp
  real(8) :: pwmat_gmax
  real(8) :: vecql(3), vecqc(3), vecgc(3)
  integer :: vecgl(3)

  integer, allocatable :: nlam(:,:), nlmlam(:), lam2apwlo(:,:,:), idxlmlam(:,:,:), lm2l(:), idxgnt(:,:,:)
  real(8), allocatable :: listgnt(:,:,:)
  complex(8), allocatable :: rignt(:,:,:)

  public :: pwmat_init, pwmat_init_qg, pwmat_genpwmat, pwmat_destroy

! variable

! methods
  contains
      subroutine pwmat_init( pw_lmaxapw_, pw_lmaxexp_)
          integer, intent( in) :: pw_lmaxapw_, pw_lmaxexp_

          integer :: l1, l2, l3, m1, m2, m3, o1, lm1, lm2, lm3, i, is, ilo, lam
          real(8) :: gnt, gaunt

          external :: gaunt

          pw_lmaxapw = min( pw_lmaxapw_, input%groundstate%lmaxapw)
          pw_lmaxexp = min( pw_lmaxexp_, pw_lmaxapw)
          ! count combined (l,m,o) indices and build index maps
          allocate( nlam( 0:pw_lmaxapw, nspecies))
          allocate( lam2apwlo( apwordmax+nlomax, 0:pw_lmaxapw, nspecies))
          nlammax = 0
          do is = 1, nspecies
            do l1 = 0, pw_lmaxapw
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

          allocate( lm2l( (pw_lmaxapw + 1)**2))
          allocate( nlmlam( nspecies))
          allocate( idxlmlam( (pw_lmaxapw + 1)**2, nlammax, nspecies)) 
          idxlmlam(:,:,:) = 0
          nlmlammax = 0
          nlmlam(:) = 0
          do l1 = 0, pw_lmaxapw
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

          ! build non-zero Gaunt list
          allocate( idxgnt( pw_lmaxapw + 1, (pw_lmaxapw + 1)**2, (pw_lmaxapw + 1)**2), &
                    listgnt( pw_lmaxapw + 1, (pw_lmaxapw + 1)**2, (pw_lmaxapw + 1)**2))
          idxgnt(:,:,:) = 0
          listgnt(:,:,:) = 0.d0
          do l3 = 0, pw_lmaxapw
            do m3 = -l3, l3
              lm3 = idxlm( l3, m3)

              do l1 = 0, pw_lmaxapw
                do m1 = -l1, l1
                  lm1 = idxlm( l1, m1)

                  i = 0
                  do l2 = 0, pw_lmaxexp
                    do m2 = -l2, l2
                      lm2 = idxlm( l2, m2)
                      gnt = gaunt( l1, l2, l3, m1, m2, m3)
                      if( abs( gnt) .gt. 1.d-20) then
                        i = i + 1
                        listgnt( i, lm1, lm3) = gnt
                        idxgnt( i, lm1, lm3) = lm2
                      end if
                    end do
                  end do

                end do
              end do
            
            end do
          end do

          return
      end subroutine pwmat_init
      
      subroutine pwmat_init_qg( vecql_, vecgl_)
          real(8), intent( in) :: vecql_(3)
          integer, intent( in) :: vecgl_(3)

          integer :: l1, l2, l3, m1, m2, m3, o1, o2, lm1, lm2, lm3, is, ia, ias, i
          integer :: lam1, lam2, ilo1, ilo2, idxlo1, idxlo2, idxlostart
          real(8) :: qgc, tp(2)
          complex(8) :: strf, ylm( (pw_lmaxexp + 1)**2)

          real(8), allocatable :: radint(:,:,:)

          vecql = vecql_ + dble( vecgl_)
          call r3frac( input%structure%epslat, vecql, vecgl)
          !write(*,'("vq and vg: ",3f13.6,3i)') vecql, vecgl
          call r3mv( bvec, vecql, vecqc)
          call r3mv( bvec, dble( vecgl), vecgc)
          
          ! compute radial integrals times Gaunt and expansion prefactor
          if( allocated( rignt)) deallocate( rignt)
          allocate( rignt( nlmlammax, nlmlammax, natmtot))
          rignt(:,:,:) = zzero
          allocate( radint( nlammax, nlammax, 0:pw_lmaxexp))
          ! generate radial functions
          !call readstate
          call readfermi
          !call linengy
          !call genapwfr
          !call genlofr
          ! multiply radial integral with Gaunt coefficients and expansion prefactor
          idxlostart = 0
          ! generate spherical harmonics
          call sphcrd( vecqc(:) + vecgc(:), qgc, tp)
          call genylm( pw_lmaxexp, tp, ylm)

          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              strf = fourpi*exp( -zi*dot_product( vecqc(:) + vecgc(:), atposc( :, ia, is)))
              ! generate radial integral
              do l1 = 0, pw_lmaxapw
                do l2 = 0, pw_lmaxapw
                  call pwmat_genri( is, ia, l1, l2, radint)
                  !write(*,'("init_qg: ",3I)') ias, l1, l2
                  do m1 = -l1, l1
                    lm1 = idxlm( l1, m1)
                    do m2 = -l2, l2
                      lm2 = idxlm( l2, m2)
                      
                      do lam1 = 1, nlam( l1, is)
                        do lam2 = 1, nlam( l2, is)
                        
                          i = 1
                          do while( idxgnt( i, lm1, lm2) .ne. 0)
                            lm3 = idxgnt( i, lm1, lm2)
                            l3 = lm2l( lm3)
                            rignt( idxlmlam( lm1, lam1, is), idxlmlam( lm2, lam2, is), ias) = & 
                                 rignt( idxlmlam( lm1, lam1, is), idxlmlam( lm2, lam2, is), ias) + &
                                 conjg( zil( l3))*conjg( ylm( lm3))*cmplx( listgnt( i, lm1, lm2)*radint( lam1, lam2, l3), 0, 8)
                            i = i + 1
                          end do
                          
                        end do
                      end do

                    end do
                  end do
                end do
              end do

              rignt( :, :, ias) = strf*rignt( :, :, ias)
            end do ! atoms
          end do !species

          deallocate( radint)
          return
      end subroutine pwmat_init_qg
      
      subroutine pwmat_genri( is, ia, l1, l2, radint)
          integer, intent( in) :: is, ia, l1, l2
          real(8), intent( out) :: radint( nlammax, nlammax, 0:pw_lmaxexp)

          integer :: ias, ir, nr, l3, o1, o2, ilo1, ilo2, lam1, lam2
          real(8) :: x, qgc, tp(2)

          real(8), allocatable :: jlqgr(:,:), fr(:), gf(:), cf(:,:)
          
          radint(:,:,:) = 0.d0
          
          ias = idxas( ia, is)
          call sphcrd( vecqc(:) + vecgc(:), qgc, tp)
          nr = nrmt( is)
          ! generate spherical Bessel functions
          allocate( jlqgr( 0:pw_lmaxexp, nr))
          do ir = 1, nr
            x = qgc*spr( ir, is)
            call sbessel( pw_lmaxexp, x, jlqgr( :, ir))
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
                  do l3 = 0, pw_lmaxexp
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
                  do l3 = 0, pw_lmaxexp
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
                  do l3 = 0, pw_lmaxexp
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
                  do l3 = 0, pw_lmaxexp
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
      
      subroutine pwmat_genpwmat( veckl, veckc, pw_fst1, pw_lst1, pw_fst2, pw_lst2, evec1, evec2, pwmat)
          real(8), intent( in) :: veckl(3), veckc(3)
          integer, intent( in) :: pw_fst1, pw_lst1, pw_fst2, pw_lst2
          complex(8), intent( in) :: evec1( nmatmax, pw_fst1:pw_lst1), evec2( nmatmax, pw_fst2:pw_lst2)
          complex(8), intent( out) :: pwmat( pw_fst1:pw_lst1, pw_fst2:pw_lst2)

          integer :: ngknr, ngkq, i, is, ia, ias, l, m, lm, o, ilo, lam
          integer :: pw_nst1, pw_nst2
          real(8) :: t1, veckql(3), veckqc(3)
          
          integer :: shift(3), igk, igq, ig1, ig2, g(3)

          integer, allocatable :: igkignr(:), igkqig(:) 
          real(8), allocatable :: vecgkql(:,:), vecgkqc(:,:), gkqc(:), tpgkqc(:,:)
          complex(8), allocatable :: sfacgkq(:,:), apwalm(:,:,:,:)
          complex(8), allocatable :: auxmat(:,:), auxvec(:), evecshort1(:,:), evecshort2(:,:,:), cfunmat(:,:)

          complex(8) :: zdotc
          
          pw_nst1 = pw_lst1-pw_fst1+1
          pw_nst2 = pw_lst2-pw_fst2+1
          pwmat(:,:) = zzero

          ! check if q-vector is zero
          t1 = norm2( vecql + dble( vecgl))
          if( t1 .lt. input%structure%epslat) then
            do i = 0, min( pw_nst1, pw_nst2)-1
              pwmat( pw_fst1+i, pw_fst2+i) = zone
            end do
            return
          end if

          ! k+q-vector in lattice coordinates
          veckql = veckl + vecql
          ! map vector components to [0,1) interval
          call r3frac( input%structure%epslat, veckql, g)
          ! k+q-vector in Cartesian coordinates
          call r3mv( bvec, veckql, veckqc)
          
          !--------------------------------------!
          !      muffin-tin matrix elements      !
          !--------------------------------------!

          allocate( igkignr( ngkmax), igkqig( ngkmax), vecgkql( 3, ngkmax), vecgkqc( 3, ngkmax), gkqc( ngkmax), tpgkqc( 2, ngkmax))
          allocate( sfacgkq( ngkmax, natmtot))
          allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
          allocate( evecshort2( pw_fst2:pw_lst2, nlmlammax, natmtot))
          allocate( auxvec( max( pw_nst1, pw_nst2)))

          ! generate the G+k+q-vectors
          call gengpvec( veckql, veckqc, ngkq, igkqig, vecgkql, vecgkqc, gkqc, tpgkqc)
          ! generate the structure factors
          call gensfacgp( ngkq, vecgkqc, ngkmax, sfacgkq)
          ! find matching coefficients for k-point k+q
          call match( ngkq, gkqc, tpgkqc, sfacgkq, apwalm)

          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              evecshort2( :, :, ias) = zzero
              do lm = 1, (pw_lmaxapw + 1)**2
                l = lm2l( lm)
                do lam = 1, nlam( l, is)
                  ! lam belongs to apw
                  if( lam2apwlo( lam, l, is) .gt. 0) then
                    o = lam2apwlo( lam, l, is)
                    call zgemv( 't', ngkq, pw_nst2, zone, &
                           evec2, nmatmax, &
                           apwalm( :, o, lm, ias), 1, zzero, &
                           auxvec, 1)
                    evecshort2( :, idxlmlam( lm, lam, is), ias) = auxvec( 1:pw_nst2)
                  end if
                  ! lam belongs to lo
                  if( lam2apwlo( lam, l, is) .lt. 0) then
                    ilo = -lam2apwlo( lam, l, is)
                    evecshort2( :, idxlmlam( lm, lam, is), ias) = evec2( ngkq+idxlo( lm, ilo, ias), :)
                  end if
                end do
              end do
            end do
          end do

          ! generate the G+k-vectors
          call gengpvec( veckl, veckc, ngknr, igkignr, vecgkql, vecgkqc, gkqc, tpgkqc)
          ! generate the structure factors
          call gensfacgp( ngknr, vecgkqc, ngkmax, sfacgkq)
          ! find matching coefficients for k-point k
          call match( ngknr, gkqc, tpgkqc, sfacgkq, apwalm)
          
          allocate( evecshort1( pw_fst1:pw_lst1, nlmlammax))
          allocate( auxmat( pw_fst1:pw_lst1, nlmlammax))
          
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)

              evecshort1(:,:) = zzero
              do lm = 1, (pw_lmaxapw + 1)**2
                l = lm2l( lm)
                do lam = 1, nlam( l, is)
                  ! lam belongs to apw
                  if( lam2apwlo( lam, l, is) .gt. 0) then
                    o = lam2apwlo( lam, l, is)
                    call zgemv( 'T', ngknr, pw_nst1, zone, &
                         evec1, nmatmax, &
                         apwalm( :, o, lm, ias), 1, zzero, &
                         auxvec, 1)
                    evecshort1( :, idxlmlam( lm, lam, is)) = conjg( auxvec( 1:pw_nst1))
                  end if
                  ! lam belongs to lo
                  if( lam2apwlo( lam, l, is) .lt. 0) then
                    ilo = -lam2apwlo( lam, l, is)
                    evecshort1( :, idxlmlam( lm, lam, is)) = conjg( evec1( ngknr+idxlo( lm, ilo, ias), :))
                  end if
                end do
              end do

              call zgemm( 'N', 'N', pw_nst1, nlmlam( is), nlmlam( is), zone, &
                   evecshort1, pw_nst1, &
                   rignt( :, :, ias), nlmlammax, zzero, &
                   auxmat, pw_nst1)
              call zgemm( 'N', 'T', pw_nst1, pw_nst2, nlmlam( is), zone, &
                   auxmat, pw_nst1, &
                   evecshort2( :, :, ias), pw_nst2, zone, &
                   pwmat, pw_nst1)
            end do
          end do

          deallocate( vecgkql, vecgkqc, gkqc, tpgkqc, sfacgkq, apwalm) 
          deallocate( evecshort1, evecshort2, auxmat, auxvec)
          
          !--------------------------------------!
          !     interstitial matrix elements     !
          !--------------------------------------!
 
          shift = g
          
          allocate( cfunmat( ngknr, ngkq))
          allocate( auxmat( pw_fst1:pw_lst1, ngkq))
          cfunmat(:,:) = zzero
          do igk = 1, ngknr
            do igq = 1, ngkq
              ig1 = igkignr( igk)
              ig2 = igkqig( igq)
              g = ivg( :, ig1) - ivg( :, ig2) + shift - vecgl
              if( (g(1) .ge. intgv(1,1)) .and. (g(1) .le. intgv(1,2)) .and. &
                  (g(2) .ge. intgv(2,1)) .and. (g(2) .le. intgv(2,2)) .and. &
                  (g(3) .ge. intgv(3,1)) .and. (g(3) .le. intgv(3,2))) then
                cfunmat( igk, igq) = cfunig( ivgig( g(1), g(2), g(3)))
              end if
            end do
          end do

          call zgemm( 'c', 'n', pw_nst1, ngkq, ngknr, zone, &
                 evec1, nmatmax, &
                 cfunmat, ngknr, zzero, &
                 auxmat, pw_nst1)
          call zgemm( 'n', 'n', pw_nst1, pw_nst2, ngkq, zone, &
                 auxmat, pw_nst1, &
                 evec2, nmatmax, zone, &
                 pwmat, pw_nst1)

          deallocate( auxmat, cfunmat, igkignr, igkqig)

          return
      end subroutine pwmat_genpwmat

      subroutine pwmat_destroy
          if( allocated( nlam)) deallocate( nlam)
          if( allocated( nlmlam)) deallocate( nlmlam)
          if( allocated( lam2apwlo)) deallocate( lam2apwlo)
          if( allocated( idxlmlam)) deallocate( idxlmlam)
          if( allocated( lm2l)) deallocate( lm2l)
          if( allocated( rignt)) deallocate( rignt)
          if( allocated( idxgnt)) deallocate( idxgnt)
          if( allocated( listgnt)) deallocate( listgnt)
          
          return
      end subroutine pwmat_destroy
end module mod_pwmat
