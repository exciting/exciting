module mod_symmat
  use modinput
  use mod_atoms
  use mod_APW_LO
  use mod_muffin_tin
  use mod_constants
  use mod_eigensystem
  use mod_symmetry
  use mod_lattice, only : bvec
  use mod_Gvector, only : ngvec, gc, cfunig, ivg, ivgig, ngrtot, ngrid, intgv
  use mod_Gkvector, only : gkmax, ngkmax
  use mod_kpointset
  use m_plotmat

  implicit none
  private

  integer :: nlammax, nlmlammax, sym_lmaxapw
  real(8) :: sym_gmax

  integer, allocatable :: nlam(:,:), nlmlam(:), lam2apwlo(:,:,:), apw2lam(:,:,:), lo2lam(:,:,:), idxlmlam(:,:,:), lm2l(:)
  real(8), allocatable :: sym_radolp(:,:,:,:)

  public :: symmat_init, symmat_gensymmat, symmat_destroy

! variable

! methods
  contains
      subroutine symmat_init( sym_lmaxapw_)
          integer, intent( in) :: sym_lmaxapw_

          integer :: l1, m1, o1, lm1, i, is, ia, ias, ilo1, ilo2, lam

          sym_lmaxapw = min( sym_lmaxapw_, input%groundstate%lmaxapw)
          ! count combined (l,m,o) indices and build index maps
          allocate( nlam( 0:sym_lmaxapw, nspecies))
          allocate( lam2apwlo( apwordmax+nlomax, 0:sym_lmaxapw, nspecies))
          allocate( apw2lam( apwordmax, 0:sym_lmaxapw, nspecies))
          allocate( lo2lam( nlomax, 0:sym_lmaxapw, nspecies))
          nlammax = 0
          do is = 1, nspecies
            lam2apwlo( :, :, is) = 0
            apw2lam( :, :, is) = 0
            lo2lam( :, :, is) = 0
            nlam( :, is) = 0
            do l1 = 0, sym_lmaxapw
              do o1 = 1, apword( l1, is)
                nlam( l1, is) = nlam( l1, is)+1
                apw2lam( o1, l1, is) = nlam( l1, is)
                lam2apwlo( nlam( l1, is), l1, is) = o1
              end do
            end do
            do ilo1 = 1, nlorb( is)
              l1 = lorbl( ilo1, is)
              if( l1 .le. sym_lmaxapw) then
                nlam( l1, is) = nlam( l1, is) + 1
                lo2lam( ilo1, l1, is) = nlam( l1, is)
                lam2apwlo( nlam( l1, is), l1, is) = -ilo1
              end if
              nlammax = max( nlammax, nlam( l1, is))
            end do
          end do

          allocate( lm2l( (sym_lmaxapw + 1)**2))
          allocate( nlmlam( nspecies))
          allocate( idxlmlam( (sym_lmaxapw + 1)**2, nlammax, nspecies)) 
          idxlmlam(:,:,:) = 0
          nlmlammax = 0
          nlmlam(:) = 0
          do l1 = 0, sym_lmaxapw
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
          write(*,'("symmat_init: nlmlammax = :",I4)') nlmlammax

          ! radial overlap-intergrals
          ! generate radial functions
          !call readstate
          call readfermi
          !call linengy
          !call genapwfr
          !call genlofr
          !call olprad

          if( allocated( sym_radolp)) deallocate( sym_radolp)
          allocate( sym_radolp( nlammax, nlammax, 0:sym_lmaxapw, natmtot))
          sym_radolp(:,:,:,:) = 0.d0
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              ! APW-APW integral
              do l1 = 0, sym_lmaxapw
                do o1 = 1, apword( l1, is)
                  sym_radolp( apw2lam( o1, l1, is), apw2lam( o1, l1, is), l1, ias) = 1.d0
                end do
              end do
              do ilo1 = 1, nlorb( is)
                l1 = lorbl( ilo1, is)
                if( l1 .le. sym_lmaxapw) then
                  do o1 = 1, apword( l1, is)
                    if( (apw2lam( o1, l1, is) .gt. 0) .and. (lo2lam( ilo1, l1, is) .gt. 0)) then
                      sym_radolp( apw2lam( o1, l1, is), lo2lam( ilo1, l1, is), l1, ias) = oalo( o1, ilo1, ias)
                      sym_radolp( lo2lam( ilo1, l1, is), apw2lam( o1, l1, is), l1, ias) = oalo( o1, ilo1, ias)
                    end if
                  end do
                  do ilo2 = 1, nlorb( is)
                    if( lorbl( ilo2, is) .eq. l1) then
                      if( (lo2lam( ilo1, l1, is) .gt. 0) .and. (lo2lam( ilo2, l1, is) .gt. 0)) then
                        sym_radolp( lo2lam( ilo1, l1, is), lo2lam( ilo2, l1, is), l1, ias) = ololo( ilo1, ilo2, ias)
                        sym_radolp( lo2lam( ilo2, l1, is), lo2lam( ilo1, l1, is), l1, ias) = ololo( ilo2, ilo1, ias)
                      end if
                    end if
                  end do
                end if
              end do
            end do
          end do

          return
      end subroutine symmat_init
      
      subroutine symmat_gensymmat( vkl1, vkl2, isym, sym_fst1, sym_lst1, sym_fst2, sym_lst2, evec1, evec2, symmat)
          real(8), intent( in) :: vkl1(3), vkl2(3)
          integer, intent( in) :: isym, sym_fst1, sym_lst1, sym_fst2, sym_lst2
          complex(8), intent( in) :: evec1( nmatmax, sym_fst1:sym_lst1), evec2( nmatmax, sym_fst2:sym_lst2)
          complex(8), intent( out) :: symmat( sym_fst1:sym_lst1, sym_fst2:sym_lst2)

          integer :: ngk1, ngk2, i, is, ia, ias, jas, l, m, lm, m2, o, ilo, lam, lspl, ilspl
          integer :: sym_nst1, sym_nst2
          real(8) :: t1, veckql(3), veckqc(3), vkc1(3), vkc2(3), rotl(3,3), rotc(3,3), tau(3), vreal(3)
          
          integer :: shift(3), igk, igq, ig1, ig2, g(3)

          integer, allocatable :: igkig1(:), igkig2(:) 
          real(8), allocatable :: vecgkl2(:,:), vecgkc2(:,:), gkc2(:), tpgkc2(:,:)
          complex(8), allocatable :: sfacgk2(:,:), apwalm(:,:,:,:), rint(:,:)
          complex(8), allocatable :: auxmat(:,:), auxvec(:), evecshort1(:,:), evecshort2(:,:,:), cfunmat(:,:)

          complex(8) :: zdotc
          complex(8), external :: getdlmm
          
          sym_nst1 = sym_lst1-sym_fst1+1
          sym_nst2 = sym_lst2-sym_fst2+1
          symmat(:,:) = zzero
 
          lspl = lsplsymc( isym)
          ilspl = isymlat( lspl)
          rotl = dble( symlat( :, :, ilspl))
          rotc = dble( symlatc( :, :, ilspl))
          tau = vtlsymc( :, isym)
          !write(*,'(i3,3f13.6)') isym, tau

          call r3mv( bvec, vkl1, vkc1)
          call r3mv( bvec, vkl2, vkc2)

          ! check if symmetry is identity
          !if( isym .eq. 1) then
          !  do i = 0, min( sym_nst1, sym_nst2)-1
          !    symmat( sym_fst1+i, sym_fst2+i) = zone
          !  end do
          !  return
          !end if

          !--------------------------------------!
          !      muffin-tin matrix elements      !
          !--------------------------------------!

          allocate( igkig1( ngkmax), igkig2( ngkmax), vecgkl2( 3, ngkmax), vecgkc2( 3, ngkmax), gkc2( ngkmax), tpgkc2( 2, ngkmax))
          allocate( sfacgk2( ngkmax, natmtot))
          allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot))
          allocate( evecshort2( sym_fst2:sym_lst2, nlmlammax, natmtot))
          allocate( auxvec( max( sym_nst1, sym_nst2)))

          ! generate the G+k2-vectors
          call gengpvec( vkl2, vkc2, ngk2, igkig2, vecgkl2, vecgkc2, gkc2, tpgkc2)
          ! generate the structure factors
          call gensfacgp( ngk2, vecgkc2, ngkmax, sfacgk2)
          ! find matching coefficients for k-point k2
          call match( ngk2, gkc2, tpgkc2, sfacgk2, apwalm)

          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)
              jas = idxas( ieqatom( ia, is, isym), is)
              !write(*,*) "muffin tin", ias, jas
              evecshort2( :, :, ias) = zzero
              do l = 0, sym_lmaxapw
                !write(*,*) ias, l
                !call plotmat( zone*sym_radolp( :, :, l, ias))
                !write(*,*)
                do m = -l, l
                  lm = idxlm( l, m)
                  do lam = 1, nlam( l, is)
                    ! lam belongs to apw
                    if( lam2apwlo( lam, l, is) .gt. 0) then
                      o = lam2apwlo( lam, l, is)
                      call zgemv( 't', ngk2, sym_nst2, zone, &
                             evec2, nmatmax, &
                             apwalm( :, o, lm, jas), 1, zzero, &
                             auxvec, 1)
                      do m2 = -l, l
                        evecshort2( :, idxlmlam( idxlm( l, m2), lam, is), ias) = evecshort2( :, idxlmlam( idxlm( l, m2), lam, is), ias) + &
                            getdlmm( rotc, l, m, m2)*auxvec( 1:sym_nst2)
                      end do
                    end if
                    ! lam belongs to lo
                    if( lam2apwlo( lam, l, is) .lt. 0) then
                      ilo = -lam2apwlo( lam, l, is)
                      do m2 = -l, l
                        evecshort2( :, idxlmlam( idxlm( l, m2), lam, is), ias) = evecshort2( :, idxlmlam( idxlm( l, m2), lam, is), ias) + &
                            getdlmm( rotc, l, m, m2)*evec2( ngk2+idxlo( lm, ilo, jas), :)
                      end do
                    end if
                  end do
                end do
              end do
            end do
          end do

          ! generate the G+k1-vectors
          call gengpvec( vkl1, vkc1, ngk1, igkig1, vecgkl2, vecgkc2, gkc2, tpgkc2)
          ! generate the structure factors
          call gensfacgp( ngk1, vecgkc2, ngkmax, sfacgk2)
          ! find matching coefficients for k-point k1
          call match( ngk1, gkc2, tpgkc2, sfacgk2, apwalm)
          
          write(*,*) "ngk", ngk1, ngk2

          allocate( evecshort1( sym_fst1:sym_lst1, nlmlammax))
          allocate( auxmat( sym_fst1:sym_lst1, nlmlammax))
          allocate( rint( nlmlammax, nlmlammax))
          
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)

              evecshort1(:,:) = zzero
              rint = zzero
              do l = 0, sym_lmaxapw
                do m = -l, l
                  lm = idxlm( l, m)
                  do lam = 1, nlam( l, is)
                    ! lam belongs to apw
                    if( lam2apwlo( lam, l, is) .gt. 0) then
                      o = lam2apwlo( lam, l, is)
                      call zgemv( 't', ngk1, sym_nst1, zone, &
                             evec1, nmatmax, &
                             apwalm( :, o, lm, ias), 1, zzero, &
                             auxvec, 1)
                      evecshort1( :, idxlmlam( lm, lam, is)) = conjg( auxvec( 1:sym_nst1))
                    end if
                    ! lam belongs to lo
                    if( lam2apwlo( lam, l, is) .lt. 0) then
                      ilo = -lam2apwlo( lam, l, is)
                      evecshort1( :, idxlmlam( lm, lam, is)) = conjg( evec1( ngk1+idxlo( lm, ilo, ias), :))
                    end if
                    ! set up integral matrix
                    do m2 = 1, nlam( l, is)
                      rint( idxlmlam( lm, lam, is), idxlmlam( lm, m2, is)) = cmplx( sym_radolp( lam, m2, l, ias), 0, 8)
                    end do
                  end do
                end do
              end do

              call zgemm( 'n', 'n', sym_nst1, nlmlam( is), nlmlam( is), zone, &
                     evecshort1, sym_nst1, &
                     rint, nlmlammax, zzero, &
                     auxmat, sym_nst1)
              call zgemm( 'n', 't', sym_nst1, sym_nst2, nlmlam( is), zone, &
                     auxmat, sym_nst1, &
                     evecshort2( :, :, ias), sym_nst2, zone, &
                     symmat, sym_nst1)
            end do
          end do
          
          deallocate( evecshort1, evecshort2, auxmat, auxvec, rint)
          
          !--------------------------------------!
          !     interstitial matrix elements     !
          !--------------------------------------!
 
          allocate( cfunmat( ngk1, ngk2))
          allocate( auxmat( sym_fst1:sym_lst1, ngk2))
          cfunmat(:,:) = zzero
          do igq = 1, ngk2
            t1 = -twopi*dot_product( vecgkl2( :, igq), tau)
            ig2 = igkig2( igq)
            call r3mtv( rotl, dble( ivg( :, ig2)), vreal)
            do igk = 1, ngk1
              ig1 = igkig1( igk)
              g = ivg( :, ig1) - nint( vreal)
              if( (g(1) .ge. intgv(1,1)) .and. (g(1) .le. intgv(1,2)) .and. &
                  (g(2) .ge. intgv(2,1)) .and. (g(2) .le. intgv(2,2)) .and. &
                  (g(3) .ge. intgv(3,1)) .and. (g(3) .le. intgv(3,2))) then
                cfunmat( igk, igq) = cfunig( ivgig( g(1), g(2), g(3)))*cmplx( cos( t1), sin( t1), 8)
              end if
            end do
          end do

          call zgemm( 'c', 'n', sym_nst1, ngk2, ngk1, zone, &
                 evec1, nmatmax, &
                 cfunmat, ngk1, zzero, &
                 auxmat, sym_nst1)
          call zgemm( 'n', 'n', sym_nst1, sym_nst2, ngk2, zone, &
                 auxmat, sym_nst1, &
                 evec2, nmatmax, zone, &
                 symmat, sym_nst1)

          deallocate( vecgkl2, vecgkc2, gkc2, tpgkc2, sfacgk2, apwalm) 
          deallocate( auxmat, cfunmat, igkig1, igkig2)

          return
      end subroutine symmat_gensymmat

      subroutine symmat_destroy
          if( allocated( nlam)) deallocate( nlam)
          if( allocated( nlmlam)) deallocate( nlmlam)
          if( allocated( lam2apwlo)) deallocate( lam2apwlo)
          if( allocated( idxlmlam)) deallocate( idxlmlam)
          if( allocated( lm2l)) deallocate( lm2l)
          
          return
      end subroutine symmat_destroy
end module mod_symmat
