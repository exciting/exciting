module m_ematqk
  use modmpi
  use modinput
  use mod_constants
  use mod_atoms
  use mod_APW_LO
  use mod_muffin_tin
  use mod_ematgrids
  use mod_eigenvalue_occupancy
  use modxs, only : fftmap_type

use m_writecmplxparts

  implicit none
  private

  logical :: initialized = .false.
  logical :: fftsaved = .false.

  integer(4) :: nlmomax, lmaxapw, lmaxexp
  integer(4) :: gs_lmaxapw_save

  ! Index maps
  integer(4), allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:), lm2l(:)
  integer(4), allocatable :: nidxlo(:), idxlostart(:), idxlomap(:,:)
  integer(4), allocatable :: idxgnt(:,:,:), nnz(:,:)
  real(8), allocatable :: listgnt(:,:,:)

  ! FFT grid 
  type(fftmap_type) :: fftmap
  ! \Theta(r)
  complex(8), allocatable :: zfftcf(:)
  ! FFT of evecks
  complex(8), allocatable :: fftevecsk(:,:,:)

  ! Radial integrals
  complex(8), allocatable, dimension(:,:,:,:) :: rigntaa, rigntal, rigntla, rigntll

  public :: emat_write_plain
  public :: emat_initialize, emat_gen, emat_finalize

  contains

    subroutine emat_initialize(lmaxapw_, lmaxexp_)
      integer(4), intent(in), optional :: lmaxapw_, lmaxexp_

      ! Check if k-space grids are initialized
      if( .not. ematgrids_initialized()) then
        write(*,*) "Error(m_ematqk::emat_initialize): Initialize k-space grids beforehand&
          & with ematgrids_init"
        call terminate
      end if

      ! Deallocate module variables
      call emat_module_free()

      ! Save l_max setting pertaining to the APWs
      ! Note: Used via input%groundstate%lmaxapw in 
      !       match routine. It is reset when
      !       emat_finalize is called.
      gs_lmaxapw_save = input%groundstate%lmaxapw
      lmaxapw = gs_lmaxapw_save
      if(present(lmaxapw_)) then
        if( lmaxapw_ >= 0 ) then
          lmaxapw = lmaxapw_
        end if
      end if
      input%groundstate%lmaxapw = lmaxapw

      ! Default l_max for expansion of exponential
      lmaxexp = 3
      if(present(lmaxexp_)) then
        if( lmaxexp_ >= 0 ) then
          lmaxexp = lmaxexp_
        end if
      else
        if(associated(input%xs)) then 
          lmaxexp = input%xs%lmaxemat
        end if
      end if

      ! Write info
      write(*,*) "Info(emat_initialize):"
      write(*,*) "  lmaxapw_orig =", gs_lmaxapw_save
      write(*,*) "  input%groundstate%lmaxmat = ", input%groundstate%lmaxmat
      write(*,*) "  lmaxapw =", lmaxapw
      write(*,*) "  lmaxexp =", lmaxexp

      ! Make combinded index maps for lmo and LO's 
      call emat_make_index_maps()

      ! Make list of nonzero gaunts
      call emat_make_gaunt()

      ! Note: Assume init2 to be called
      ! Generate radial functions
      ! stored in mod_apw_lo
      ! Note: apwfr are depenent on the effective muffin tin 
      !       potential stored in STATE.OUT, so needs a call to readstate
      !       e.g. from init2
      !       (entails compilcations if k and k+q grids are not the same)
      !call genapwfr
      !call genlofr

      ! Make fft grid and get realspcae
      ! representation of G-space characteristic function
      call emat_fft_init()

      ! Allocate radial integrals for
      ! one G+q 
      call emat_alloc_radial()

      ! Set module flag. 
      initialized = .true.

    end subroutine emat_initialize

    subroutine emat_finalize()

      if( .not. initialized) then
        write(*,*) "Warning(m_ematqk::emat_finalize): Module was not initialize."
        return
      end if

      ! Reset lmax values
      input%groundstate%lmaxapw = gs_lmaxapw_save

      ! Free index maps, fft grid and radial integrals
      call emat_module_free()

      ! Set module flag 
      initialized = .false.

      fftsaved = .false.

    end subroutine emat_finalize

    subroutine emat_module_free()

      ! Deallocate combinded index maps
      if(allocated(nlmo)) deallocate(nlmo)
      if(allocated(lmo2l)) deallocate(lmo2l)
      if(allocated(lmo2m)) deallocate(lmo2m)
      if(allocated(lmo2o)) deallocate(lmo2o)
      if(allocated(lm2l)) deallocate(lm2l)
      if(allocated(nidxlo)) deallocate(nidxlo)
      if(allocated(idxlostart)) deallocate(idxlostart)
      if(allocated(idxlomap)) deallocate(idxlomap)

      ! Deallocate gaunt list and index map
      if(allocated(idxgnt)) deallocate(idxgnt)
      if(allocated(listgnt)) deallocate(listgnt)
      if(allocated(nnz)) deallocate(nnz)

      ! Deallocate radial arrays
      !if(allocated(riaa)) deallocate(riaa)
      !if(allocated(rial)) deallocate(rial)
      !if(allocated(rill)) deallocate(rill)
      if(allocated(rigntaa)) deallocate(rigntaa)
      if(allocated(rigntal)) deallocate(rigntal)
      if(allocated(rigntla)) deallocate(rigntla)
      if(allocated(rigntll)) deallocate(rigntll)

      ! Deallocate fft grid and G-space characteristic fuctions
      if(associated(fftmap%igfft)) deallocate(fftmap%igfft)
      if(allocated(zfftcf)) deallocate(zfftcf)

      ! Deallocate saved FFTs
      if(allocated(fftevecsk)) deallocate(fftevecsk)

    end subroutine emat_module_free

    subroutine emat_make_index_maps()

      integer(4) :: is, l1, o1, m1, ia, ias, ilo1, idxlo1
      integer(4) :: lm1

      ! Count combined (l,m,o) indices and build index maps
      allocate(nlmo(nspecies))
      allocate(lmo2l((lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2m((lmaxapw + 1)**2*apwordmax, nspecies), &
                lmo2o((lmaxapw + 1)**2*apwordmax, nspecies))

      nlmomax = 0
      do is = 1, nspecies
        nlmo(is) = 0
        do l1 = 0, lmaxapw
          do o1 = 1, apword(l1, is)
            do m1 = -l1, l1
              nlmo(is) = nlmo(is) + 1
              lmo2l(nlmo(is), is) = l1
              lmo2m(nlmo(is), is) = m1
              lmo2o(nlmo(is), is) = o1
            end do
          end do
        end do
        nlmomax = max(nlmomax, nlmo(is))
      end do

      allocate(lm2l((lmaxexp + 1)**2))
      do l1 = 0, lmaxexp
        do m1 = -l1, l1
          lm1 = idxlm(l1, m1)
          lm2l(lm1) = l1
        end do
      end do

      ! Make LO-index map
      !   Number of LO's per species
      allocate(nidxlo(nspecies))
      do is = 1, nspecies
        nidxlo(is) = 0
        do ilo1 = 1, nlorb(is)
          l1 = lorbl(ilo1, is)
          nidxlo(is) = nidxlo(is) + 2*l1+1
        end do
      end do


      !   Start index of the LO's of atom ia of species is
      !   in the total indexlist of LO's
      allocate(idxlostart(natmtot))
      idxlostart(1) = 1
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          if(ias < natmtot) then
            idxlostart(ias+1) = idxlostart(ias)+nidxlo(is)
          end if
        end do
      end do

      !   Full LO index map 
      allocate(idxlomap(nlotot, 2))
      idxlo1 = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          do ilo1 = 1, nlorb(is)
            l1 = lorbl(ilo1, is)
            do m1 = -l1, l1
              lm1 = idxlm(l1, m1)
              idxlo1 = idxlo1 + 1
              idxlomap(idxlo1, :) = [ lm1, ilo1 ]
            end do
          end do
        end do
      end do
      
    end subroutine emat_make_index_maps

    subroutine emat_make_gaunt()

      integer(4) :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3, i
      real(8) :: gnt

      real(8), external :: gaunt

      ! Build non-zero Gaunt list
      allocate(idxgnt(lmaxapw + 1, (lmaxapw + 1)**2, (lmaxapw + 1)**2))
      allocate(listgnt(lmaxapw + 1, (lmaxapw + 1)**2, (lmaxapw + 1)**2))
      allocate(nnz((lmaxapw + 1)**2, (lmaxapw + 1)**2))

      idxgnt = 0
      listgnt = 0.d0
      nnz = 0

      do l3 = 0, lmaxapw
        do m3 = -l3, l3
          lm3 = idxlm(l3, m3)

          do l1 = 0, lmaxapw
            do m1 = -l1, l1
              lm1 = idxlm(l1, m1)

              i = 0
              do l2 = 0, lmaxexp
                do m2 = -l2, l2
                  lm2 = idxlm(l2, m2)

                  gnt = gaunt(l1, l2, l3, m1, m2, m3)

                  if(gnt .ne. 0.d0) then
                    i = i + 1
                    listgnt(i, lm1, lm3) = gnt
                    idxgnt(i, lm1, lm3) = lm2
                  end if

                end do
              end do

              ! Remember how many non zero gaunts exsist at (lm1,lm3)
              nnz(lm1, lm3) = i

            end do
          end do
        
        end do
      end do

    end subroutine emat_make_gaunt

    subroutine emat_fft_init()
      use mod_Gvector, only: cfunig
      real(8) :: emat_gmax

      integer(4) :: ig

      ! Setup fft G grid
      emat_gmax = 2*gkset%gkmax + gqset%gkmax
      call genfftmap(fftmap, emat_gmax)

      ! Inplace 3d inverse fourier transform of the
      ! G-space characteristic function 
      ! \Theta(G) --> \Theta(r) = { 0 inside MTs, 1 outside MTs}
      allocate(zfftcf(fftmap%ngrtot+1))
      zfftcf = zzero
      do ig = 1, gset%ngvec
        if(gset%gc(ig) .lt. emat_gmax) then
         zfftcf(fftmap%igfft(ig)) = cfunig(ig)
        end if
      end do
      call zfftifc(3, fftmap%ngrid, 1, zfftcf)

    end subroutine emat_fft_init

    subroutine emat_alloc_radial()
      
      ! Allocate radial integral arrays
      allocate(rigntaa(nlmomax, nlmomax, natmtot, gqset%ngkmax), &
                rigntal(nlmomax, nlotot, natmtot, gqset%ngkmax), &
                rigntla(nlotot, nlmomax, natmtot, gqset%ngkmax), &
                rigntll(nlotot, nlotot, natmtot, gqset%ngkmax))

      !allocate(riaa(0:lmaxexp, 0:lmaxapw, apwordmax, 0:lmaxapw, apwordmax), &
      !          rial(0:lmaxexp, 0:lmaxapw, apwordmax, nlomax), &
      !          rill(0:lmaxexp, nlomax, nlomax))

    end subroutine emat_alloc_radial

    subroutine emat_radial_init(iq)
      integer(4), intent(in) :: iq

      integer(4) :: igq, ngq
      integer(4) :: l1, l2, l3, m1, m2, o1, o2, lm1, lm2, lm3, is, ia, ias, i
      integer(4) :: lmo1, lmo2, ilo1, ilo2, idxlo1, idxlo2
      real(8), allocatable :: riaa(:,:,:,:,:), rial(:,:,:,:), rill(:,:,:)

      integer(4) :: ispin
      real(8), parameter :: epsdiff = 1.0d-13
      complex(8) :: strf
      complex(8), allocatable :: ylm(:)

      ! Only one spin component supported so far
      ispin = 1

      ! Get number of G+q vectors
      ngq = gqset%ngk(ispin, iq)


      ! Do for all G+q
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(strf, ylm, riaa, rial, rill), &
!$OMP& PRIVATE(i, igq, is, ia, ias, idxlo1, idxlo2), &
!$OMP& PRIVATE( lm3, l3, lmo2, l2, m2, o2, lm2, ilo2, ilo1, lmo1, l1, m1, o1, lm1)
#endif
      allocate(ylm((lmaxexp + 1)**2))
      allocate(riaa(0:lmaxexp, 0:lmaxapw, apwordmax, 0:lmaxapw, apwordmax), &
                rial(0:lmaxexp, 0:lmaxapw, apwordmax, nlomax), &
                rill(0:lmaxexp, nlomax, nlomax))
#ifdef USEOMP
!$OMP DO
#endif
      do igq = 1, ngq
        
        ! Get sph. harmonics up to l_max used in 
        ! expansion of exponential
        call genylm(lmaxexp, gqset%tpgkc(1:2,igq,ispin,iq), ylm)

        ! Radial integrals times Gaunt and expansion prefactor
        rigntaa(:,:,:,igq) = zzero
        rigntal(:,:,:,igq) = zzero
        rigntla(:,:,:,igq) = zzero
        rigntll(:,:,:,igq) = zzero

        !++++++++++++++++++++++++++++!
        ! Compute the radial integal !
        !++++++++++++++++++++++++++++!
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)

            ! e^{-i (q+G) \cdot R_\alpha}
            strf = conjg(gqset%sfacgk(igq, ias, ispin, iq))

            ! Generate radial integral
            call emat_genri(gqset%gkc(igq, ispin, iq), is, ia, riaa, rial, rill)

            ! APW-APW
            do lmo2 = 1, nlmo(is)
              l2 = lmo2l(lmo2, is)
              m2 = lmo2m(lmo2, is)
              o2 = lmo2o(lmo2, is)
              lm2 = idxlm(l2, m2)

              do lmo1 = 1, nlmo(is)
                l1 = lmo2l(lmo1, is)
                m1 = lmo2m(lmo1, is)
                o1 = lmo2o(lmo1, is)
                lm1 = idxlm(l1, m1)

                do i = 1, nnz(lm1, lm2) 
                  lm3 = idxgnt(i, lm1, lm2)
                  l3 = lm2l(lm3)

                  rigntaa(lmo1, lmo2, ias, igq) = rigntaa(lmo1, lmo2, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm1, lm2)*riaa(l3, l1, o1, l2, o2)
                end do

              end do
            end do

            ! APW-LO and LO-APW
            do lmo1 = 1, nlmo(is)
              l1 = lmo2l(lmo1, is)
              m1 = lmo2m(lmo1, is)
              o1 = lmo2o(lmo1, is)
              lm1 = idxlm(l1, m1)

              do idxlo1 = idxlostart(ias), idxlostart(ias)+nidxlo(is)-1
                ilo1 = idxlomap(idxlo1, 2)
                lm2  = idxlomap(idxlo1, 1)

                ! APW-LO
                do i = 1, nnz(lm1, lm2)
                  lm3 = idxgnt(i, lm1, lm2)
                  l3 = lm2l(lm3)

                  rigntal(lmo1, idxlo1, ias, igq) =&
                    & rigntal(lmo1, idxlo1, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm1, lm2)*rial(l3, l1, o1, ilo1)
                end do

                ! LO-APW
                do i = 1, nnz(lm1, lm2) 
                  lm3 = idxgnt(i, lm2, lm1)
                  l3 = lm2l(lm3)

                  rigntla(idxlo1, lmo1, ias, igq) =&
                    & rigntla(idxlo1, lmo1, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm2, lm1)*rial(l3, l1, o1, ilo1)
                end do

              end do
            end do

            ! LO-LO
            do idxlo2 = idxlostart(ias), idxlostart(ias)+nidxlo(is)-1
              ilo2 = idxlomap(idxlo2, 2)
              lm2  = idxlomap(idxlo2, 1)

              do idxlo1 = idxlostart(ias), idxlostart(ias)+nidxlo(is)-1
                ilo1 = idxlomap(idxlo1, 2)
                lm1  = idxlomap(idxlo1, 1)

                do i = 1, nnz(lm1, lm2) 
                  lm3 = idxgnt(i, lm1, lm2)
                  l3 = lm2l(lm3)

                  rigntll(idxlo1, idxlo2, ias, igq) =&
                    & rigntll(idxlo1, idxlo2, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm1, lm2)*rill(l3, ilo1, ilo2)
                end do

              end do
            end do
            
          end do ! atoms
        end do !species

      end do ! G+q
#ifdef USEOMP
!$OMP END DO  
#endif
      deallocate(riaa, rial, rill)
      deallocate(ylm)
#ifdef USEOMP
!$OMP END PARALLEL
#endif


      contains

        subroutine emat_genri(gqc, is, ia, riaa, rial, rill)
          integer(4), intent(in) :: is, ia
          real(8), intent(in) :: gqc
          real(8), intent(out) :: riaa(0:lmaxexp, 0:lmaxapw, apwordmax,&
                                       & 0:lmaxapw, apwordmax)
          real(8), intent(out) :: rial(0:lmaxexp, 0:lmaxapw, apwordmax, nlomax)
          real(8), intent(out) :: rill(0:lmaxexp, nlomax, nlomax)

          integer(4) :: ias, ir, nr, l1, l2, l3, o1, o2, ilo1, ilo2
          real(8) :: x

          real(8), allocatable :: jlqgr(:,:), fr(:), gf(:), cf(:,:)
          
          ! Get combined atom-species index
          ias = idxas(ia, is)

          ! Number of radial points per species
          nr = nrmt(is)

          ! Generate spherical Bessel functions
          allocate(jlqgr(0:lmaxexp, nr))
          do ir = 1, nr
            x = gqc*spr(ir, is)
            call sbessel(lmaxexp, x, jlqgr(:, ir))
          end do

          allocate(fr(nr), gf(nr), cf(3, nr))

          ! APW-APW
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(l1, l2, l3, o1, o2, ir, fr, gf, cf)
!!$OMP DO  
!#endif
          do l2 = 0, lmaxapw
            do o2 = 1, apword(l2, is)

              do l1 = 0, lmaxapw
                do o1 = 1, apword(l1, is)

                  do l3 = 0, lmaxexp

                    do ir = 1, nr
                      fr(ir) = apwfr(ir, 1, o1, l1, ias)*jlqgr(l3, ir)&
                             &* apwfr(ir, 1, o2, l2, ias)*spr(ir, is)**2
                    end do
                    call fderiv(-1, nr, spr(:, is), fr, gf, cf)
                    riaa(l3, l1, o1, l2, o2) = gf(nr)

                  end do

                end do
              end do

            end do
          end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif

          ! APW-LO
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ilo1, l1, l3, o1, ir, fr, gf, cf)
!!$OMP DO  
!#endif
          do ilo1 = 1, nlorb(is)
            do l1 = 0, lmaxapw
              do o1 = 1, apword(l1, is)

                do l3 = 0, lmaxexp

                  do ir = 1, nr
                    fr(ir) = apwfr(ir, 1, o1, l1, ias)*jlqgr(l3, ir)&
                           &* lofr(ir, 1, ilo1, ias)*spr(ir, is)**2
                  end do
                  call fderiv(-1, nr, spr(:, is), fr, gf, cf)
                  rial(l3, l1, o1, ilo1) = gf(nr)

                end do

              end do
            end do
          end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif

          ! LO-LO
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ilo1, ilo2, l3, ir, fr, gf, cf)
!!$OMP DO  
!#endif
          do ilo2 = 1, nlorb(is)
            do ilo1 = 1, nlorb(is)
              do l3 = 0, lmaxexp
                do ir = 1, nr
                  fr(ir) = lofr(ir, 1, ilo1, ias)*jlqgr(l3, ir)&
                         &* lofr(ir, 1, ilo2, ias)*spr(ir, is)**2
                end do
                call fderiv(-1, nr, spr(:, is), fr, gf, cf)
                rill(l3, ilo1, ilo2) = gf(nr)
              end do
            end do
          end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif

          deallocate(jlqgr, fr, gf, cf)
          return
        end subroutine emat_genri
      
    end subroutine emat_radial_init

    subroutine emat_buffer_ffts(evecks)
      complex(8), intent(in) :: evecks(:,:,:)

      integer(4) :: igp, ik, ist, nk, nstk, ispin

      ispin = 1
      nk = size(evecks,1)
      nstk = size(evecks,2)

      if(.not. initialized) then
        write(*,*) "Error(emat_save_ffts): Module not initialized"
        return
      end if

      ! 3d inverse FFT of evec(G,k,nstk) for all passed bands and k
      ! evec(G+k,i) --> evec(r,i)
      ! and then multiply by characteristic lattice function
      ! \Theta(r) and take the complex conjugate
      if(allocated(fftevecsk)) deallocate(fftevecsk)
      allocate(fftevecsk(fftmap%ngrtot+1, nstk, nk))
      fftevecsk = zzero
#ifdef USEOMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ik, ist, igp) COLLAPSE(2)
#endif
      do ik = 1, nk
        do ist = 1, nstk

          ! Map evec to FFT grid
          do igp = 1, gkset%ngk(ispin, ik) 
            fftevecsk(fftmap%igfft(gkset%igkig(igp, ispin, ik)), ist, ik) =&
              & evecks(igp, ist, ik)
          end do
          ! Do the inverse FFT
          call zfftifc(3, fftmap%ngrid, 1, fftevecsk(:, ist, ik))
          ! (evec(r,i)*\Theta(r))^*
          fftevecsk(:, ist, ik) = conjg(fftevecsk(:, ist, ik))*zfftcf(:) 

        end do
      end do
#ifdef USEOMP
!$OMP END PARALLEL DO
#endif
      fftsaved = .true.
    end subroutine emat_buffer_ffts

    subroutine emat_gen(ik, iq, eveck, eveckq, emat)
      use mod_Gkvector, only: ngkmax

      ! I/O
      integer(4), intent(in) :: ik, iq
      complex(8), intent(in) :: eveck(:,:), eveckq(:,:)
      complex(8), intent(out) :: emat(:,:,:)

      ! Local
      integer(4) :: ikq, ivgs(3), igs
      integer(4) :: nstk, nstkq, ngq, ngkmax_save
      integer(4) :: ngk, ngkq
      integer(4) :: ispin
      complex(8), allocatable :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)

      ! Spin index not yet supported
      ispin = 1
      
      ! Get corresponding kq grid point k+q = kq + G, with k,q and kq in [0,1)
      ikq = qset%ikiq2ikp_nr(ik,iq)
      igs = qset%ikiq2ig_nr(ik,iq)
      ivgs = gset%ivg(1:3,igs)

      ! Sanity checks
      if(size(eveck,1) /= size(eveckq,1)) then
        write(*,*) "Error(m_ematqk::emat_gen):&
          & Differing basis size (size(eveck,1) /= size(eveckq,1))"
        call terminate
      end if
      nstk = size(eveck,2)
      nstkq = size(eveckq,2)
      if(size(emat,1) /= nstk .or. size(emat,2) /= nstkq) then
        write(*,*) "Error(m_ematqk::emat_gen):&
          & emat array shape is not compatible with eigenvector arrays&
          & shape(emat) =", shape(emat), " shape(eveck) =", shape(eveck),&
          & " shape(eveckq) =", shape(eveckq)
        call terminate
      end if

      ngq = gqset%ngk(ispin, iq)

      if(size(emat,3) /= ngq) then
        write(*,*) "Error(m_ematqk::emat_gen):&
          & emat array shape is not compatible with ngq&
          & size(emat,3) =", size(emat,3), " ngq =", ngq
        call terminate
      end if

      emat = zzero

      ngk = gkset%ngk(ispin, ik)
      ngkq = gkqmtset%ngk(ispin, ikq)
      
      ! Find matching coefficients for k-point ik
      !   Note: Uses groundstate lmaxapw
      !   Note: Uses mod_Gkvector::ngkmax
      ngkmax_save = ngkmax
      ngkmax = gkset%ngkmax
      allocate( apwalmk( gkset%ngkmax, apwordmax, lmmaxapw, natmtot))
      call match(ngk, gkset%gkc(:, ispin, ik),&
       & gkset%tpgkc(1:2, :, ispin, ik), gkset%sfacgk(:,:,ispin, ik), apwalmk)

      ! Find the matching coefficients for k-point k+q
      ngkmax = gkqmtset%ngkmax
      allocate( apwalmkq( gkqmtset%ngkmax, apwordmax, lmmaxapw, natmtot))
      call match(ngkq, gkqmtset%gkc(:, ispin, ikq),&
       & gkqmtset%tpgkc(1:2, :, ispin, ikq), gkset%sfacgk(:,:,ispin, ikq), apwalmkq)

      !----------------------------------------------!
      !  Compute PW matrix elements for G+q < gqmax  !
      !  (or exactly G=0, if gqmax = 0 was set)      !
      !----------------------------------------------!
      !--------------------------------------!
      !     Muffin tin matrix elements       !
      !--------------------------------------!
      ! Dependes on G+q

      ! Generate radial integrals for G+q, where q is in [0,1) cell
      ! And store them in the module variables
      call emat_radial_init(iq)
      ! Build block matrix (for each G+q)
      call emat_mt_part(emat)
       
      !--------------------------------------!
      !     Interstitial matrix elements     !
      !--------------------------------------!
      ! Dependes on k, q, k+q, G_shift
      !call emat_inter_part(emat)
      call emat_inter_part(emat)

      deallocate(apwalmk, apwalmkq)

      contains 

        subroutine emat_mt_part(emat)

          complex(8), intent(inout) :: emat(:,:,:)

          integer(4) :: is, ia, ias
          integer(4) :: l, m, o, lmo, lm
          integer(4) :: igq, i

          ! Block matrix
          ! [_AA_|_AL_]
          ! [ LA | LL ]
          complex(8), allocatable, dimension(:,:) :: blockmt
          ! Helper matrices
          complex(8), allocatable, dimension(:,:) :: auxmat, auxmat2
          complex(8), allocatable, dimension(:,:) :: match_combined1, match_combined2
          complex(8), allocatable, dimension(:,:) :: aamat, almat, lamat


#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP& PRIVATE( blockmt, aamat, almat, lamat), &
!$OMP& PRIVATE( match_combined1, match_combined2), &
!$OMP& PRIVATE( auxmat, auxmat2), &
!$OMP& PRIVATE( is, ia, ias, igq), &
!$OMP& PRIVATE( lmo, l, m, o, lm)
#endif
          allocate(blockmt(ngk+nlotot,ngkq+nlotot))
          allocate(auxmat(nlmomax, ngkq))
          allocate(auxmat2(ngk+nlotot, nstkq))
          allocate(match_combined1(nlmomax, ngk))
          allocate(match_combined2(nlmomax, ngkq))
          allocate(aamat(ngk, ngkq), almat(ngk, nlotot), lamat(nlotot, ngkq))

#ifdef USEOMP
!$OMP DO 
#endif
          ! Do for all G+q vectors for current q
          Gqloop: do igq = 1, ngq

            ! Build block matrix for current G+q
            blockmt(:,:) = zzero

! If nias >> nigq
!#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, ia, ias), &
!!$OMP& PRIVATE( match_combined1, match_combined2), &
!!$OMP& PRIVATE( lmo, l, m, o, lm, auxmat, aamat, almat, lamat)
!#endif
!            allocate(auxmat(nlmomax, ngkq))
!            allocate(match_combined1(nlmomax, ngk))
!            allocate(match_combined2(nlmomax, ngkq))
!            allocate(aamat(ngk, ngkq), almat(ngk, nlotot), lamat(nlotot, ngkq))
!
!#ifdef USEOMP
!!$OMP DO COLLAPSE(2) REDUCTION(+:blockmt)
!#endif

            do is = 1, nspecies
              do ia = 1, natoms(is)
                ias = idxas(ia, is)

                ! Write combined matching coefficient matrices
                ! Note: zeroing unnecessary
                !match_combined1(:,:) = zzero
                !match_combined2(:,:) = zzero

                do lmo = 1, nlmo(is)
                  l = lmo2l(lmo, is)
                  m = lmo2m(lmo, is)
                  o = lmo2o(lmo, is)
                  lm = idxlm(l, m)
                  match_combined1(lmo, :) = apwalmk(1:ngk, o, lm, ias)
                  match_combined2(lmo, :) = apwalmkq(1:ngkq, o, lm, ias)
                end do

                ! Sum up block matrix

                ! APW-APW
                ! Note: zeroing unnecessary
                !auxmat(:,:) = zzero
                call ZGEMM('N', 'N', nlmo(is), ngkq, nlmo(is), zone, &
                    rigntaa(1:nlmo(is), 1:nlmo(is), ias, igq), nlmo(is), &
                    match_combined2(1:nlmo(is), :), nlmo(is), zzero, &
                    auxmat(1:nlmo(is), :), nlmo(is))
                call ZGEMM('C', 'N', ngk, ngkq, nlmo(is), zone, &
                    match_combined1(1:nlmo(is), :), nlmo(is), &
                    auxmat(1:nlmo(is), :), nlmo(is), zzero, &
                    aamat, ngk)

                blockmt(1:ngk, 1:ngkq) = blockmt(1:ngk, 1:ngkq) + aamat(:,:)

                ! APW-LO
                call ZGEMM('C', 'N', ngk, nlotot, nlmo(is), zone, &
                    match_combined1(1:nlmo(is), :), nlmo(is), &
                    rigntal(1:nlmo(is), :, ias, igq), nlmo(is), zzero, &
                    almat, ngk)

                blockmt(1:ngk, (ngkq+1):(ngkq+nlotot)) =&
                  & blockmt(1:ngk, (ngkq+1):(ngkq+nlotot)) + almat(:,:)

                ! LO-APW
                call ZGEMM('N','N', nlotot, ngkq, nlmo(is), zone, &
                    rigntla(:, 1:nlmo(is), ias, igq), nlotot, &
                    match_combined2(1:nlmo(is), :), nlmo(is), zzero, &
                    lamat, nlotot)

                blockmt((ngk+1):(ngk+nlotot), 1:ngkq) =&
                  & blockmt((ngk+1):(ngk+nlotot), 1:ngkq) + lamat(:,:)

                ! LO-LO
                blockmt((ngk+1):(ngk+nlotot), (ngkq+1):(ngkq+nlotot)) =&
                  & blockmt((ngk+1):(ngk+nlotot), (ngkq+1):(ngkq+nlotot))&
                  &+ rigntll(:, :, ias, igq)

              end do
            end do
! IF nias >> nigq
!#ifdef USEOMP
!!$OMP END DO
!#endif
!            deallocate(auxmat)
!            deallocate(match_combined1)
!            deallocate(match_combined2)
!            deallocate(aamat, almat, lamat)

!#ifdef USEOMP
!!$OMP END PARALLEL
!#endif

            ! Compute final total muffin tin contribution
            call ZGEMM('N', 'N', ngk+nlotot, nstkq, ngkq+nlotot, zone, &
                blockmt(:,:), ngk+nlotot, &
                eveckq(1:(ngkq+nlotot), :), ngkq+nlotot, zzero, &
                auxmat2(:,:), ngk+nlotot)
            call ZGEMM('C', 'N', nstk, nstkq, ngk+nlotot, cmplx(fourpi, 0.d0, 8), &
                eveck(1:(ngk+nlotot), :), ngk+nlotot, &
                auxmat2(:,:), ngk+nlotot, zzero, &
                emat(:,:, igq), nstk)

          end do Gqloop
#ifdef USEOMP
!$OMP END DO
#endif
          deallocate(auxmat, auxmat2)
          deallocate(match_combined1)
          deallocate(match_combined2)
          deallocate(aamat, almat, lamat)
          deallocate(blockmt)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
      
        end subroutine emat_mt_part

        subroutine emat_inter_part(emat)
          complex(8), intent(inout) :: emat(:,:,:)

          integer(4) :: igp
          integer(4) :: igs, ivg(3)
          integer(4) :: ist1, ist2
          complex(8), allocatable :: zfftk(:,:), zfftkq(:), zfftres(:)

          if( .not. fftsaved) then
            ! 3d inverse FFT of evec(G+k,nstk) for all passed bands
            ! evec(G+k,i) --> evec(r,i)
            ! and then multiply by characteristic lattice function
            ! \Theta(r) and take the complex conjugate
            allocate(zfftk(fftmap%ngrtot+1, nstk))
            zfftk = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,igp)
!$OMP DO
#endif
            do ist1 = 1, nstk
              ! Map evec to FFT grid
              do igp = 1, ngk 
                zfftk(fftmap%igfft(gkset%igkig(igp, ispin, ik)), ist1) = eveck(igp, ist1)
              end do
              ! Do the inverse FFT
              call zfftifc(3, fftmap%ngrid, 1, zfftk(:, ist1))
              ! (evec(r,i)*\Theta(r))^*
              zfftk(:, ist1) = conjg(zfftk(:, ist1))*zfftcf(:) 
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
          end if
               
          ! Generate the 3d fourier transform of the real space integrand 
          ! zfftres = eveck^*_ist1(r)*\Theta(r)*eveckq_ist2(r)
          ! for all ist1/2 combinations and add this to emat.
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,ist2,igp,igs,ivg,zfftkq,zfftres)
#endif
          allocate(zfftres(fftmap%ngrtot+1))
          allocate(zfftkq(fftmap%ngrtot+1))
#ifdef USEOMP
!$OMP DO
#endif
          do ist2 = 1, nstkq

            ! 3d inverse FFT of evec(G+k+q,i) for current band index
            ! evec(G+k+q,i) --> evec(r,i)
            zfftkq = zzero
            ! Map to FFT grid, respect G shift if k+q has lattice component
            if( any(ivgs /= 0) ) then
              do igp = 1, ngkq 
                ivg = gset%ivg(1:3, gkqmtset%igkig(igp,ispin,ikq)) - ivgs
                igs = gset%ivgig( ivg(1), ivg(2), ivg(3))
                zfftkq(fftmap%igfft(igs)) = eveckq(igp, ist2)
              end do
            else
              do igp = 1, ngkq 
                zfftkq(fftmap%igfft(gkqmtset%igkig(igp,ispin,ikq))) = eveckq(igp, ist2)
              end do
            end if
            call zfftifc(3, fftmap%ngrid, 1, zfftkq)

            do ist1 = 1, nstk
              ! Generate real space integrand
              ! eveck^*_ist1(r)*\Theta(r)*eveckq_ist2(r)
              ! and do a 3d fft to G space
              do igp = 1, fftmap%ngrtot
                if(.not. fftsaved) then 
                  zfftres(igp) = zfftkq(igp)*zfftk(igp, ist1) 
                else
                  zfftres(igp) = zfftkq(igp)*fftevecsk(igp, ist1, ik) 
                end if
              end do
              zfftres(fftmap%ngrtot+1) = zzero
              call zfftifc(3, fftmap%ngrid, -1, zfftres)

              ! Add result to emat
              do igp = 1, ngq
                emat(ist1, ist2, igp) = emat(ist1, ist2, igp)&
                  &+ zfftres(fftmap%igfft(gqset%igkig(igp, ispin, iq)))
              end do

            end do

          end do
#ifdef USEOMP
!$OMP END DO
#endif
          deallocate(zfftres, zfftkq)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
          deallocate(zfftk)
        end subroutine emat_inter_part

    end subroutine emat_gen

    subroutine emat_write_plain(fname, emat, igq)
      use m_writecmplxparts
      character(*), intent(in) :: fname
      complex(8), intent(in) :: emat(:,:,:)
      integer(4), intent(in), optional :: igq

      integer(4) :: i

      if(present(igq)) then 
        call writecmplxparts(trim(fname), remat=dble(emat(:,:,igq)),&
          & immat=aimag(emat(:,:,igq))) 
      else
        do i = 1, size(emat,3)
          call writecmplxparts(trim(fname), remat=dble(emat(:,:,i)),&
            & immat=aimag(emat(:,:,i)), ik1=i) 
        end do
      end if

    end subroutine emat_write_plain

end module m_ematqk
