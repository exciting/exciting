module m_ematqk
  use modmpi
  use modinput
  use mod_constants
  use mod_atoms, only: nspecies, natmtot, natoms, idxas
  use mod_APW_LO, only: apwordmax, apword, nlorb, lorbl, nlotot, nlomax, lolmax
  use mod_muffin_tin

  !use mod_ematgrids
  use mod_xsgrids, only: g, g_k, g_qq, xsgrids_initialized

  use mod_eigenvalue_occupancy
  use modxs, only : fftmap_type

  implicit none
  private

  logical :: initialized = .false.
  logical :: fftsaved = .false.

  integer(4) :: iqprev

  integer(4) :: nlmomax, nidxlomax, lmaxapw, lmaxexp, lmaxapwlo, lmaxmat, lmaxmatlo
  integer(4) :: lmaxapw_save

  integer(4) :: maxl_apwlo, maxl_mat, maxl_e

  ! Index maps
  integer(4), allocatable :: nlmo(:), lmo2l(:,:), lmo2m(:,:), lmo2o(:,:), lm2l(:)
  integer(4), allocatable :: nidxlo(:), idxlostart(:), idxlomap(:,:), myidxlm(:,:)
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
  real(8), allocatable :: riaa(:,:,:,:,:,:,:), rial(:,:,:,:,:,:), rill(:,:,:,:,:)

  ! Module flags
  logical :: precompute_radial, precompute_ffts, samestates

  public :: emat_write_plain
  public :: emat_initialize, emat_gen, emat_finalize, emat_precal_fft

  contains

    subroutine emat_initialize(lmaxapw_, lmaxexp_, pre_rad_, pre_fft_)
      integer(4), intent(in), optional :: lmaxapw_, lmaxexp_
      logical, intent(in), optional :: pre_rad_, pre_fft_

      ! Check if k-space grids are initialized
      if( .not. xsgrids_initialized()) then
        write(*,*) "Error(m_ematqk::emat_initialize):&
          & Initialize k-space grids beforehand&
          & with xsgrids_init"
        call terminate
      end if

      ! Deallocate module variables
      call emat_module_free()

      ! lmaxapw: Maximal l value considered to be included in 
      ! the LAPW part of the basis.
      ! Note: input%groundstate%lmaxapw defaults to 8, BUT xs%lmaxapw defaults to 10 
      !       in xsinit input%groundstate%lmaxapw is overwritten with the
      !       input%xs%lmaxapw value.
      ! Note: Used via input%groundstate%lmaxapw in 
      !       match routine.
      lmaxapw_save = input%groundstate%lmaxapw
      lmaxapw = lmaxapw_save
      if(present(lmaxapw_)) then
        if( lmaxapw_ >= 0 ) then
          lmaxapw = lmaxapw_
        end if
      end if
      ! Overwrite input variable to be used in match routine.
      ! Will be reset with finalization routine.
      input%groundstate%lmaxapw = lmaxapw

      ! Maximum l over LAPW and LO part of the basis 
      lmaxapwlo = max(lmaxapw, lolmax)

      ! lmaxmat: Maximal l value considered in the ground state Hamiltonian.
      ! Note: input%groundstate%lmaxmat defaluts to 8,
      !       so it is by defalut the same as input%groundstate%lmaxapw.
      !       In the XS case input%groundstate%lmaxmat is overwritten by
      !       input%xs%lmaxmat, which defaults to 5.
      !       There is also input%xs%lmaxapwwf which defaults to input%xs%lmaxmat.
      !       lmaxapwwf is used in the old implementation of ematqk.
      lmaxmat = input%groundstate%lmaxapw
      ! Maximum l over LAPW and LO part of the basis used in the GS Hamiltoninan.
      lmaxmatlo = max(lmaxmat, lolmax)

      ! lmaxexp: Maximal l value included in the Reighley expansion
      ! of the exponential e^{-i(G+q)r}.
      ! Note: input%xs%lmaxemat defaluts to 3.
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
      write(*,*) "  lmaxapw_orig =", lmaxapw_save
      write(*,*) "  lmaxapw =", lmaxapw
      write(*,*) "  lolmax =", lolmax
      write(*,*) "  lmaxapwlo =", lmaxapwlo
      write(*,*) "  lmaxmat =", lmaxmat
      write(*,*) "  lmaxmatlo =", lmaxmatlo
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
      !call genapwfr
      !call genlofr

      ! Make fft grid and get realspcae
      ! representation of G-space characteristic function
      call emat_fft_init()

      ! Allocate radial integrals
      call emat_alloc_radial()

      ! Precalculate radial integrals for reduced q-grid
      if(present(pre_rad_)) then 
        precompute_radial = pre_rad_
      else
        precompute_radial = .false.
      end if

      ! Set iqprev to unset
      iqprev = -1

      if(precompute_radial) then 
        call emat_precal_rad()
      end if

      ! Precalculate fft's of eigenvectors
      if(present(pre_fft_)) then
        precompute_ffts = pre_fft_
      else
        precompute_ffts = .false.
      end if

      ! Set module flag. 
      initialized = .true.

    end subroutine emat_initialize

    subroutine emat_finalize()

      if( .not. initialized) then
        write(*,*) "Warning(m_ematqk::emat_finalize): Module was not initialize."
        return
      end if

      ! Reset lmax values
      input%groundstate%lmaxapw = lmaxapw_save

      ! Free index maps, fft grid and radial integrals
      call emat_module_free()

      ! Set module flag 
      initialized = .false.

      precompute_radial = .false.

      ! Set iqprev to unset
      iqprev = -1

    end subroutine emat_finalize

    subroutine emat_module_free()

      ! Deallocate combinded index maps
      if(allocated(myidxlm)) deallocate(myidxlm)
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
      if(allocated(riaa)) deallocate(riaa)
      if(allocated(rial)) deallocate(rial)
      if(allocated(rill)) deallocate(rill)
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
      use mod_atoms, only: nspecies, natmtot, natoms, idxas
      use mod_APW_LO, only: apwordmax, apword, nlorb, lorbl, nlotot

      integer(4) :: is, l1, o1, m1, ia, ias, ilo1, idxlo1
      integer(4) :: lm1

      ! Count combined (l,m,o) indices and build index maps
      ! Maps for LAPW part of the basis form combined
      ! l,m,o index to the individual one l, m and o.
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
              ! Count number of LAPW basis functions for each species
              nlmo(is) = nlmo(is) + 1
              ! lmo --> l
              lmo2l(nlmo(is), is) = l1
              ! lmo --> m
              lmo2m(nlmo(is), is) = m1
              ! lmo --> o
              lmo2o(nlmo(is), is) = o1
            end do
          end do
        end do
        ! Maximum number of LAPW basis functions over all species
        nlmomax = max(nlmomax, nlmo(is))
      end do

      ! Build lm combinded index.
      ! Note: idxlm is computed in init0 using the value of
      !       input%groundstate%lmaxapw
      allocate(myidxlm(0:lmaxapwlo,-lmaxapwlo:lmaxapwlo))
      lm1=0
      do l1=0, lmaxapwlo
        do m1= -l1,l1
          lm1 = lm1+1
          myidxlm(l1, m1) = lm1
        end do
      end do

      ! Mapping for Reighley expanded exponential
      ! lm --> l
      allocate(lm2l((lmaxexp + 1)**2))
      do l1 = 0, lmaxexp
        do m1 = -l1, l1
          lm1 = myidxlm(l1, m1)
          lm2l(lm1) = l1
        end do
      end do

      ! Make LO-index map
      !   Count LO's
      nidxlomax = 0 
      allocate(nidxlo(nspecies))
      do is = 1, nspecies
        nidxlo(is) = 0
        do ilo1 = 1, nlorb(is)
          l1 = lorbl(ilo1, is)
          ! Number of LO's per species
          nidxlo(is) = nidxlo(is) + 2*l1+1
        end do
        ! Maximum number of LO's over all species
        nidxlomax = max(nidxlomax, nidxlo(is))
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
              lm1 = myidxlm(l1, m1)
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
      allocate(idxgnt(lmaxapwlo + 1, (lmaxapwlo + 1)**2, (lmaxapwlo + 1)**2))
      allocate(listgnt(lmaxapwlo + 1, (lmaxapwlo + 1)**2, (lmaxapwlo + 1)**2))
      allocate(nnz((lmaxapwlo + 1)**2, (lmaxapwlo + 1)**2))

      idxgnt = 0
      listgnt = 0.d0
      nnz = 0

      do l2 = 0, lmaxapwlo
        do m2 = -l2, l2
          lm2 = myidxlm(l2, m2)

          do l1 = 0, lmaxapwlo
            do m1 = -l1, l1
              lm1 = myidxlm(l1, m1)

              i = 0
              do l3 = 0, lmaxexp
                do m3 = -l3, l3
                  lm3 = myidxlm(l3, m3)

                  gnt = gaunt(l1, l3, l2, m1, m3, m2)

                  if(gnt .ne. 0.d0) then
                    i = i + 1
                    listgnt(i, lm1, lm2) = gnt
                    idxgnt(i, lm1, lm2) = lm3
                  end if

                end do
              end do

              ! Remember how many non zero gaunts exsist at (lm1,lm2)
              nnz(lm1, lm2) = i

            end do
          end do
        
        end do
      end do

    end subroutine emat_make_gaunt

    subroutine emat_fft_init()
      use mod_Gvector, only: cfunig
      use mod_xsgrids, only: g, g_k, g_qq
      real(8) :: emat_gmax

      integer(4) :: ig

      ! Setup fft G grid
      emat_gmax = 2*g_k%gkmax + g_qq%gkmax
      call genfftmap(fftmap, emat_gmax)

      ! Inplace 3d inverse fourier transform of the
      ! G-space characteristic function 
      ! \Theta(G) --> \Theta(r) = { 0 inside MTs, 1 outside MTs}
      allocate(zfftcf(fftmap%ngrtot+1))
      zfftcf = zzero
      do ig = 1, g%ngvec
        if(g%gc(ig) .lt. emat_gmax) then
         zfftcf(fftmap%igfft(ig)) = cfunig(ig)
        end if
      end do
      call zfftifc(3, fftmap%ngrid, 1, zfftcf)

    end subroutine emat_fft_init

    subroutine emat_alloc_radial()
      use mod_atoms, only: natmtot
      use mod_APW_LO, only: apwordmax, nlomax
      use mod_xsgrids, only: g_qq

      integer(4) :: ngqmax

      ngqmax = g_qq%ngknrmax

      ! Allocate radial integral arrays
      allocate(riaa(0:lmaxexp, 0:lmaxapw, apwordmax, 0:lmaxapw, apwordmax, natmtot, ngqmax), &
                rial(0:lmaxexp, 0:lmaxapw, apwordmax, nlomax, natmtot, ngqmax), &
                rill(0:lmaxexp, nlomax, nlomax, natmtot, ngqmax))
      
      ! Allocate radial integral times gaunt coefficients and bessel arrays
      allocate(rigntaa(nlmomax, nlmomax, natmtot, ngqmax), &
                rigntal(nlmomax, nidxlomax, natmtot, ngqmax), &
                rigntla(nidxlomax, nlmomax, natmtot, ngqmax), &
                rigntll(nidxlomax, nidxlomax, natmtot, ngqmax))

    end subroutine emat_alloc_radial

    subroutine emat_precal_rad()
      use m_genfilname
      use m_getunit
      use mod_xsgrids, only: q_q, g_qq
      use mod_atoms, only: nspecies, natoms, idxas

      character(256) :: fname
      integer(4) :: un
      integer(4) :: iq, iqnr, igq, is, ia, ias, ispin, ngq

      ispin = 1

      call system('[[ ! -e EMAT ]] && mkdir EMAT')

      ! Loop over reduced q-points (if q-set was reduced)
      do iq = 1, q_q%qset%nkpt

        ! Get eqvialent index of non-reduced set
        iqnr = q_q%qset%ikp2ik(iq)

        call genfilname(basename='EMAT/EMAT_RAD', iq=iq, filnam=fname)

        ! Get number of G+q vectors
        ngq = g_qq%ngknr(ispin, iq)

#ifdef USEOMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(igq, is, ia, ias)
#endif
        do igq = 1, ngq

          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)

              call emat_genri(g_qq%gknrc(igq, ispin, iq), is, ia,&
                & riaa(:,:,:,:,:,ias,igq), rial(:,:,:,:,ias,igq), rill(:,:,:,ias,igq))

            end do
          end do

        end do
#ifdef USEOMP
!$OMP END PARALLEL DO
#endif

        ! Save to file
        call getunit(un)
        open(un, file=trim(fname), form='unformatted', action='write', status='replace')
        write(un) riaa, rial, rill
        close(un)

      end do

    end subroutine emat_precal_rad

    subroutine emat_fetch_rad(iqnr)
      use m_genfilname
      use m_getunit
      use mod_symmetry, only: maxsymcrys
      use mod_xsgrids, only: q_q, g_qq
      use modinput

      integer(4), intent(in) :: iqnr

      character(256) :: fname
      integer(4) :: un
      integer(4) :: iq, numgq, ispin

      real(8) :: vqnr(3), vq(3)
      integer(4) :: jsym, jsymi
      integer(4) :: nsc, ivgsym(3)
      integer(4) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
      integer(4), allocatable :: igqmap(:)

      ! Fetch only if not already in module variable
      if(iqnr == iqprev) return

      ! Only one spin component supported so far
      ispin = 1

      ! Get equivalent reduced q-point
      iq = q_q%qset%ik2ikp(iqnr)
      numgq = g_qq%ngknr(ispin, iqnr)

      ! Read radial integals for reduced q-point from file
      call genfilname(basename='EMAT/EMAT_RAD', iq=iq, filnam=fname)
      call getunit(un)
      open(un, file=trim(fname), form='unformatted', action='read', status='old')
      read(un) riaa, rial, rill
      close(un)

      ! Get lattice coordiates of q vectors
      vq = q_q%qset%vkl(1:3,iq)
      vqnr = q_q%qset%vklnr(1:3,iqnr)

      ! Retreive radial integrals for non reduced q-point form reduced q-point integrals
      ! by remapping the G+q vector indices
      if(q_q%qset%isreduced .and. any(vq /= vqnr)) then 

        ! Find symmetry operations that map the reduced q-point to the non reduced one
        call findsymeqiv(input%xs%bse%fbzq, vqnr, vq, nsc, sc, ivgsc)

        ! Find the map that rotates the G+q-vectors
        allocate(igqmap(numgq))
        call b_findgqmap(iqnr, iq, nsc, sc, ivgsc, numgq, jsym, jsymi, ivgsym, igqmap)

        ! Rotate radial integrals using the remapped G+q indices
        riaa(:,:,:,:,:,:,1:numgq) = riaa(:,:,:,:,:,:,igqmap)
        rial(:,:,:,:,:,1:numgq) = rial(:,:,:,:,:,igqmap)
        rill(:,:,:,:,1:numgq) = rill(:,:,:,:,igqmap)

        deallocate(igqmap)

      end if

      iqprev = iqnr

    end subroutine emat_fetch_rad

    subroutine emat_precal_fft(ik, eveck, eveckq_)
      use m_genfilname
      use m_getunit
      use mod_xsgrids, only: g_k, g_kqmtp

      integer(4), intent(in) :: ik
      complex(8), intent(in) :: eveck(:,:)
      complex(8), intent(in), optional :: eveckq_(:,:)

      character(256) :: fname
      integer(4) :: un, nstk, nstkq, ist1, igp, ngk, ngkq, ispin
      complex(8), allocatable :: zfftk(:,:)

      ispin = 1

      call system('[[ ! -e EMAT ]] && mkdir EMAT')

      call genfilname(basename='EMAT/EMAT_FFT_EVECK', iq=ik, filnam=fname)

      nstk = size(eveck,2)
      ngk = g_k%ngk(ispin, ik)

      allocate(zfftk(fftmap%ngrtot+1, nstk))
      zfftk = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,igp)
!$OMP DO
#endif
      do ist1 = 1, nstk
        ! Map evec to FFT grid
        do igp = 1, ngk 
          zfftk(fftmap%igfft(g_k%igkig(igp, ispin, ik)), ist1) = eveck(igp, ist1)
        end do
        ! Do the inverse FFT
        call zfftifc(3, fftmap%ngrid, 1, zfftk(:, ist1))
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      ! Save to file
      call getunit(un)
      open(un, file=trim(fname), form='unformatted', action='write', status='replace')
      write(un) zfftk
      close(un)
      deallocate(zfftk)

      if(present(eveckq_)) then 
        nstkq = size(eveckq_,2)
        ngkq = g_kqmtp%ngk(ispin, ik)

        call genfilname(basename='EMAT/EMAT_FFT_EVECKQ', iq=ik, filnam=fname)
        write(*,*) "Filename=", trim(fname)

        allocate(zfftk(fftmap%ngrtot+1, nstkq))
        zfftk = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,igp)
!$OMP DO
#endif
        do ist1 = 1, nstkq
          ! Map evec to FFT grid
          do igp = 1, ngkq 
            zfftk(fftmap%igfft(g_kqmtp%igkig(igp, ispin, ik)), ist1) = eveckq_(igp, ist1)
          end do
          ! Do the inverse FFT
          call zfftifc(3, fftmap%ngrid, 1, zfftk(:, ist1))
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        ! Save to file
        call getunit(un)
        open(un, file=trim(fname), form='unformatted', action='write', status='replace')
        write(un) zfftk
        close(un)
        deallocate(zfftk)

      else

        samestates = .true.

      end if

    end subroutine emat_precal_fft

    subroutine emat_fetch_fft(zfft, ik_, ikq_)
      use m_genfilname
      use m_getunit

      complex(8), intent(out) :: zfft(:,:)
      integer(4), intent(in), optional :: ik_, ikq_

      character(256) :: fname
      integer(4) :: un

      if( .not. present(ik_) .and. .not. present(ikq_)) then
        write(*,*) "Specifiy ik or ikq"
        call terminate
      end if
      if( present(ik_) .and. present(ikq_)) then
        write(*,*) "Specifiy ik or ikq"
        call terminate
      end if
      
      if(present(ik_)) then 
        call genfilname(basename='EMAT/EMAT_FFT_EVECK', iq=ik_, filnam=fname)
        !write(*,*) "Getting FFT of eveck @ ik=", ik_, " Filename=", trim(fname)
      else
        if(samestates) then 
          call genfilname(basename='EMAT/EMAT_FFT_EVECK', iq=ikq_, filnam=fname)
        else
          call genfilname(basename='EMAT/EMAT_FFT_EVECKQ', iq=ikq_, filnam=fname)
        end if
        !write(*,*) "Getting FFT of eveckq @ ikq=", ikq_, " Filename=", trim(fname)
      end if

      call getunit(un)
      open(un, file=trim(fname), form='unformatted', action='read', status='old')
      read(un) zfft
      close(un)

    end subroutine emat_fetch_fft

    subroutine emat_radial_init(iq)
      use mod_atoms, only: nspecies, natoms, idxas
      use mod_xsgrids, only: g_qq

      integer(4), intent(in) :: iq

      integer(4) :: igq, ngq
      integer(4) :: l1, l2, l3, m1, m2, o1, o2, lm1, lm2, lm3, is, ia, ias, i
      integer(4) :: lmo1, lmo2, ilo1, ilo2, idxlo1, idxlo2

      integer(4) :: ispin
      real(8), parameter :: epsdiff = 1.0d-13
      complex(8) :: strf
      complex(8), allocatable :: ylm(:)

      ! Compute if not already in module variable
      if(iq == iqprev) return

      ! Only one spin component supported so far
      ispin = 1

      ! Get number of G+q vectors
      ngq = g_qq%ngknr(ispin, iq)

      ! Read radial integrals from file (riaa,rial,rill)
      if( precompute_radial ) then 
        call emat_fetch_rad(iq)
      end if

      ! Do for all G+q
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(strf, ylm), &
!$OMP& PRIVATE(i, igq, is, ia, ias, idxlo1, idxlo2), &
!$OMP& PRIVATE( lm3, l3, lmo2, l2, m2, o2, lm2, ilo2, ilo1, lmo1, l1, m1, o1, lm1)
#endif
      allocate(ylm((lmaxexp + 1)**2))
#ifdef USEOMP
!$OMP DO
#endif
      do igq = 1, ngq
        
        ! Get sph. harmonics up to l_max used in 
        ! expansion of exponential
        call genylm(lmaxexp, g_qq%tpgknrc(1:2,igq,ispin,iq), ylm)

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
            strf = conjg(g_qq%sfacgknr(igq, ias, ispin, iq))

            ! Generate radial integral
            if( .not. precompute_radial) then 
              call emat_genri(g_qq%gknrc(igq, ispin, iq), is, ia,&
                & riaa(:,:,:,:,:,ias,igq), rial(:,:,:,:,ias,igq), rill(:,:,:,ias,igq))
            end if

            ! APW-APW
            do lmo2 = 1, nlmo(is)
              l2 = lmo2l(lmo2, is)
              m2 = lmo2m(lmo2, is)
              o2 = lmo2o(lmo2, is)
              lm2 = myidxlm(l2, m2)

              do lmo1 = 1, nlmo(is)
                l1 = lmo2l(lmo1, is)
                m1 = lmo2m(lmo1, is)
                o1 = lmo2o(lmo1, is)
                lm1 = myidxlm(l1, m1)

                do i = 1, nnz(lm1, lm2) 
                  lm3 = idxgnt(i, lm1, lm2)
                  l3 = lm2l(lm3)

                  rigntaa(lmo1, lmo2, ias, igq) = rigntaa(lmo1, lmo2, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm1, lm2)*riaa(l3, l1, o1, l2, o2, ias, igq)
                end do

              end do
            end do

            ! APW-LO and LO-APW

            do lmo1 = 1, nlmo(is)
              l1 = lmo2l(lmo1, is)
              m1 = lmo2m(lmo1, is)
              o1 = lmo2o(lmo1, is)
              lm1 = myidxlm(l1, m1)

              do idxlo1 = 1, nidxlo(is)
                ilo1 = idxlomap(idxlo1+idxlostart(ias)-1, 2)
                lm2  = idxlomap(idxlo1+idxlostart(ias)-1, 1)

                ! APW-LO
                do i = 1, nnz(lm1, lm2)
                  lm3 = idxgnt(i, lm1, lm2)
                  l3 = lm2l(lm3)

                  rigntal(lmo1, idxlo1, ias, igq) =&
                    & rigntal(lmo1, idxlo1, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm1, lm2)*rial(l3, l1, o1, ilo1, ias, igq)
                end do

                ! LO-APW
                do i = 1, nnz(lm1, lm2) 
                  lm3 = idxgnt(i, lm2, lm1)
                  l3 = lm2l(lm3)

                  rigntla(idxlo1, lmo1, ias, igq) =&
                    & rigntla(idxlo1, lmo1, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm2, lm1)*rial(l3, l1, o1, ilo1, ias, igq)
                end do

              end do
            end do

            ! LO-LO
            do idxlo2 = 1, nidxlo(is)
              ilo2 = idxlomap(idxlo2+idxlostart(ias)-1, 2)
              lm2  = idxlomap(idxlo2+idxlostart(ias)-1, 1)

              do idxlo1 = 1, nidxlo(is)
                ilo1 = idxlomap(idxlo1+idxlostart(ias)-1, 2)
                lm1  = idxlomap(idxlo1+idxlostart(ias)-1, 1)

                do i = 1, nnz(lm1, lm2) 
                  lm3 = idxgnt(i, lm1, lm2)
                  l3 = lm2l(lm3)

                  rigntll(idxlo1, idxlo2, ias, igq) =&
                    & rigntll(idxlo1, idxlo2, ias, igq)&
                    &+ strf*conjg(zil(l3))*conjg(ylm(lm3))&
                    &* listgnt(i, lm1, lm2)*rill(l3, ilo1, ilo2, ias, igq)
                end do

              end do
            end do
            
          end do ! atoms
        end do !species

      end do ! G+q
#ifdef USEOMP
!$OMP END DO  
#endif
      deallocate(ylm)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
      
    end subroutine emat_radial_init

    subroutine emat_genri(gqc, is, ia, riaa, rial, rill)
      use mod_APW_LO, only: apwordmax, nlomax, apwfr, lofr
      use mod_atoms, only: idxas, spr

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
          do igp = 1, g_k%ngk(ispin, ik) 
            fftevecsk(fftmap%igfft(g_k%igkig(igp, ispin, ik)), ist, ik) =&
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
      use mod_xsgrids, only: q_q, g_kqmtp

      ! I/O
      integer(4), intent(in) :: ik, iq
      complex(8), intent(in) :: eveck(:,:), eveckq(:,:)
      complex(8), intent(out) :: emat(:,:,:)

      ! Local
      integer(4) :: ikq
      integer(4) :: nstk, nstkq, ngq, ngkmax_save, lmmaxapw_save, nbasis
      integer(4) :: ngk, ngkq
      integer(4) :: ispin
      complex(8), allocatable :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)

      ! Spin index not yet supported
      ispin = 1
      
      ! Get corresponding kq grid point k+q = kq + G, with k,q and kq in [0,1)
      ikq = q_q%ikiq2ikp_nr(ik,iq)

      ! Sanity checks
      if(size(eveck,1) /= size(eveckq,1)) then
        write(*,*) "Error(m_ematqk::emat_gen):&
          & Differing basis size (size(eveck,1) /= size(eveckq,1))"
        call terminate
      end if
      nbasis = size(eveck,1)
      nstk = size(eveck,2)
      nstkq = size(eveckq,2)
      if(size(emat,1) /= nstk .or. size(emat,2) /= nstkq) then
        write(*,*) "Error(m_ematqk::emat_gen):&
          & emat array shape is not compatible with eigenvector arrays&
          & shape(emat) =", shape(emat), " shape(eveck) =", shape(eveck),&
          & " shape(eveckq) =", shape(eveckq)
        call terminate
      end if

      ngq = g_qq%ngknr(ispin, iq)

      if(size(emat,3) /= ngq) then
        write(*,*) "Error(m_ematqk::emat_gen):&
          & emat array shape is not compatible with ngq&
          & size(emat,3) =", size(emat,3), " ngq =", ngq
        call terminate
      end if

      emat = zzero

      ngk = g_k%ngk(ispin, ik)
      ngkq = g_kqmtp%ngk(ispin, ikq)

      !-----------------------------------------------------------------------------!
      ! Generating matching coefficients                                            !
      !-----------------------------------------------------------------------------!
      ! Find matching coefficients for k-point ik
      ! Note, match uses:
      !   input%groundstate%lmaxapw
      !   mod_APW_LO::lmmaxapw
      !   mod_Gkvector::ngkmax
      lmmaxapw_save = lmmaxapw
      lmmaxapw = (lmaxapw+1)**2
      ngkmax_save = ngkmax
      ngkmax = g_k%ngkmax

      allocate( apwalmk( g_k%ngkmax, apwordmax, (lmaxapw+1)**2, natmtot))
      call match(ngk, g_k%gkc(:, ispin, ik),&
       & g_k%tpgkc(1:2, :, ispin, ik), g_k%sfacgk(:,:,ispin, ik), apwalmk)

      ! Find the matching coefficients for k-point k+q
      ngkmax = g_kqmtp%ngkmax
      allocate( apwalmkq( g_kqmtp%ngkmax, apwordmax, (lmaxapw+1)**2, natmtot))
      call match(ngkq, g_kqmtp%gkc(:, ispin, ikq),&
       & g_kqmtp%tpgkc(1:2, :, ispin, ikq), g_k%sfacgk(:,:,ispin, ikq), apwalmkq)

      ! Resorte external module variables
      ngkmax = ngkmax_save
      lmmaxapw = lmmaxapw_save
      !-----------------------------------------------------------------------------!

      !-----------------------------------------------------------------------------!
      !  Compute PW matrix elements for G+q < gqmax                                 !
      !  (or exactly G=0, if gqmax = 0 was set)                                     !
      !-----------------------------------------------------------------------------!
      !--------------------------------------!
      !     Muffin tin matrix elements       !
      !--------------------------------------!
      ! Dependes on G+q

      ! Generate radial integrals for G+q, where q is in [0,1) cell
      ! And store them in the module variables
      call emat_radial_init(iq)
      
      ! Build block matrix (for each G+q)

      !call emat_mt_part(emat)
      call emat_mt_part_new(emat)

      iqprev = iq
       
      !--------------------------------------!
      !     Interstitial matrix elements     !
      !--------------------------------------!
      ! Dependes on k, q, k+q, G_shift
      call emat_inter_part(emat)
      !-----------------------------------------------------------------------------!

      deallocate(apwalmk, apwalmkq)

      contains 


        subroutine emat_mt_part(emat)

          complex(8), intent(inout) :: emat(:,:,:)

          integer(4) :: is, ia, ias
          integer(4) :: l, m, o, lmo, lm
          integer(4) :: igq, i
          integer(4) :: ilok1, ilok2, ilokq1, ilokq2

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
!$OMP& PRIVATE( ilok1, ilok2, ilokq1, ilokq2), &
!$OMP& PRIVATE( lmo, l, m, o, lm)
#endif
          allocate(blockmt(ngk+nlotot,ngkq+nlotot))
          allocate(auxmat(nlmomax, ngkq))
          allocate(auxmat2(ngk+nlotot, nstkq))
          allocate(match_combined1(nlmomax, ngk))
          allocate(match_combined2(nlmomax, ngkq))
          allocate(aamat(ngk, ngkq), almat(ngk, nidxlomax), lamat(nidxlomax, ngkq))

#ifdef USEOMP
!$OMP DO 
#endif
          ! Do for all G+q vectors for current q
          Gqloop: do igq = 1, ngq

            ! Build block matrix for current G+q
            blockmt(:,:) = zzero

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
                  lm = myidxlm(l, m)
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

                ! Index range of local orbitals of atom ias in eigenvectors
                ilok1 = ngk+idxlostart(ias)
                ilok2 = ngk+idxlostart(ias)+nidxlo(is)-1
                ilokq1 = ngkq+idxlostart(ias)
                ilokq2 = ngkq+idxlostart(ias)+nidxlo(is)-1

                ! APW-LO
                call ZGEMM('C', 'N', ngk, nidxlo(is), nlmo(is), zone, &
                    match_combined1(1:nlmo(is), :), nlmo(is), &
                    rigntal(1:nlmo(is), 1:nidxlo(is), ias, igq), nlmo(is), zzero, &
                    almat(1:ngk,1:nidxlo(is)), ngk)

                blockmt(1:ngk,ilokq1:ilokq2) = blockmt(1:ngk,ilokq1:ilokq2)+almat(1:ngk,1:nidxlo(is))

                ! LO-APW
                call ZGEMM('N','N', nidxlo(is), ngkq, nlmo(is), zone, &
                    rigntla(1:nidxlo(is), 1:nlmo(is), ias, igq), nidxlo(is), &
                    match_combined2(1:nlmo(is), :), nlmo(is), zzero, &
                    lamat(1:nidxlo(is),1:ngkq), nidxlo(is))

                blockmt(ilok1:ilok2, 1:ngkq) =&
                  & blockmt(ilok1:ilok2, 1:ngkq) + lamat(1:nidxlo(is),1:ngkq)

                ! LO-LO
                blockmt(ilok1:ilok2, ilokq1:ilokq2) =&
                  & blockmt(ilok1:ilok2, ilokq1:ilokq2)&
                  &+ rigntll(1:nidxlo(is), 1:nidxlo(is), ias, igq)

              end do
            end do

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

          integer(4) :: ig, igp, igk, igkq
          integer(4) :: igs, ivg(3), ivgs(3)
          integer(4) :: ist1, ist2
          complex(8), allocatable :: zfftk(:,:), zfftkq(:,:), zfftres(:)
          complex(8), allocatable :: zfftktheta(:,:)

          ! What is done here:
          ! Consider the integral over the interstitial region of <mk|e^{-i(G+q)r}|nk>.
          ! Writing it using the Lattice theta function as an integral over
          ! the whole volume: \int_V \theta(r) \Phi^*_{mk}(r) e^{-i(G+q)r} \Phi_{nk}(r)
          !
          ! \int_V \theta(r) \Phi^*_{mk}(r) e^{-i(G+q)r} \Phi_{nk}(r) 
          !
          ! Insert lattice FT expressions -->
          ! =\int_V (\Sum_{G1} \theta(G1) e^{iG1r}) (\Sum_{G2} \Phi_{mk}(G2) e^{i(G2+k)r})^*
          !        e^{-i(G+q)r}(\Sum_{G3} \Phi_{nk}(G3) e^{i(G3+k')r})
          ! =\int_V (\Sum_{G1} \theta(G1) e^{iG1r}) (\Sum_{G2} \Phi_{mk}(G2) e^{iG2r})^*
          !         (\Sum_{G3} \Phi_{nk}(G3) e^{iG3)r}) e^{-i(k+q-k')r} e^{-iGr}
          ! =\int_V (\Sum_{G1} \theta(G1) e^{iG1r}) (\Sum_{G2} \Phi_{mk}(G2) e^{iG2r})^*
          !         (\Sum_{G3} \Phi_{nk}(G3) e^{iG3)r}) e^{-i(G+Gshift)r}
          !
          ! By design it is k+q-k' = Gshift with some Gshift lattice vector.
          ! All sums in the round brackets are directly computable as
          ! inverse FFT's using the eigenvector coefficients (LAPW part) and \theta(G).
          ! The corresponding real space quantities are then multiplied and subsequently
          ! transformed back to reciprocal space.
          ! If Gshift is non zero the final back transform needs to be shifted.
          ! 
          ! So 3 inverse FFT's are computed, the real space results multiplied and
          ! subsequently back transformed with another FFT to get the integral at all G

          ! Get Gshift index of k + q = k' + Gshift
          igs = q_q%ikiq2ig_nr(ik,iq)
          ! Get Gshift vector
          ivgs = g%ivg(1:3,igs)

          allocate(zfftk(fftmap%ngrtot+1, nstk))
          allocate(zfftkq(fftmap%ngrtot+1, nstkq))
          allocate(zfftktheta(fftmap%ngrtot+1, nstk))
          zfftk = zzero
          zfftktheta = zzero
          zfftkq= zzero

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,ist2,igk,igkq,igp,igs,ig,ivg,zfftres)
#endif
          if( .not. precompute_ffts) then
#ifdef USEOMP
!$OMP DO
#endif
            do ist1 = 1, nstk
              ! Map evec to FFT grid
              do igk = 1, ngk 
                zfftk(fftmap%igfft(g_k%igkig(igk, ispin, ik)), ist1) = eveck(igk, ist1)
              end do
              ! Do the inverse FFT
              call zfftifc(3, fftmap%ngrid, 1, zfftk(:, ist1))
            end do
#ifdef USEOMP
!$OMP END DO NOWAIT
#endif
#ifdef USEOMP
!$OMP DO
#endif
            do ist2 = 1, nstkq
              ! Map evec to FFT grid
              do igkq = 1, ngkq 
                zfftkq(fftmap%igfft(g_kqmtp%igkig(igkq, ispin, ikq)), ist2) = eveckq(igkq, ist2)
              end do
              ! Do the inverse FFT
              call zfftifc(3, fftmap%ngrid, 1, zfftkq(:, ist2))
            end do

#ifdef USEOMP
!$OMP END DO
#endif
          else
#ifdef USEOMP
!$OMP MASTER
#endif
            call emat_fetch_fft(zfftk, ik_=ik)
            call emat_fetch_fft(zfftkq,ikq_=ikq)
#ifdef USEOMP
!$OMP END MASTER
#endif
          end if
#ifdef USEOMP
!$OMP DO
#endif
          do ist1 = 1, nstk
            ! (evec(r,i)*\Theta(r))^*
            zfftktheta(:, ist1) = conjg(zfftk(:, ist1))*zfftcf(:) 
          end do
#ifdef USEOMP
!$OMP END DO
#endif
          ! Generate the 3d fourier transform of the real space integrand 
          ! zfftres = eveck^*_ist1(r)*\Theta(r)*eveckq_ist2(r)
          ! for all ist1/2 combinations and add this to emat.
          allocate(zfftres(fftmap%ngrtot+1))
#ifdef USEOMP
!$OMP DO
#endif
          do ist2 = 1, nstkq
            do ist1 = 1, nstk
              ! Generate real space integrand
              ! eveck^*_ist1(r)*\Theta(r)*eveckq_ist2(r)
              ! and do a 3d fft to G space
              do igp = 1, fftmap%ngrtot
                zfftres(igp) = zfftktheta(igp, ist1)*zfftkq(igp, ist2) 
              end do
              zfftres(fftmap%ngrtot+1) = zzero
              call zfftifc(3, fftmap%ngrid, -1, zfftres)

              ! Add result to emat
              do igp = 1, ngq
                if(any(ivgs /= 0)) then 
                  ! Get G vector index corresponding to G+q
                  ig=g_qq%igknrig(igp, ispin, iq)
                  ! Apply shift 
                  ivg=g%ivg(:,ig) + ivgs
                  ! Get G vector index of shifted G
                  igs=g%ivgig(ivg(1),ivg(2),ivg(3))
                  ! Map shifted FFT results
                  emat(ist1, ist2, igp) = emat(ist1, ist2, igp)&
                    &+ zfftres(fftmap%igfft(igs))
                else
                  emat(ist1, ist2, igp) = emat(ist1, ist2, igp)&
                    &+ zfftres(fftmap%igfft(g_qq%igknrig(igp, ispin, iq)))
                end if
              end do

            end do
          end do
#ifdef USEOMP
!$OMP END DO
#endif
          deallocate(zfftres)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
          deallocate(zfftk, zfftktheta, zfftkq)
        end subroutine emat_inter_part

        subroutine emat_mt_part_new(emat)

          complex(8), intent(inout) :: emat(:,:,:)

          integer(4) :: is, ia, ias
          integer(4) :: l, m, o, lmo, lm
          integer(4) :: igq, i, nlmolo

          ! Block matrix
          ! [_AA_|_AL_]
          ! [ LA | LL ]
          complex(8), allocatable, dimension(:,:) :: blockmt
          ! Helper matrices
          complex(8), allocatable, dimension(:,:) :: auxmat
          complex(8), allocatable, dimension(:,:) :: matchk, matchkq
          complex(8), allocatable, dimension(:,:,:) :: acveck, acveckq

          ! Shared arrays
          allocate(acveck(nlmomax+nidxlomax, nstk, natmtot))
          allocate(acveckq(nlmomax+nidxlomax, nstkq, natmtot))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP& PRIVATE( matchk, matchkq), &
!$OMP& PRIVATE( is, ia, ias), &
!$OMP& PRIVATE( lmo, l, m, o, lm)
#endif
          allocate(matchk(nlmomax, ngk))
          allocate(matchkq(nlmomax, ngkq))
#ifdef USEOMP
!$OMP DO COLLAPSE(2)
#endif
          ! Make product of Matching coefficients A and 
          ! eigenvector coefficients C, as computation helper
          ! AC_{lmo,ist;ias} = \Sum_G+k A_{lmo,G+k;ias}*C_{G+k,ist}
          ! Also append the coefficients of the corresponding local orbitals.
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)

              do lmo = 1, nlmo(is)
                l = lmo2l(lmo, is)
                m = lmo2m(lmo, is)
                o = lmo2o(lmo, is)
                lm = myidxlm(l, m)
                matchk(lmo, :) = apwalmk(1:ngk, o, lm, ias)
                matchkq(lmo, :) = apwalmkq(1:ngkq, o, lm, ias)
              end do 

              ! AC_{lmo,ist;ias} = \Sum_G+k A_{lmo,G+k;ias}*C_{G+k,ist}
              call ZGEMM('N', 'N', nlmo(is), nstk, ngk, zone, &
                  matchk(1:nlmo(is), 1:ngk), nlmo(is), &
                  eveck(1:ngk, 1:nstk), ngk, zzero, &
                  acveck(1:nlmo(is), 1:nstk, ias), nlmo(is))
              ! Append eigenvector coefficients corresponding to the local 
              ! orbitals of the atom ias
              acveck((nlmo(is)+1):(nlmo(is)+nidxlo(is)), 1:nstk, ias) = &
                & eveck((ngk+idxlostart(ias)):(ngk+idxlostart(ias)+nidxlo(is)-1), 1:nstk)

              ! Same for the k+q eigenvectors
              call ZGEMM('N', 'N', nlmo(is), nstkq, ngkq, zone, &
                  matchkq(1:nlmo(is), 1:ngkq), nlmo(is), &
                  eveckq(1:ngkq, 1:nstkq), ngkq, zzero, &
                  acveckq(1:nlmo(is), 1:nstkq, ias), nlmo(is))
              acveckq((nlmo(is)+1):(nlmo(is)+nidxlo(is)), 1:nstkq, ias) = &
                & eveckq((ngkq+idxlostart(ias)):(ngkq+idxlostart(ias)+nidxlo(is)-1), 1:nstkq)

            end do
          end do
#ifdef USEOMP
!$OMP END DO
#endif
          deallocate(matchk, matchkq)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED), &
!$OMP& PRIVATE( blockmt), &
!$OMP& PRIVATE( auxmat), &
!$OMP& PRIVATE( is, ia, ias, igq, nlmolo)
#endif
          allocate(blockmt(nlmomax+nidxlomax,nlmomax+nidxlomax))
          allocate(auxmat(nlmomax+nidxlomax, nstkq))
#ifdef USEOMP
!$OMP DO 
#endif
          ! Do for all G+q vectors for current q
          Gqloop: do igq = 1, ngq

            do is = 1, nspecies
              do ia = 1, natoms(is)
                ias = idxas(ia, is)

                ! Size of block matrix
                nlmolo = nlmo(is)+nidxlo(is)

                ! Setup block matrix
                ! APW-APW
                blockmt(1:nlmo(is), 1:nlmo(is)) = &
                  & rigntaa(1:nlmo(is), 1:nlmo(is), ias, igq)
                ! APW-LO
                blockmt(1:nlmo(is), (nlmo(is)+1):nlmolo) = &
                  & rigntal(1:nlmo(is), 1:nidxlo(is), ias, igq)
                ! LO-APW
                blockmt((nlmo(is)+1):nlmolo, 1:nlmo(is)) =&
                  & rigntla(1:nidxlo(is), 1:nlmo(is), ias, igq)
                ! LO-LO
                blockmt((nlmo(is)+1):nlmolo, (nlmo(is)+1):nlmolo) =&
                  & rigntll(1:nidxlo(is), 1:nidxlo(is), ias, igq)

                ! Compute muffin tin contribution from atom ias and add it to emat
                call ZGEMM('N', 'N', nlmolo, nstkq, nlmolo,&
                    zone, &
                    blockmt(1:nlmolo, 1:nlmolo), nlmolo, &
                    acveckq(1:nlmolo, 1:nstkq, ias), nlmolo, zzero, &
                    auxmat(1:nlmolo, 1:nstkq), nlmolo)

                call ZGEMM('C', 'N', nstk, nstkq, nlmolo, cmplx(fourpi, 0.d0, 8), &
                    acveck(1:nlmolo, 1:nstk, ias), nlmolo, &
                    auxmat(1:nlmolo,1:nstkq), nlmolo, zone, &
                    emat(1:nstk,1:nstkq, igq), nstk)

              end do
            end do

          end do Gqloop
#ifdef USEOMP
!$OMP END DO
#endif
          deallocate(auxmat)
          deallocate(blockmt)
#ifdef USEOMP
!$OMP END PARALLEL
#endif
          deallocate(acveck)
          deallocate(acveckq)
      
        end subroutine emat_mt_part_new

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
