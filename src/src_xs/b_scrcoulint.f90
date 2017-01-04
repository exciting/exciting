! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_scrcoulint
! !INTERFACE:
subroutine b_scrcoulint(iqmt, fra, fti)
! !USES:
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi
  use mod_constants, only: zzero, zone, fourpi
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: iqmap, ngridq, vql
  use mod_kpoint, only: nkptnr, ivknr
  use mod_lattice, only: omega
  use mod_symmetry, only: maxsymcrys
  use modxs, only: xsgnt, unitout,&
                 & ngqmax,&
                 & nqptr, qpari, qparf, ivqr,&
                 & ngq, ppari, pparf, iqmapr,&
                 & vqlr,&
                 & bsedl, bsedu, bsedd,&
                 & bsed, bcbs
  use m_xsgauntgen
  use m_findgntn0
  use m_writevars
  use m_genfilname
  use m_getunit
  use m_b_ematqk
  use m_putgetbsemat
  use modbse

use m_writecmplxparts
! !DESCRIPTION:
!   Calculates the direct term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!      level for the treatment of core excitations (using local orbitals).
!      October 2010 (Weine Olovsson)
!   Forked from scrcoulint.F90 and adapted for non TDA BSE. (Aurich)
!EOP
!BOC      

  implicit none

  ! I/O
  integer, intent(in) :: iqmt ! Index of momentum transfer Q
  logical, intent(in) :: fra  ! Construct RA coupling block
  logical, intent(in) :: fti  ! Use time inverted anti-resonant basis

  ! Local variables
  character(*), parameter :: thisnam = 'b_scrcoulint'

  ! ik,jk block of W matrix (final product)
  complex(8), allocatable :: sccli(:,:)

  ! Auxilliary arrays for the construction of sccli
  complex(8), allocatable :: sccli_t1(:, :), sccli_t2(:, :, :, :)
  complex(8), allocatable :: zm(:,:)
  ! Plane wave arrays
  complex(8), allocatable :: muu(:, :, :), cmuu(:, :)
  complex(8), allocatable :: moo(:, :, :), cmoo(:, :)
  complex(8), allocatable :: mou(:, :, :), cmou(:, :)
  complex(8), allocatable :: muo(:, :, :), cmuo(:, :)
  ! The fourier coefficients of the screend coulomb potential
  complex(8), allocatable :: wfc(:, :)
  ! Auxilliary arrays/vars for the construction of wfc
  complex(8), allocatable :: scieffg(:, :, :)
  complex(8), allocatable :: phf(:, :)
  ! Symmerty maps creation 
  integer(4) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer(4), allocatable :: igqmap(:)
  integer(4) :: jsym, jsymi
  integer(4) :: nsc, ivgsym(3)
  logical :: tphf
  ! Mappings of jk ik combinations to q points
  integer(4) :: ikkp, iknr, jknr, ik, jk
  integer(4) :: iqrnr, iqr, iq
  real(8) :: vqr(3), vq(3)
  integer(4) :: numgq
  integer(4) :: iv(3)
  logical :: tq0
  ! Number of occupied/unoccupied states at ik and jk
  integer(4) :: ino, inu, jno, jnu
  ! Number of transitions at ik and jk
  integer(4) :: inou, jnou
  ! Number of (o/u)_i (o/u)_j combinations
  integer(4) :: noo, nuu, nou, nuo
  ! State loop indices
  integer(4) :: io, jo, iu, ju
  ! Combinded loop indices 
  integer(4) :: jaoff, iaoff, ia, ja
  ! Aux.
  integer(4) :: j1, j2
  complex(8) :: pref 
  ! Timing vars
  real(8) :: tscc1, tscc0

  ! External functions
  integer, external :: idxkkp
  logical, external :: tqgamma

  logical :: fwp

  fwp = input%xs%bse%writeparts

  !---------------!
  !   main part   !
  !---------------!

write(*,*) "Hello, this is b_scrcoulint at rank:", rank

  ! General setup
  call init0
  ! k-point setup
  call init1
  ! q-point and qmt-point setup
  call init2

  ! Check number of empty states
  if(input%xs%screening%nempty .lt. input%groundstate%nempty) then
    write(*,*)
    write(*, '("Error(",a,"): Too few empty states in&
      & screening eigenvector file - the screening should&
      & include many empty states (BSE/screening)", 2i8)')&
      & trim(thisnam), input%groundstate%nempty, input%xs%screening%nempty
    write(*,*)
    call terminate
  end if

  ! Save variables for the gamma q-point, i.e. the unshifted (apart from gs:vkloff)
  ! k set.
  call xssave0

  ! Generate gaunt coefficients used in the construction of 
  ! the plane wave matrix elements in ematqk.
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
    & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  ! Find indices for non-zero gaunt coefficients
  call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
    & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

  if(mpiglobal%rank == 0) then
    write(unitout, '(a,3i8)') 'Info(' // thisnam // '):&
      & Gaunt coefficients generated within lmax values:', input%groundstate%lmaxapw,&
      & input%xs%lmaxemat, input%groundstate%lmaxapw
    write(unitout, '(a, i6)') 'Info(' // thisnam // '): Number of reduced q-points: ', nqptr
    call flushifc(unitout)
  end if

  ! Read Fermi energy from file EFERMI
  ! Use EFERMI_QMT000.OUT
  call genfilname(iqmt=0, setfilext=.true.)
  call readfermi

  ! Set EVALSV_QMTXXX.OUT as basis for the occupation limits search
  !   Note: EVALSV_QMT000.OUT is always produced,
  !         EVALSV_QMT001.OUT has the same content, when
  !         the first entry in the q-point list is set 0 0 0
  ! To be exact the following genfilname set filext to _QMTXXX.OUT
  call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

  ! Set ist* variables and ksgap in modxs using findocclims
  ! This also reads in 
  ! (QMTXXX)
  ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
  ! (QMT000)
  ! modxs:evalsv0, modxs:occsv0
  call setranges_modxs(iqmt)

  ! Select relevant transitions for the construction
  ! of the BSE hamiltonian
  ! Also sets nkkp_bse, nk_bse 
  call select_transitions(iqmt, serial=.false.)

  ! Write support information to file
  if(mpiglobal%rank == 0) then
    if(fra) then
      if(fti) then 
        call genfilname(basename=trim(infofbasename)//'_'//trim(scclictifbasename),&
          & iqmt=iqmt, filnam=infofname)
        call b_putbseinfo(infofname, iqmt)
      else
        call genfilname(basename=trim(infofbasename)//'_'//trim(scclicfbasename),&
          & iqmt=iqmt, filnam=infofname)
        call b_putbseinfo(infofname, iqmt)
      end if
    else
      call genfilname(basename=trim(infofbasename)//'_'//trim(scclifbasename),&
        & iqmt=iqmt, filnam=infofname)
      call b_putbseinfo(infofname, iqmt)
    end if
  end if

  ! Set output file name
  if(fra) then
    if(fti) then
      call genfilname(basename=scclictifbasename, iqmt=iqmt, filnam=scclifname)
    else
      call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=scclifname)
    end if
  else
    call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=scclifname)
  end if

  ! Change file extension and write out k an q points
  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(mpiglobal%rank == 0) then
    call writekpts
    call writeqpts
  end if
  ! Revert to previous file extension
  call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

  ! Allocate local arrays for screened coulomb interaction and
  ! W(G,G',qr)
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  scieffg(:, :, :) = zzero
  ! Phases for transformations form reduced q points to non reduced ones.
  allocate(phf(ngqmax, ngqmax))

  !------------------------------------!
  ! GENERATE FOURIER COEFFICIENTS OF W !     
  ! (and radial integrals for emat)    !
  !------------------------------------!
  ! Parallelize over reduced q-point set
  call genparidxran('q', nqptr)

  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_scrcoulint):&
      & Calculating W(G1,G2,qr) fourier coefficients")')
    call timesec(tscc0)
  end if

  do iqr = qpari, qparf ! Reduced q

    ! Locate reduced q-point in non-reduced set
    iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))

    ! Get number of G+q vectors for current q
    numgq = ngq(iqrnr)

    ! Calculate effective screened coulomb interaction
    ! by inverting the symmetrized RPA dielectric matrix for a given q and
    ! 0 frequency and then multiplying
    ! it with v^{1/2} from both sides.
    call genscclieff(iqr, ngqmax, numgq, scieffg(:,:,iqr))

    ! Generate radial integrals for matrix elements of plane wave
    ! and save them to disk.
    call putematrad(iqr, iqrnr)

  end do

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(set=nqptr, rlen=ngqmax*ngqmax,&
    & zbuf=scieffg, inplace=.true., comm=mpiglobal)

  if(mpiglobal%rank == 0) then
    call timesec(tscc1)
    write(unitout, '("  Timing (in seconds):", f12.3)') tscc1 - tscc0
  end if

  !-------------------------------!
  ! CONSTRUCT W MATRIX ELEMENTS   !
  !-------------------------------!

  ! Normalization factor 1/V and per k point
  pref=1.0d0/(omega*dble(nk_bse))

  ! Allocate arrays used in ematqk (do not change in the following loop)
  call ematqalloc

  ! Work arrays (allocate for maximal size over all participating k points)
  allocate(sccli_t2(nu_bse_max, no_bse_max, nu_bse_max, no_bse_max))
  allocate(sccli_t1(no_bse_max**2, nu_bse_max**2))
  allocate(sccli(nou_bse_max, nou_bse_max))

  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_scrcoulint): W matrix elements")')
    call timesec(tscc0)
  end if

  ! Distributed loop over combinations of non-reduced k-point
  ! that contribute to the desired energy range (or bands).
  call genparidxran('p', nkkp_bse)

  kkploop: do ikkp = ppari, pparf

    ! Get individual k-point indices from combined kk' index.
    !   ik runs from 1 to nk_bse, jk from ik to nk_bse
    call kkpmap(ikkp, nk_bse, ik, jk)
    !! Note: If the screened coulomb interaction matrix 
    !! has the elements W_{i,j} and the indices enumerate the
    !! states according to
    !! i = {u1o1k1, u2o1k1, ..., uMo1k1,
    !!      u1o2k1, ..., uMoNk1, u1o1k2, ..., uMoNkO} -> {1,...,M*N*O}
    !! then because of W_{j,i} = W^*_{i,j} only kj = ki,..,kN is 
    !! needed (the diagonal blocks where k'=k will be fully computed 
    !! but only the upper triangle will be needed in the end).
    !! (The RA part in the standard basis is symmetric instead of 
    !!  hermitian, but one still just needs the upper triangle)

    ! Get total k point indices
    iknr = kmap_bse_rg(ik)
    jknr = kmap_bse_rg(jk) 

    !! Get corresponding q-point for ki,kj combination.
    ! Note: Why do we discard the Potential lattice component of jk-ik = iq + G ?
    if(fti .and. fra) then 
      ! RA^{ti}
      ! k-point difference -(k_j+k_i) on integer grid.
      iv(:) = -(ivknr(:,jknr) + ivknr(:,iknr))
    else
      ! RR and RA
      ! k-point difference k_j-k_i on integer grid.
      iv(:) = ivknr(:,jknr) - ivknr(:,iknr)
    end if
    ! Map to reciprocal unit cell [01) 
    iv(:) = modulo(iv(:), ngridq(:))
    ! Find corresponding q-point index (reduced)
    iqr = iqmapr(iv(1), iv(2), iv(3))
    ! Get corresponding vector in lattice coordinated
    vqr(:) = vqlr(:, iqr)
    ! q-point (non-reduced)
    iq = iqmap(iv(1), iv(2), iv(3))
    ! Get lattice coordinates
    vq(:) = vql(:, iq)

    ! Check if iq is Gamma point (mod_qpoint::vqc(:,iq) has length < 1d-12)
    tq0 = tqgamma(iq)

    ! Local field effects size (Number of G+q vectors)
    numgq = ngq(iq)

    allocate(igqmap(numgq))
    allocate(wfc(numgq, numgq))

    !! Find results  
    !! for a non reduced q-point with the help of the result 
    !! for the corresponding reduced q-point using symmetry operations.
    !! (Radial emat integrals and screened coulomb potential Fourier coefficients)
    !!<--
    !! RR & RA
    ! Find symmetry operations that map the reduced q-point to the non reduced one
    call findsymeqiv(input%xs%bse%fbzq, vq, vqr, nsc, sc, ivgsc)
    ! Find a crystal symmetry operation that rotates the G+q-vectors onto G'+q_r-vectors
    ! and generate a Map G' --> G
    call findgqmap(iq, iqr, nsc, sc, ivgsc, numgq, jsym, jsymi, ivgsym, igqmap)

    ! Get radial integrals for q_r (previously calculated for reduced q set)
    call getematrad(iqr, iq)
    ! Rotate radial integrals calculated for the reduced q to get those for non-reduced q
    call rotematrad(numgq, igqmap)

    ! Generate phase factor for dielectric matrix due to non-primitive
    ! translations
    call genphasedm(iq, jsym, ngqmax, numgq, phf, tphf)
    ! W(G,G',q) <-- W(\tilde{G},\tilde{G}',qr)
    wfc(:,:) = phf(:numgq, :numgq) * scieffg(igqmap, igqmap, iqr)
    !!-->

    ! Get ik & jk dependent band ranges for 
    ! plane wave matrix calculation
    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1
    jnu = koulims(2,jknr) - koulims(1,jknr) + 1
    jno = koulims(4,jknr) - koulims(3,jknr) + 1
    noo = ino*jno
    nuu = inu*jnu
    nou = ino*jnu
    nuo = inu*jno
    ! Number of transitions at ik & jk
    inou = kousize(iknr)
    jnou = kousize(jknr)

    ! W^RR
    if(.not. fra) then 
      !-------------------------------!
      ! Resonant-Resonant Part        !
      !-------------------------------!

      allocate(cmoo(noo, numgq), cmuu(nuu, numgq))

      ! Calculate M_{io jo ik}(G, q) and M_{iu ju ik+qmt}(G, q)
      ! for current current q=jk-ik and ik
      ! NOTE: only qmt=0 supported currently
      call getpwesrr(iq, iknr, moo, muu)

      ! Combine indices for matrix elements of plane wave.
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,jo,j1)
      do jo = 1, jno   ! jo
        do io = 1, ino ! io
          j1 = io + (jo-1)*ino ! iojo
          ! cmoo_j = M_o1o1, M_o2o1, ..., M_oNo1, M_o1o2, ..., M_oNoM
          cmoo(j1, :) = moo(io, jo, :)
        end do
      end do
      !$OMP END PARALLEL DO
      deallocate(moo)

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(iu,ju,j2)
      do ju = 1, jnu   ! ju
        do iu = 1, inu ! iu
          j2 = iu + (ju-1)*inu ! iuju
          ! cmuu_j = M_u1u1, M_u2u1, ..., M_uNu1, M_u1u2, ..., M_uNuM
          cmuu(j2, :) = muu(iu, ju, :)
        end do
      end do
      !$OMP END PARALLEL DO
      deallocate(muu)

      ! M_{iojo} -> M^*_{iojo}
      cmoo = conjg(cmoo)

      ! Allocate helper array
      allocate(zm(noo,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclit_{j1, j2}(q)
      ! zm = cmoo * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmoo_{j1,G} wfc_{G,G'}
      !   i.e. zm_{io_j1,jo_j1}(G',q) = \Sum_{G} M^*_{io_j1,jo_j1}(G,q) W(G,G', q)
      call zgemm('n', 'n', noo, numgq, numgq, zone, cmoo, noo,&
        & wfc, numgq, zzero, zm, noo)
      ! scclit = pref * zm * cmuu^T
      !   i.e. scclit_{j1, j2} = \Sum_{G'} zm_{j1,G'} (cmuu^T)_{G',j2}
      !   i.e. scclit_{io_j1 jo_j1, iu_j2 ju_j2} = 
      !          \Sum{G,G'} M^*_{io_j1,jo_j1}(G,q) W(G,G',q) M_{iu_j2 ju_j2}(G',q)
      call zgemm('n', 't', noo, nuu, numgq, pref, zm, noo,&
        & cmuu, nuu, zzero, sccli_t1(1:noo,1:nuu), noo)

      deallocate(zm)        
      deallocate(cmoo, cmuu)

      ! Map back to individual band indices
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(4),&
      !$OMP& DEFAULT(SHARED), PRIVATE(iu,ju,io,jo,j1,j2)
      do ju = 1, jnu    ! ju
        do iu = 1, inu  ! iu
          do jo = 1, jno   ! jo
            do io = 1, ino ! io
              j2 = iu + (ju-1)*inu
              j1 = io + (jo-1)*ino
              ! scclit_{io_j1 jo_j1, iu_j2 ju_j2} -> sccli_{iu_j2 io_j1, ju_j2 jo_j1}
              sccli_t2(iu, io, ju, jo) = sccli_t1(j1, j2)
            end do
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      ! W^RR matrix element arrays for one jk-ik=q

      ! Save only the selected transitions
      jaoff = sum(kousize(1:jknr-1))
      iaoff = sum(kousize(1:iknr-1))

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(iu,ju,io,jo,ia,ja)
      do ja = 1, jnou
        do ia = 1, inou
          ju = smap_rel(1,ja+jaoff)
          jo = smap_rel(2,ja+jaoff)
          iu = smap_rel(1,ia+iaoff)
          io = smap_rel(2,ia+iaoff)
          ! scclit_{io_j1 jo_j1, iu_j2 ju_j2} -> sccliab_{iu_a io_a, ju_b jo_b}
          sccli(ia, ja) = sccli_t2(iu, io, ju, jo)
        end do
      end do
      !$OMP END PARALLEL DO

      if(fwp) then 
        if(ikkp == 1) then 
          call writecmplxparts('ikkp1_Wrr',dble(sccli(1:inou,1:jnou)), aimag(sccli(1:inou,1:jnou)))
        end if
      end if

      ! Parallel write
      call b_putbsemat(scclifname, 77, ikkp, iqmt, sccli)

    ! W^{RA} or W^{RA,ti}
    else
      !-------------------------------!
      ! Resonant-Anti-Resonant Part   !
      !-------------------------------!

      allocate(cmou(nou, numgq), cmuo(nuo, numgq))

      ! Calculate M_{io ju ik}(G, q-qmt) and M_{iu jo ik+qmt}(G', q-qmt)
      ! for current current q=jk-ik (or q = -jk-ik for ti) and ik 
      ! NOTE: only qmt=0 supported currently
      call getpwesra(iq, iknr, mou, muo)

      ! Combine indices for matrix elements of plane wave.
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,ju,j1)
      do ju = 1, jnu   ! ju
        do io = 1, ino ! io
          j1 = io+(ju-1)*ino
          ! cmou_j = M_o1u1, M_o2u1, ..., M_oNu1, M_o1u2, ..., M_oNuM
          cmou(j1, :) = mou(io, ju, :)
        end do
      end do
      !$OMP END PARALLEL DO
      deallocate(mou)

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,ju,j1)
      do jo = 1, jno   ! jo
        do iu = 1, inu ! iu
          j2 = iu+(jo-1)*inu
          ! cmuo_j = M_u1o1, M_u2o1, ..., M_uMo1, M_u1o2, ..., M_uMoN
          cmuo(j2, :) = muo(iu, jo, :)
        end do
      end do
      !$OMP END PARALLEL DO
      deallocate(muo)

      ! M_{ioju} -> M^*_{ioju}
      cmou = conjg(cmou)

      ! Allocate helper array
      allocate(zm(nou,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclict_{j1, j2}(q)
      ! zm = cmou * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmou_{j1,G} wfc_{G,G'}
      !   i.e. zm_{io_j1,ju_j1}(G',q) =
      !          \Sum_{G} M^*_{io_j1,ju_j1,ik}(G,q-qmt) W(G,G',q-qmt)
      ! NOTE: only qmt=0 supported currently
      call zgemm('n', 'n', nou, numgq, numgq, zone, cmou, nou,&
        & wfc, numgq, zzero, zm, nou)
      ! scclit = pref * zm * cmuo^T
      !   i.e. scclict(j1, j2) = \Sum_{G'} zm_{j1,G'} (cmuo^T)_{G',j2}
      !   i.e. scclict_{io_j1 ju_j1, iu_j2 jo_j2} =
      !   \Sum{G,G'} 
      !   M^*_{io_j1,ju_j1,ik}(G,q-qmt) W(G, G',q-qmt) M_{iu_j2,jo_j2,ik+qmt}(G',q-qmt)
      ! NOTE: only qmt=0 supported currently
      call zgemm('n', 't', nou, nuo, numgq, pref, zm,&
        & nou, cmuo, nuo, zzero, sccli_t1(1:nou,1:nuo), nou)
      deallocate(zm)        
      deallocate(cmou, cmuo)

      ! Map back to individual band indices
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(4),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,jo,iu,ju,j1,j2)
      do ju = 1, jnu   ! ju
        do io = 1, ino ! io
          do jo = 1, jno   ! jo
            do iu = 1, inu ! iu
              j1 = io+(ju-1)*ino !ioju
              j2 = iu+(jo-1)*inu ! iujo
              ! scclict_{io_j1 ju_j1, iu_j2 jo_j2}(ik,jk)(qmt)
              ! -> scclic(iu_j2, io_j1, ju_j1, jo_j2)(ik,jk)(qmt)
              sccli_t2(iu, io, ju, jo) = sccli_t1(j1, j2)
            end do
          end do
        end do
      end do 
      !$OMP END PARALLEL DO

      ! W matrix elements

      ! Save only the selected transitions
      jaoff = sum(kousize(1:jknr-1))
      iaoff = sum(kousize(1:iknr-1))
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,jo,iu,ju,ia,ja)
      do ja = 1, jnou
        do ia = 1, inou
          ju = smap_rel(1,ja+jaoff)
          jo = smap_rel(2,ja+jaoff)
          iu = smap_rel(1,ia+iaoff)
          io = smap_rel(2,ia+iaoff)
          ! scclit_{io_j1 jo_j1, iu_j2 ju_j2} -> sccliab_{iu_a io_a, ju_b jo_b}
          sccli(ia, ja) = sccli_t2(iu, io, ju, jo)
        end do
      end do
      !$OMP END PARALLEL DO

      if(fwp) then 
        if(ikkp == 1 .and. fti ) then 
          call writecmplxparts('ikkp1_Wra_ti',dble(sccli(1:inou,1:jnou)), aimag(sccli(1:inou,1:jnou)))
        else if(ikkp == 1) then 
          call writecmplxparts('ikkp1_Wra',dble(sccli(1:inou,1:jnou)), aimag(sccli(1:inou,1:jnou)))
        end if
      end if

      ! Parallel write
      call b_putbsemat(scclifname, 78, ikkp, iqmt, sccli)

    end if

    ! Deallocate G+q dependent work arrays
    deallocate(igqmap)
    deallocate(wfc)

    !if(rank == 0) then
    !  write(6, '(a,"Scrcoulint progess:", f10.3)', advance="no")&
    !    & achar( 13), 100.0d0*dble(ikkp-ppari+1)/dble(pparf-ppari+1)
    !  flush(6)
    !end if

  ! End loop over(k,kp)-pairs
  end do kkploop

  !if(rank == 0) then
  !  write(*,*)
  !end if

  if(mpiglobal%rank == 0) then
    call timesec(tscc1)
    write(unitout, '("  Timing (in seconds):", f12.3)') tscc1 - tscc0
  end if

  ! Deallocate helper array 
  deallocate(sccli)
  deallocate(sccli_t1)
  deallocate(sccli_t2)

  call barrier

  call findgntn0_clear

  if(mpiglobal%rank .eq. 0) then
    write(unitout, '("Info(b_scrcoulint): Screened coulomb interaction&
      & finished for iqmt=",i8)') iqmt
  end if

  contains

    subroutine getpwesrr(iqnr, iknr, moo, muu)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: moo(:,:,:), muu(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the oo plane wave elements
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=jno
      ematbc%il2=koulims(3,jknr)
      ematbc%iu2=koulims(4,jknr)
      ! Allocate space for M_{o1o2,G} at fixed (k, q)
      if(allocated(moo)) deallocate(moo)
      allocate(moo(ino,jno,numgq))
      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, moo, ematbc)

      !! Calculate the uu plane wave elements
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=jnu
      ematbc%il2=koulims(1,jknr)
      ematbc%iu2=koulims(2,jknr)
      ! Allocate space for M_{u1u2,G} at fixed (k, q)
      if(allocated(muu)) deallocate(muu)
      allocate(muu(inu,jnu,numgq))
      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, muu, ematbc)

    end subroutine getpwesrr

    subroutine getpwesra(iqnr, iknr, mou, muo)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: mou(:,:,:), muo(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the o-u plane wave elements
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=jnu
      ematbc%il2=koulims(1,jknr)
      ematbc%iu2=koulims(2,jknr)
      ! Allocate space for M_{ou,G} at fixed (k, q)
      if(allocated(mou)) deallocate(mou)
      allocate(mou(ino,jnu,numgq))
      ! Calculate M_{ou,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, mou, ematbc)

      !! Calculate the u-o plane wave elements
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=jno
      ematbc%il2=koulims(3,jknr)
      ematbc%iu2=koulims(4,jknr)
      ! Allocate space for M_{uo,G} at fixed (k, q)
      if(allocated(muo)) deallocate(muo)
      allocate(muo(inu,jno,numgq))
      ! Calculate M_{uo,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, muo, ematbc)

    end subroutine getpwesra

end subroutine b_scrcoulint
!EOC

