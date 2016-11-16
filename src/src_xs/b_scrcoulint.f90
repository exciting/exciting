! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_scrcoulint
! !INTERFACE:
subroutine b_scrcoulint
! !USES:
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi
  use mod_constants, only: zzero, zone, fourpi
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: nqpt, iqmap, ngridq, vql
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

  ! Local variables
  character(*), parameter :: thisnam = 'b_scrcoulint'

  complex(8) :: zt1, pref
  complex(8), allocatable :: sccliab(:,:), scclit(:, :), sccli(:, :, :, :), scclid(:, :)
  complex(8), allocatable :: sccliabc(:,:), scclitc(:, :), scclic(:, :, :, :)
  complex(8), allocatable :: scieffg(:, :, :), wfc(:, :), bsedt(:, :),zm(:,:)
  complex(8), allocatable :: phf(:, :)
  real(8) :: vqr(3), vq(3)
  integer(4) :: ik, jk, noo, nuu, nou, nuo, ino, inu, jno, jnu, inou, jnou
  integer(4) :: io, jo, iu, ju
  integer(4) :: jaoff, iaoff, ia, ja
  integer :: ikkp, iknr, jknr, iqr, iq, iqrnr, jsym, jsymi, numgq
  integer(4) :: iqmt
  integer :: nsc, iv(3), ivgsym(3), j1, j2
  integer :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer, allocatable :: igqmap(:)
  logical :: tq0, tphf
  logical :: fcoup

  ! Plane wave arrays
  complex(8), allocatable :: muu(:, :, :), cmuu(:, :)
  complex(8), allocatable :: moo(:, :, :), cmoo(:, :)
  complex(8), allocatable :: mou(:, :, :), cmou(:, :)
  complex(8), allocatable :: muo(:, :, :), cmuo(:, :)

  ! External functions
  integer, external :: idxkkp
  logical, external :: tqgamma

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

  if(rank .eq. 0) then
    write(unitout, '(a,3i8)') 'info(' // thisnam // '):&
      & Gaunt coefficients generated within lmax values:', input%groundstate%lmaxapw,&
      & input%xs%lmaxemat, input%groundstate%lmaxapw
    write(unitout, '(a, i6)') 'info(' // thisnam // '): number of q-points: ', nqpt
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
  iqmt = 0
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
  call select_transitions(iqmt)

  ! Include coupling terms
  fcoup = input%xs%bse%coupling

  ! Write support information to file
  call genfilname(basename=trim(infofbasename)//'_'//trim(scclifbasename),&
    & iqmt=iqmt, filnam=infofname)
  call b_putbseinfo(infofname, iqmt)
  if(fcoup) then
    call genfilname(basename=trim(infofbasename)//'_'//trim(scclicfbasename),&
      & iqmt=iqmt, filnam=infocfname)
    call b_putbseinfo(infocfname, iqmt)
  end if

  ! Set output file names
  call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=scclifname)
  if(fcoup) then
    call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=scclicfname)
  end if

  ! Change file extension and write out k an q points
  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if
  ! Revert to previous file extension
  call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

  ! Allocate local arrays for screened coulomb interaction
  ! Phases for transformations form reduced q points to non reduced ones.
  allocate(phf(ngqmax, ngqmax))
  ! W(G,G',q)
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  scieffg(:, :, :) = zzero

  !------------------------------------!
  ! GENERATE FOURIER COEFFICIENTS OF W !     
  ! (and radial integrals for emat)    !
  !------------------------------------!
  !!<--
  ! Parallelize over reduced q-point set
  call genparidxran('q', nqptr)

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
    call putematrad(iqr, iqrnr)

  end do

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(set=nqptr, rlen=ngqmax*ngqmax,&
    & zbuf=scieffg, inplace=.true., comm=mpiglobal)

  call barrier
  !!-->

  !-------------------------------!
  ! CONSTRUCT W MATRIX ELEMENTS   !
  !-------------------------------!

  pref=1.0d0/(omega*dble(nk_bse))

! Needs to be adapted for beyondtd
!  allocate(bsedt(3, 0:procs-1))
!  bsedt(1, :) = 1.d8
!  bsedt(2, :) = -1.d8
!  bsedt(3, :) = zzero

  ! Allocate arrays used in ematqk (do not change in the following loop)
  call ematqalloc

  ! Distributed loop over combinations of non-reduced k-point combinations
  ! that contribute to the desired energy range.
  call genparidxran('p', nkkp_bse)

  kkploop: do ikkp = ppari, pparf

!write(*,*) "ikkp:", ikkp, " of", pparf-ppari+1

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

    ! Get total k point indices
    iknr = kmap_bse_rg(ik)
    jknr = kmap_bse_rg(jk) 

!write(*,*) "ik:", ik, " jk", jk
!write(*,*) "iknr:", iknr, " jknr", jknr

    !! Get corresponding q-point for ki,kj combination.
    ! K-point difference k_j-k_i on integer grid.
    iv(:) = ivknr(:,jknr) - ivknr(:,iknr)
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
    ! Find symmetry operations that reduce the q-point to the irreducible
    ! part of the Brillouin zone
    call findsymeqiv(input%xs%bse%fbzq, vq, vqr, nsc, sc, ivgsc)
    ! Find the map that rotates the G-vectors
    call findgqmap(iq, iqr, nsc, sc, ivgsc, numgq, jsym, jsymi, ivgsym, igqmap)
    ! Get radial integrals (previously calculated for reduced q set)
    call getematrad(iqr, iq)
    ! Rotate radial integrals
    ! (i.e. apply symmetry transformation to the result for non reduced q point)
    call rotematrad(numgq, igqmap)
    ! Generate phase factor for dielectric matrix due to non-primitive
    ! translations
    call genphasedm(iq, jsym, ngqmax, numgq, phf, tphf)
    ! Rotate inverse of screening, coulomb potential and radial integrals
    ! (i.e. apply symmetry transformation to the result for non reduced q point)
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

!write(*,*) "inu,ino,inou", inu, ino, inou
!write(*,*) "jnu,jno,jnou", jnu, jno, jnou

    ! Work arrays
    if(allocated(scclit)) deallocate(scclit)
    allocate(scclit(noo, nuu))
    if(allocated(sccli)) deallocate(sccli)
    allocate(sccli(inu, ino, jnu, jno))

    if(fcoup) then
      ! Work array
      if(allocated(scclitc)) deallocate(scclitc)
      allocate(scclitc(nou, nuo))
      if(allocated(scclic)) deallocate(scclic)
      allocate(scclic(inu, ino, jnu, jno))
    end if

! Needs do be adapetd for beyondtd
!    if(allocated(scclid)) deallocate(scclid)
!    allocate(scclid(ino, inu))

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
   ! ! Does the same as the above (test for speed)
   ! cmoo = reshape(moo, [noo, numgq])
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
   ! ! Does the same as the above (test for speed)
   ! cmuu = reshape(muu, [nuu, numgq])
    deallocate(muu)

    ! M_{iojo} -> M^*_{iojo}
    cmoo = conjg(cmoo)

    ! Allocate helper array
    allocate(zm(noo,numgq))

    ! Calculate matrix elements of screened coulomb interaction scclit_{j1, j2}(q)
    ! zm = cmoo * wfc
    !   i.e. zm_{j1,G'} = \Sum_{G} cmoo_{j1,G} wfc_{G,G'}
    !   i.e. zm_{io_j1,jo_j1}(G',q) = \Sum_{G} M^*_{io_j1,jo_j1}(G,q) W(G,G', q)
    call zgemm('n', 'n', noo, numgq, numgq, zone, cmoo, noo, wfc, numgq, zzero, zm, noo)
    ! scclit = pref * zm * cmuu^T
    !   i.e. scclit_{j1, j2} = \Sum_{G'} zm_{j1,G'} (cmuu^T)_{G',j2}
    !   i.e. scclit_{io_j1 jo_j1, iu_j2 ju_j2} = 
    !          \Sum{G,G'} M^*_{io_j1,jo_j1}(G,q) W(G,G',q) M_{iu_j2 ju_j2}(G',q)
    call zgemm('n', 't', noo, nuu, numgq, pref, zm,&
      & noo, cmuu, nuu, zzero, scclit, noo)
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
            sccli(iu, io, ju, jo) = scclit(j1, j2)
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    deallocate(scclit)

    ! W matrix element arrays for one jk-ik=q
    if(allocated(sccliab)) deallocate(sccliab)
    allocate(sccliab(inou, jnou))

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
        sccliab(ia, ja) = sccli(iu, io, ju, jo)
      end do
    end do
    !$OMP END PARALLEL DO

! Needs do be adapetd for beyondtd
   ! ! Analyze BSE diagonal
   ! if(iknr .eq. jknr) then
   !   ! Selected occupied
   !   do io = 1, ino
   !     ! Selected unoccupied
   !     do iu = 1, inu
   !       zt1 = sccli(iu, io, iu, io)
   !       scclid(io, iu) = zt1
   !       bsedt(1, rank) = min(dble(bsedt(1, rank)), dble(zt1))
   !       bsedt(2, rank) = max(dble(bsedt(2, rank)), dble(zt1))
   !       bsedt(3, rank) = bsedt(3, rank) + zt1 / (ino*inu)
   !     end do
   !   end do
   ! end if

    ! Parallel write
    call b_putbsemat(scclifname, 77, ikkp, iqmt, sccliab)
    
    !-------------------------------!
    ! Resonant-Anti-Resonant Part   !
    !-------------------------------!
    if(fcoup) then

      allocate(cmou(nou, numgq), cmuo(nuo, numgq))

      ! Calculate M_{io ju ik}(G, q-qmt) and M_{iu jo ik+qmt}(G', q-qmt)
      ! for current current q=jk-ik and ik 
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
     ! cmou = reshape(mou,[nou, numgq])
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
     ! cmuo = reshape(muo,[nuo, numgq])
      deallocate(muo)

      ! M_{ioju} -> M^*_{ioju}
      cmou = conjg(cmou)

      ! Allocate helper array
      allocate(zm(nou,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclitc_{j1, j2}(q)
      ! zm = cmou * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmou_{j1,G} wfc_{G,G'}
      !   i.e. zm_{io_j1,ju_j1}(G',q) = \Sum_{G} M^*_{io_j1,ju_j1,ik}(G,q) W(G,G',q-qmt)
      ! NOTE: only qmt=0 supported currently
      call zgemm('n', 'n', nou, numgq, numgq, zone, cmou, nou, wfc, numgq, zzero, zm, nou)
      ! scclit = pref * zm * cmuo^T
      !   i.e. scclitc(j1, j2) = \Sum_{G'} zm_{j1,G'} (cmuo^T)_{G',j2}
      !   i.e. scclitc_{io_j1 ju_j1, iu_j2 jo_j2} =
      !   \Sum{G,G'} 
      !   M^*_{io_j1,ju_j1,ik}(G,q-qmt) W(G, G',q-qmt) M_{iu_j2,jo_j2,ik+qmt}(G',q-qmt)
      ! NOTE: only qmt=0 supported currently
      call zgemm('n', 't', nou, nuo, numgq, pref, zm,&
        & nou, cmuo, nuo, zzero, scclitc, nou)
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
              ! scclitc_{io_j1 ju_j1, iu_j2 jo_j2}(ik,jk)(qmt)
              ! -> scclic(iu_j2, io_j1, ju_j1, jo_j2)(ik,jk)(qmt)
              scclic(iu, io, ju, jo) = scclitc(j1, j2)
            end do
          end do
        end do
      end do 
      !$OMP END PARALLEL DO
      deallocate(scclitc)

      ! W matrix elements
      if(allocated(sccliabc)) deallocate(sccliabc)
      allocate(sccliabc(inou, jnou))

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
          sccliabc(ia, ja) = scclic(iu, io, ju, jo)
        end do
      end do
      !$OMP END PARALLEL DO

      ! Parallel write
      call b_putbsemat(scclicfname, 78, ikkp, iqmt, sccliabc)

    end if
    
    ! Deallocate G dependent work arrays
    deallocate(igqmap)
    deallocate(wfc)

  ! End loop over(k,kp)-pairs
  end do kkploop

  ! Deallocate helper array 
  deallocate(sccli)
  deallocate(sccliab)
  if(fcoup) then
    deallocate(scclic)
    deallocate(sccliabc)
  end if

  call barrier

! Needs to be adapted for beyondtd
!  ! Communicate array-parts wrt. q-points
!  call mpi_allgatherv_ifc(set=mpiglobal%procs,rlen=3,&
!    & zbuf=bsedt,inplace=.true.,comm=mpiglobal)
!
!  ! BSE kernel diagonal parameters
!  bsedl = minval(dble(bsedt(1, :)))
!  bsedu = maxval(dble(bsedt(2, :)))
!  bsedd = bsedu - bsedl
!  bsed = sum(bsedt(3, :)) / nkptnr
!  deallocate(bsedt, scclid)

!  ! Write BSE kernel diagonal parameters
!  if(rank .eq. 0) call putbsediag('BSEDIAG.OUT')

  call findgntn0_clear

  if(rank .eq. 0) then
    write(unitout, '("Info(b_scrcoulint): Screened coulomb interaction&
      & finished")')
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

