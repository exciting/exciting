! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: scrcoulint
! !INTERFACE:
subroutine b_scrcoulint
! !USES:
  use modinput, only: input
  use modmpi, only: procs, rank, mpi_allgatherv_ifc, barrier
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
  use m_ematqk
  use modbse
! !DESCRIPTION:
!   Calculates the direct term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!      level for the treatment of core excitations (using local orbitals).
!      October 2010 (Weine Olovsson)
!EOP
!BOC      

  implicit none

  ! Local variables
  character(*), parameter :: thisnam = 'scrcoulint'

  complex(8) :: zt1, pref
  complex(8), allocatable :: scclit(:, :), sccli(:, :, :, :), scclid(:, :)
  complex(8), allocatable :: scclitc(:, :), scclic(:, :, :, :)
  complex(8), allocatable :: scieffg(:, :, :), wfc(:, :), bsedt(:, :),zm(:,:)
  complex(8), allocatable :: phf(:, :)
  real(8) :: vqr(3), vq(3)
  integer :: ikkp, iknr, jknr, iqr, iq, iqrnr, jsym, jsymi, numgq, reclen
  integer :: nsc, iv(3), ivgsym(3), j1, j2, nkkp
  integer(4) :: io1, io2, iu1, iu2
  integer :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer, allocatable :: igqmap(:)
  logical :: tq0, tphf
  logical :: fcoup

  ! Plane wave arrays
  complex(8), allocatable :: muu(:, :, :), cmuu(:, :)
  complex(8), allocatable :: moo(:, :, :), cmoo(:, :)
  complex(8), allocatable :: mou(:, :, :), cmuo(:, :)
  complex(8), allocatable :: mou(:, :, :), cmuo(:, :)

  ! External functions
  integer, external :: idxkkp
  logical, external :: tqgamma

  !---------------!
  !   main part   !
  !---------------!

  ! General setup
  call init0
  ! K-point setup
  call init1
  ! Q-point setup
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

  ! Read Fermi energy from file
  call readfermi

  ! Save variables for the gamma q-point
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

  ! Set EVALSV_SCR.OUT as basis for the occupation limits search
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Set ist* variables in modxs using findocclims
  call setranges_modxs(0)

  ! Set band combinations (modbse:bcou & modbse:bcouabs)
  call setbcbs_bse

  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if


  ! Allocate local arrays for screened coulomb interaction
  ! Phased for transformations form reduced q points to non reduced ones.
  allocate(phf(ngqmax, ngqmax))
  ! W(G,G',q)
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  scieffg(:, :, :) = zzero

  ! Set file extension (for what ? It isn't used...)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

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
  call mpi_allgatherv_ifc(nqptr,ngqmax*ngqmax,zbuf=scieffg)
  call barrier
  !!-->

  !-------------------------------!
  ! CONSTRUCT W MATRIX ELEMENTS   !
  !-------------------------------!

  pref=1.0d0/(omega*dble(nkptnr))

  allocate(bsedt(3, 0:procs-1))
  bsedt(1, :) = 1.d8
  bsedt(2, :) = -1.d8
  bsedt(3, :) = zzero

  ! Allocate arrays used in ematqk (do not change in the following loop)
  call ematqalloc

  ! Include coupling terms
  fcoup = input%xs%bse%coupling

  ! Allocate arrays which do not change size in the following loop
  ! i.e. arrays that do not depend on numgq
  allocate(scclit(noo, nuu))
  scclit = zzero
  if(fcoup) then
    allocate(scclitc(nou, nou))
    scclitc = zzero
  end if

  ! W matrix element arrays for one k2-k1=q
  allocate(sccli(no, nu, no, nu), scclid(no, nu))
  sccli = zzero
  if(fcoup) then 
    allocate(scclic(no, nu, no, nu)
    scclic=zzero
  end if

  ! Distributed loop over combinations of non-reduced k-point combinations
  ! Combinations of k and k', where ik'>=ik
  nkkp = (nkptnr*(nkptnr+1)) / 2
  call genparidxran('p', nkkp)
  kkploop: do ikkp = ppari, pparf

    ! Get individual k-point indices from combined kk' index.
    !   iknr runs from 1 to nkptnr, jknr from iknr to nkptnr
    call kkpmap(ikkp, nkptnr, iknr, jknr)
    !! Note: If the screened coulomb interaction matrix 
    !! has the elements W_{i,j} and the indices enumerate the
    !! states according to
    !! i = {o1u1k1, o1u2k1, ..., o1uMk1,
    !!      o2uMk1, ..., oMuMk1, o1u1k2, ..., oMuMkN} -> {1,...,M**2N}
    !! then because of W_{j,i} = W^*_{i,j} only kj = ki,..,kN is 
    !! needed (the diagonal blocks where k'=k will be fully computed 
    !! but only the upper triangle will be needed in the end).

    !! Get corresponding q-point for ki,kj combination.
    ! K-point difference k_j-k_i on integer grid.
    iv(:) = ivknr(:,jknr) - ivknr(:,iknr)
    ! Map to reciprocal unit cell 
    iv(:) = modulo(iv(:), ngridq(:))
    ! Find corresponding q-point index (reduced)
    iqr = iqmapr(iv(1), iv(2), iv(3))
    ! Get corresponding vector in lattice coordinated
    vqr(:) = vqlr(:, iqr)
    ! Q-point (non-reduced)
    iq = iqmap(iv(1), iv(2), iv(3))
    ! Get lattice coordinates
    vq(:) = vql(:, iq)

    ! Local field effects size
    tq0 = tqgamma(iq)
    numgq = ngq(iq)

    allocate(igqmap(numgq)
    allocate(wfc(numgq, numgq))

    !! Find results (radial emat integrals and screened coulomb potential Fourier coefficients) 
    !! for a non reduced q-point with the help of the result 
    !! for the corresponding reduced q-point using symmetry operations.
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

    !-------------------------------!
    ! Resonant-Resonant Part        !
    !-------------------------------!
    allocate(cmoo(noo, numgq), cmuu(nuu, numgq))
    ! Calculate Moo and Muu for current non reduced q and k 
    call getpwesrr(iq, iknr, moo, muu)

    ! Combine indices for matrix elements of plane wave.
    ! Corresponds to j1 = subhamidx(io1,io2,no)
    j1 = 0
    ! Occupied k1 = k2+q = k+q
    do io2 = 1, no
      ! Occupied k2 = k
      do io1 = 1, no
        j1 = j1 + 1
        ! cmoo_j = M_o1o1, M_o2o1, ..., M_oNo1, M_o1o2, ..., M_oNoN
        cmoo(j1, :) = moo(io1, io2, :)
      end do
    end do
    j2 = 0
    ! Corresponds to j2 = subhamidx(iu1,iu2,no)
    ! Unoccupied (k+q)
    do iu2 = 1, nu
      ! Unoccupied (k)
      do iu1 = 1, nu
        j2 = j2 + 1
        ! cmuu_j = M_u1u1, M_u2u1, ..., M_uNu1, M_u1u2, ..., M_uNuN
        cmuu(j2, :) = muu(iu1, iu2, :)
      end do
    end do

    ! M_oioj -> M^*_oioj
    cmoo(:,:)=conjg(cmoo(:,:))

    ! Allocate helper array of dimension (#o*#o,#G) (same as cmoo)
    allocate(zm(noo,numgq))
    ! Calculate matrix elements of screened coulomb interaction scclit_{o_j1 o'_j1, u_j2 u'_j2}(q)
    ! zm = cmoo * wfc
    !   i.e. zm_{j,G} = \Sum_{G'} cmoo_{j,G'} tm_{G',G}
    !   i.e. zm_{o_j,o'_j}(G,q) = \Sum_{G'} M^*_{o_j,o'_j}(G',q) W(G',G, q)
    call zgemm('n', 'n', noo, numgq, numgq, zone, cmoo, noo, wfc, numgq, zzero, zm, noo)
    ! scclit = pref * zm * cmuu^T
    !   i.e. scclit(j1, j2) = \Sum_{G} zm_{j1,G} (cmuu^T)_{G,j2}
    !   i.e. scclit_{o_j1 o'_j2, u_j2 u'_j2} = \Sum{G,G'} M^*_{o_j1,o'_j1}(G,q) W(G,G', q) M_{u_j2 u'_j2}(G',q)
    call zgemm('n', 't', noo, nuu, numgq, pref, zm,&
      & noo, cmuu, nuu, zzero, scclit, noo)
    deallocate(zm)        

    ! Map back to individual band indices
    j2 = 0
    ! Unoccupied (k+q)
    do iu2 = 1, nu
      ! Unoccupied (k)
      do iu1 = 1, nu
        j2 = j2 + 1
        j1 = 0
        ! Occupied (k+q)
        do io2 = 1, no
          ! Occupied (k)
          do io1 = 1, no
            j1 = j1 + 1
            ! scclit_{o_j1 o'_j1, u_j2 u'_j2} -> sccli_{o_j1 u_j2, o'_j1 u'_j2}
            sccli(io1, iu1, io2, iu2) = scclit(j1, j2)
          end do
        end do
      end do
    end do

    ! Analyze BSE diagonal
    if(iknr .eq. jknr) then
      ! Selected occupied
      do io1 = 1, no
        ! Selected unoccupied
        do iu1 = 1, nu
          zt1 = sccli(io1, iu1, io1, iu1)
          scclid(io1, iu1) = zt1
          bsedt(1, rank) = min(dble(bsedt(1, rank)), dble(zt1))
          bsedt(2, rank) = max(dble(bsedt(2, rank)), dble(zt1))
          bsedt(3, rank) = bsedt(3, rank) + zt1 / (no*nu)
        end do
      end do
    end if

    ! Parallel write
    call putbsemat('SCCLI.OUT', 77, sccli, ikkp, iknr, jknr,&
      & iq, iqr, no, nu, no, nu)

    deallocate(moo, muu, cmoo, cmuu)
    
    !-------------------------------!
    ! Resonant-Anit-Resonant Part   !
    !-------------------------------!
    if(fcoup) then

      allocate(cmou(nou, numgq), cmuo(nuo, numgq))
      ! Calculate Mou and Muo for current non reduced q and k 
      call getpwesra(iq, iknr, mou, muo)

      ! Combine indices for matrix elements of plane wave.
      j1 = 0
      ! Unoccupied 
      do iu1 = 1, nu ! k+q
        ! Occupied 
        do io1 = 1, no ! k
          j1 = j1 + 1
          ! cmou_j = M_o1u1, M_o2u1, ..., M_oNu1, M_o1u2, ..., M_oNuM
          cmou(j1, :) = mou(io1, iu1, :)
        end do
      end do
      j2 = 0
      ! Unoccupied 
      do iu2 = 1, nu
        ! Occupied 
        do io2 = 1, no
          j2 = j2 + 1
          ! cmuo_j = M_u1o1, M_u1o2, ..., M_u1oN, M_u2o1, ..., M_uMoN
          cmuo(j2, :) = muo(iu2, io2, :)
        end do
      end do

!!!!!!! CONTINUE HERE !!!!!!!

      ! M_oioj -> M^*_oioj
      cmoo(:,:)=conjg(cmoo(:,:))

      ! Allocate helper array of dimension (#o*#o,#G) (same as cmoo)
      allocate(zm(noo,numgq))
      ! Calculate matrix elements of screened coulomb interaction scclit_{o_j1 o'_j1, u_j2 u'_j2}(q)
      ! zm = cmoo * wfc
      !   i.e. zm_{j,G} = \Sum_{G'} cmoo_{j,G'} tm_{G',G}
      !   i.e. zm_{o_j,o'_j}(G,q) = \Sum_{G'} M^*_{o_j,o'_j}(G',q) W(G',G, q)
      call zgemm('n', 'n', noo, numgq, numgq, zone, cmoo, noo, wfc, numgq, zzero, zm, noo)
      ! scclit = pref * zm * cmuu^T
      !   i.e. scclit(j1, j2) = \Sum_{G} zm_{j1,G} (cmuu^T)_{G,j2}
      !   i.e. scclit_{o_j1 o'_j2, u_j2 u'_j2} = \Sum{G,G'} M^*_{o_j1,o'_j1}(G,q) W(G,G', q) M_{u_j2 u'_j2}(G',q)
      call zgemm('n', 't', noo, nuu, numgq, pref, zm,&
        & noo, cmuu, nuu, zzero, scclit, noo)
      deallocate(zm)        

      ! Map back to individual band indices
      j2 = 0
      ! Unoccupied (k+q)
      do iu2 = 1, nu
        ! Unoccupied (k)
        do iu1 = 1, nu
          j2 = j2 + 1
          j1 = 0
          ! Occupied (k+q)
          do io2 = 1, no
            ! Occupied (k)
            do io1 = 1, no
              j1 = j1 + 1
              ! scclit_{o_j1 o'_j1, u_j2 u'_j2} -> sccli_{o_j1 u_j2, o'_j1 u'_j2}
              sccli(io1, iu1, io2, iu2) = scclit(j1, j2)
            end do
          end do
        end do
      end do
    end if

    
    deallocate(igqmap)
    deallocate(wfc)

  ! End loop over(k,kp)-pairs
  end do kkploop

  ! Deallocate helper array (independent of number of G)
  deallocate(scclit)

  call barrier

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(procs,3,zbuf=bsedt)

  ! BSE kernel diagonal parameters
  bsedl = minval(dble(bsedt(1, :)))
  bsedu = maxval(dble(bsedt(2, :)))
  bsedd = bsedu - bsedl
  bsed = sum(bsedt(3, :)) / nkptnr
  deallocate(bsedt, scclid)

  ! Write BSE kernel diagonal parameters
  if(rank .eq. 0) call putbsediag('BSEDIAG.OUT')

  call findgntn0_clear

  if(rank .eq. 0) then
    write(unitout, '("Info(scrcoulint): Screened coulomb interaction&
      & finished")')
  end if

  contains

    subroutine getpwesrr(iqnr, iknr, moo, muu)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: moo(:,:,:), muu(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the oo plane wave elements
      ematbc%n1=no
      ematbc%il1=bcouabs%il1
      ematbc%iu1=bcouabs%iu1
      ematbc%n2=no
      ematbc%il2=bcouabs%il1
      ematbc%iu2=bcouabs%iu1
      ! Allocate space for M_{o1o2,G} at fixed (k, q)
      if(allocated(moo)) deallocate(moo)
      allocate(moo(no,no,numgq))
      ! Calculate M_{o1o2,G} at fixed (k, q)
      call ematqk(iqnr, iknr, moo, ematbc)

      !! Calculate the uu plane wave elements
      ematbc%n1=nu
      ematbc%il1=bcouabs%il2
      ematbc%iu1=bcouabs%iu2
      ematbc%n2=nu
      ematbc%il2=bcouabs%il2
      ematbc%iu2=bcouabs%iu2
      ! Allocate space for M_{u1u2,G} at fixed (k, q)
      if(allocated(muu)) deallocate(muu)
      allocate(muu(nu,nu,numgq))
      ! Calculate M_{o1o2,G} at fixed (k, q)
      call ematqk(iqnr, iknr, muu, ematbc)

    end subroutine getpwesrr

    subroutine getpwesra(iqnr, iknr, mou, muo)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: mou(:,:,:), muo(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the o-u plane wave elements
      ematbc%n1=no
      ematbc%il1=bcouabs%il1
      ematbc%iu1=bcouabs%iu1
      ematbc%n2=nu
      ematbc%il2=bcouabs%il2
      ematbc%iu2=bcouabs%iu2
      ! Allocate space for M_{ou,G} at fixed (k, q)
      if(allocated(mou)) deallocate(mou)
      allocate(mou(no,nu,numgq))
      ! Calculate M_{ou,G} at fixed (k, q)
      call ematqk(iqnr, iknr, mou, ematbc)

      !! Calculate the u-o plane wave elements
      ematbc%n1=nu
      ematbc%il1=bcouabs%il2
      ematbc%iu1=bcouabs%iu2
      ematbc%n2=no
      ematbc%il2=bcouabs%il1
      ematbc%iu2=bcouabs%iu1
      ! Allocate space for M_{uo,G} at fixed (k, q)
      if(allocated(muo)) deallocate(muo)
      allocate(muo(nu,no,numgq))
      ! Calculate M_{uo,G} at fixed (k, q)
      call ematqk(iqnr, iknr, muo, ematbc)

    end subroutine getpwesra

end subroutine b_scrcoulint
!EOC

