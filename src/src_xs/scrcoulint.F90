! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: scrcoulint
! !INTERFACE:
subroutine scrcoulint
! !USES:
  use modinput, only: input
  use modmpi
  use mod_misc, only: task
  use mod_constants, only: zzero, zone, fourpi
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: nqpt, iqmap, ngridq, vql
  use mod_kpoint, only: nkptnr, ivknr
  use mod_lattice, only: omega
  use mod_symmetry, only: maxsymcrys
  use modxs, only: sta1, sto1, sta2, sto2,&
                 & nst1, nst2, nst3, nst4,&
                 & xsgnt, unitout, istocc0, istocc,&
                 & istunocc0, istunocc, isto0, isto,&
                 & istu0, istu, ksgap, ngqmax,&
                 & nqptr, qpari, qparf, ivqr,&
                 & ngq, ppari, pparf, iqmapr,&
                 & vqlr, vgqc, dielten, xiou,&
                 & xiuo, bsedl, bsedu, bsedd,&
                 & bsed, ikmapikq
  use m_xsgauntgen
  use m_findgntn0
  use m_writevars
  use m_genfilname
  use m_getunit
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
  complex(8) :: zt1,prefactor
  complex(8), allocatable :: scclit(:, :), sccli(:, :, :, :), scclid(:, :)
  complex(8), allocatable :: scieffg(:, :, :), tm(:, :), tmi(:, :), bsedt(:, :),zm(:,:)
  complex(8), allocatable :: phf(:, :), emat12(:, :), emat34(:, :)
  real(8), parameter :: epsortho = 1.d-12
  real(8) :: vqr(3), vq(3), t1, ta
  integer :: ikkp, iknr, jknr, iqr, iq, iqrnr, jsym, jsymi, igq1, n, reclen
  integer :: nsc, iv(3), ivgsym(3), j1, j2, nkkp
  integer :: ist1, ist2, ist3, ist4, nst12, nst34, nst13, nst24
  integer :: rnst1, rnst2, rnst3, rnst4
  integer :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer, allocatable :: igqmap(:)
  logical :: tq0, tphf

  ! External functions
  integer, external :: idxkkp
  logical, external :: tqgamma

  !---------------!
  !   main part   !
  !---------------!

  ! emattype 2 selects o-o and u-u combinations
  input%xs%emattype = 2

  call init0
  call init1
  call init2

  ! Set the range of valence/core and conduction states to use
  ! Lowest occupied state (absolute index)
  sta1 = input%xs%bse%nstlbse(1) 
  ! Highest occupied state (absolute index)
  sto1 = input%xs%bse%nstlbse(2)
  ! Lowest unoccupied state (counted from first unoccupied state)
  sta2 = input%xs%bse%nstlbse(3) 
  ! Highest unoccupied state (counted from first unoccupied state)
  sto2 = input%xs%bse%nstlbse(4)      

  ! Number of occupied states
  rnst1 = sto1-sta1+1
  rnst2 = sto1-sta1+1
  ! Number of unoccupied states
  rnst3 = sto2-sta2+1
  rnst4 = sto2-sta2+1

  ! Read Fermi energy from file
  call readfermi

  ! Save variables for the gamma q-point
  call xssave0

  ! Generate gaunt coefficients
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
    & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
  ! Find indices for non-zero gaunt coefficients
  call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
    & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

  write(unitout, '(a,3i8)') 'Info(' // thisnam // '):&
    & Gaunt coefficients generated within lmax values:', input%groundstate%lmaxapw,&
    & input%xs%lmaxemat, input%groundstate%lmaxapw

  write(unitout, '(a, i6)') 'Info(' // thisnam // '): number of q-points: ', nqpt

  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Find occupation bounds for k and k+q but q=0
  call findocclims(0, ikmapikq(:,1), istocc0, istunocc0, isto0, isto, istu0, istu)
  istunocc = istunocc0
  istocc = istocc0
  ! Only for systems with a gap in energy
  if( .not. ksgap) then
    write(*,*)
    write(*,'("Warning(",a,"): There is no ks-gap present")') trim(thisnam)
    write(*,*)
  end if

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

  ! Set nstX, istlX, istuX variables for X: 1=o, 2=o, 3=u, 4=u
  call ematbdcmbs(input%xs%emattype)

  ! Number of oo-combinations
  nst12 = rnst1 * rnst2
  ! Number of uu-combinations
  nst34 = rnst3 * rnst4
  ! Number of ou-combinations
  nst13 = rnst1 * rnst3
  ! Number of ou-combinations
  nst24 = rnst2 * rnst4

  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if

  ! Allocate local arrays for screened coulomb interaction
  allocate(phf(ngqmax, ngqmax))
  allocate(sccli(rnst1, rnst3, rnst2, rnst4), scclid(rnst1, rnst3))
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  sccli(:, :, :, :) = zzero
  scieffg(:, :, :) = zzero

  ! Set file extension
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  !-----------------------------------!
  !     Loop over reduced q-points    !
  !-----------------------------------!
  ! Parallelize over reduced q-point set
  call genparidxran('q', nqptr)

  do iqr = qpari, qparf

    call chkpt(3, (/ task, 1, iqr /),&
      & 'task,sub,reduced q-point; generate effective screened Coulomb potential')

    ! Locate reduced q-point in non-reduced set
    iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))

    ! Get number of G+q vectors for current q
    n = ngq(iqrnr)

    ! Calculate effective screened coulomb interaction
    ! by inverting the symmetrized RPA dielectric matrix and multiplying
    ! it with v^{1/2} from both sides.
    call genscclieff(iqr, iqrnr, ngqmax, n, scieffg(1,1,iqr))

    ! Generate radial integrals for matrix elements of plane wave
    call putematrad(iqr, iqrnr)

  end do

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(nqptr,rlen=ngqmax*ngqmax,zbuf=scieffg)
  call barrier

  ! Information on size of output file
  nkkp = (nkptnr*(nkptnr+1)) / 2
  inquire(iolength=reclen) ikkp, iknr, jknr, iq, iqr, nst1, nst2, nst3, nst4, sccli
  write(unitout,*)
  write(unitout, '(a,f12.3)') 'File size for screened coulomb interaction (Gb):',&
    & reclen * nkkp / 1024.d0 ** 3
  write(unitout,*)

  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!

  call genparidxran('p', nkkp)

  allocate(bsedt(3, 0:procs-1))
  bsedt(1, :) = 1.d8
  bsedt(2, :) = -1.d8
  bsedt(3, :) = zzero

  ! Loop over combinations of non-reduced k-point combinations
  kkploop: do ikkp = ppari, pparf

    call timesec(ta)
    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task, sub,(k,kp)-pair; direct term of BSE Hamiltonian')

    ! Get individual k-point indices from combined kk' index.
    !   iknr runs from 1 to nkptnr, jknr from iknr to nkptnr
    call kkpmap(ikkp, nkptnr, iknr, jknr)
    !! Note: If the screened coulomb interaction matrix 
    !! has the elements W_{i,j} and the indices enumerate the
    !! states according to
    !! i = o1u1k1, o1u1k2, ..., o1u1kN, o2u1k1, ..., o2u1kN, ...,
    !!     oMu1kN, oMu2k1, ..., oMuMkN
    !! then because of W_{j,i} = W^*_{i,j} only kj = ki,..,kN is 
    !! needed.

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
    n = ngq(iq)

    allocate(igqmap(n), emat12(nst12, n), emat34(nst34, n))
    allocate(tm(n, n), tmi(n, n))
    allocate(scclit(nst12, nst34))


    !! Find results (radial emat integrals and screened coulomb potential) 
    !! for a non reduced q-point with the help of the result 
    !! for the corresponding reduced q-point using symmetry operations.
    !!<--
    ! Find symmetry operations that reduce the q-point to the irreducible
    ! part of the Brillouin zone
    call findsymeqiv(input%xs%bse%fbzq, vq, vqr, nsc, sc, ivgsc)
    ! Find the map that rotates the G-vectors
    call findgqmap(iq, iqr, nsc, sc, ivgsc, n, jsym, jsymi, ivgsym, igqmap)
    ! Generate phase factor for dielectric matrix due to non-primitive
    ! translations
    call genphasedm(iq, jsym, ngqmax, n, phf, tphf)
    ! Get radial integrals
    call getematrad(iqr, iq)
    ! Rotate radial integrals
    call rotematrad(n, igqmap)
    ! Rotate inverse of screening, coulomb potential and radial integrals
    tmi(:,:) = phf(:n, :n) * scieffg(igqmap, igqmap, iqr)
    !!-->

    ! Calculate matrix elements of the plane wave
    ! emattype = 2 corresponds to 12=oo and 34=uu
    input%xs%emattype = 2
    call ematbdcmbs(input%xs%emattype)
    ! Allocate arrays used in ematqk
    call ematqalloc
    ! Call wrapper function for ematqk
    !   Calculates o-o/u-u plane wave elements 
    !   for iknr at q=jk-ik and stores them in xiou/xiuo 
    call ematqk1(iq, iknr)
    ! Reset 12=oo and 34=uu (changed in ematqk1)
    input%xs%emattype = 2
    call ematbdcmbs(input%xs%emattype)
    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task,sub,(k,kp)-pair; direct term of BSE Hamiltonian')

    ! Select screening level (default: full)
    tm(:, :) = zzero
    select case(trim(input%xs%screening%screentype))
      case('longrange')
        ! Constant screening(q=0 average tensor)
        forall(igq1=1:n)
          tm(igq1, igq1) =&
            & fourpi * dot_product(vgqc(:, igq1, iq), matmul(dielten, vgqc(:, igq1, iq)))
        end forall
      case('diag')
        ! Only diagonal of screening
        forall(igq1=1:n)
          tm(igq1, igq1) = tmi(igq1, igq1)
        end forall
      case('full')
        ! Full screening tm is set to the phase factor corrected
        ! W_{GG'}(q,\omega=0)
        tm(:, :) = tmi(:, :)
    end select

    ! Combine indices for matrix elements of plane wave
    j1 = 0
    ! Occupied
    do ist2 = sta1, sto1
      ! Occupied
      do ist1 = sta1, sto1
        j1 = j1 + 1
        ! emat12_j = M_o1o1, M_o2o1, ..., M_oNo1, M_o1o2, ..., M_oNoN
        emat12(j1, :) = xiou(ist1-sta1+1, ist2-sta1+1, :)
      end do
    end do
    j2 = 0
    ! Unoccupied
    do ist4 = sta2, sto2
      ! Unoccupied
      do ist3 = sta2, sto2
        j2 = j2 + 1
        ! emat34_j = M_u1u1, M_u2u1, ..., M_uNu1, M_u1u2, ..., M_uNuN
        emat34(j2, :) = xiuo(ist3-sta2+1, ist4-sta2+1, :)
      end do
    end do

    ! Matrix elements of direct term (as in BSE-code of Peter and
    ! In the self-documentation of Andrea Marini)

    prefactor=1.0d0/(omega*dble(nkptnr))

    ! M_oioj -> M^*_oioj
    emat12(:,:)=conjg(emat12(:,:))

    ! Allocate helper array of dimension (#o*#o,#G) (same as emat12)
    allocate(zm(nst12,n))

    ! Calculate matrix elements of screened coulomb interaction scclit_{o_j1 o'_j1, u_j2 u'_j2}(q)
    ! zm = emat12 * tm
    !   i.e. zm_{j,G} = \Sum_{G'} emat12_{j,G'} tm_{G',G}
    !   i.e. zm_{o_j,o'_j}(G,q) = \Sum_{G'} M^*_{o_j,o'_j}(G',q) W(G',G, q)
    call zgemm('n', 'n', nst12, n, n, zone, emat12, nst12, tm, n, zzero, zm, nst12)
    ! scclit = prefactor * zm * emat34^T
    !   i.e. scclit(j1, j2) = \Sum_{G} zm_{j1,G} (emat34^T)_{G,j2}
    !   i.e. scclit_{o_j1 o'_j2, u_j2 u'_j2} = \Sum{G,G'} M^*_{o_j1,o'_j1}(G,q) W(G,G', q) M_{u_j2 u'_j2}(G',q)
    call zgemm('n', 't', nst12, nst34, n, prefactor, zm,&
      & nst12, emat34, nst34, zzero, scclit, nst12)
    deallocate(zm)        

    ! Map back to individual band indices
    j2 = 0
    ! Unoccupied
    do ist4 = 1, rnst4
      ! Unoccupied
      do ist3 = 1, rnst3
        j2 = j2 + 1
        j1 = 0
        ! Occupied
        do ist2 = 1, rnst2
          ! Occupied
          do ist1 = 1, rnst1
            j1 = j1 + 1
            ! scclit_{o_j1 o'_j1, u_j2 u'_j2} -> sccli_{o_j1 u_j2, o'_j1 u'_j2}
            sccli(ist1, ist3, ist2, ist4) = scclit(j1, j2)
          end do
        end do
      end do
    end do

    ! Analyze BSE diagonal
    if(iknr .eq. jknr) then
      ! Selected occupied
      do ist1 = 1, rnst1
        ! Selected unoccupied
        do ist3 = 1, rnst3
          zt1 = sccli(ist1, ist3, ist1, ist3)
          scclid(ist1, ist3) = zt1
          t1 = dble(zt1)
          bsedt(1, rank) = min(dble(bsedt(1, rank)), t1)
          bsedt(2, rank) = max(dble(bsedt(2, rank)), t1)
          bsedt(3, rank) = bsedt(3, rank) + zt1 / (rnst1*rnst3)
        end do
      end do
    end if

    ! Parallel write
    call putbsemat('SCCLI.OUT', 77, sccli, ikkp, iknr, jknr,&
      & iq, iqr, rnst1, rnst3, rnst2, rnst4)

    deallocate(igqmap, emat12, emat34)
    deallocate(tm, tmi)
    deallocate(scclit)

  ! End loop over(k,kp)-pairs
  end do kkploop

  call barrier

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(procs,rlen=3,zbuf=bsedt)

  ! BSE kernel diagonal parameters
  bsedl = minval(dble(bsedt(1, :)))
  bsedu = maxval(dble(bsedt(2, :)))
  bsedd = bsedu - bsedl
  bsed = sum(bsedt(3, :)) / nkptnr
  deallocate(bsedt, scclid)

  ! Write BSE kernel diagonal parameters
  if(rank .eq. 0) call putbsediag('BSEDIAG.OUT')

  call findgntn0_clear

  write(unitout, '("Info(scrcoulint): Screened coulomb interaction&
    & finished")')

end subroutine scrcoulint
!EOC

