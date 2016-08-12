! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: scrcoulint
! !INTERFACE:
subroutine scrcoulint
! !USES:
  use modinput, only: input
  use modmpi, only: procs, rank, mpi_allgatherv_ifc, barrier
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
                 & bsed
  !use summations
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
  character(256) :: fnsccli, fnscreeninv
  complex(8) :: zt1,prefactor
  complex(8), allocatable :: scclit(:, :), sccli(:, :, :, :), scclid(:, :)
  complex(8), allocatable :: scieffg(:, :, :), tm(:, :), tmi(:, :), bsedt(:, :),zm(:,:)
  complex(8), allocatable :: phf(:, :), emat12(:, :), emat34(:, :)
  real(8), parameter :: epsortho = 1.d-12
  real(8) :: vqr(3), vq(3), t1, ta, tb
  integer :: ikkp, iknr, jknr, iqr, iq, iqrnr, jsym, jsymi, igq1, igq2, n, recl, un
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

  input%xs%emattype = 2
  call init0
  call init1
  call init2

  ! Set the range of valence/core and conduction states to use
  sta1 = input%xs%bse%nstlbse(1)
  sto1 = input%xs%bse%nstlbse(2)
  sta2 = input%xs%bse%nstlbse(3)
  sto2 = input%xs%bse%nstlbse(4)
  rnst1 = sto1-sta1+1
  rnst2 = sto1-sta1+1
  rnst3 = sto2-sta2+1
  rnst4 = sto2-sta2+1

  ! Read fermi energy from file
  call readfermi

  ! Save variables for the gamma q-point
  call xssave0

  ! Generate gaunt coefficients
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
    & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))

  ! Find indices for non-zero gaunt coefficients
  call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
    & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

  write(unitout, '(a,3i8)') 'info(' // thisnam // '):&
    & Gaunt coefficients generated within lmax values:', input%groundstate%lmaxapw,&
    & input%xs%lmaxemat, input%groundstate%lmaxapw

  write(unitout, '(a, i6)') 'info(' // thisnam // '): number of q-points: ', nqpt

  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  call findocclims(0, istocc0, istocc, istunocc0, istunocc, isto0, isto, istu0, istu)

  ! Only for systems with a gap in energy
  if( .not. ksgap) then
    write(*,*)
    write(*,'("warning(",a,"): there is no ks-gap present")') trim(thisnam)
    write(*,*)
  end if

  ! Check number of empty states
  if(input%xs%screening%nempty .lt. input%groundstate%nempty) then
    write(*,*)
    write(*, '("Error(",a,"): Too few empty states in&
      & screening eigenvector file - the screening should&
      & include many empty states (bse/screening)", 2i8)')&
      & trim(thisnam), input%groundstate%nempty, input%xs%screening%nempty
    write(*,*)
    call terminate
  end if

  call ematbdcmbs(input%xs%emattype)
  nst12 = rnst1 * rnst2
  nst34 = rnst3 * rnst4
  nst13 = rnst1 * rnst3
  nst24 = rnst2 * rnst4

  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if

  ! Local arrays
  allocate(phf(ngqmax, ngqmax))
  allocate(sccli(rnst1, rnst3, rnst2, rnst4), scclid(rnst1, rnst3))
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  sccli(:, :, :, :) = zzero
  scieffg(:, :, :) = zzero

  ! Set file extension
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  !-----------------------------------!
  !     loop over reduced q-points    !
  !-----------------------------------!
  call getunit(un)
  call genparidxran('q', nqptr)

  do iqr = qpari, qparf

    call chkpt(3, (/ task, 1, iqr /),&
      & 'task,sub,reduced q-point; generate effective screened coulomb potential')

    ! Locate reduced q-point in non-reduced set
    iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))
    n = ngq(iqrnr)

    ! Calculate effective screened coulomb interaction
    call genscclieff(iqr, ngqmax, n, scieffg(1,1,iqr))

    ! Generate radial integrals for matrix elements of plane wave
    call putematrad(iqr, iqrnr)

  end do

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(nqptr,ngqmax*ngqmax,zbuf=scieffg)
  call barrier

  ! Information on size of output file
  nkkp = (nkptnr*(nkptnr+1)) / 2
  inquire(iolength=recl) ikkp, iknr, jknr, iq, iqr, nst1, nst2, nst3, nst4, sccli
  write(unitout,*)
  write(unitout, '(a,f12.3)') 'File size for screened coulomb interaction (gb):',&
    & recl * nkkp / 1024.d0 ** 3
  write(unitout,*)

  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!
  call genfilname(basename='sccli', asc=.true., filnam=fnsccli)
  call getunit(un)

  if(rank .eq. 0) open(un, file=trim(fnsccli),&
    & form='formatted', action='write', status='replace')

  call genparidxran('p', nkkp)

  allocate(bsedt(3, 0:procs-1))
  bsedt(1, :) = 1.d8
  bsedt(2, :) = -1.d8
  bsedt(3, :) = zzero

  ! Loop over combinations of k-points
  kkploop: do ikkp = ppari, pparf

    call timesec(ta)
    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task, sub,(k,kp)-pair; direct term of bse hamiltonian')

    call kkpmap(ikkp, nkptnr, iknr, jknr)
    ! K-point difference
    iv(:) = ivknr(:,jknr) - ivknr(:,iknr)
    iv(:) = modulo(iv(:), ngridq(:))
    ! Q-point(reduced)
    iqr = iqmapr(iv(1), iv(2), iv(3))
    vqr(:) = vqlr(:, iqr)
    ! Q-point(non-reduced)
    iq = iqmap(iv(1), iv(2), iv(3))
    vq(:) = vql(:, iq)
    ! Local field effects size
    tq0 = tqgamma(iq)
    n = ngq(iq)

    allocate(igqmap(n), emat12(nst12, n), emat34(nst34, n))
    allocate(tm(n, n), tmi(n, n))
    allocate(scclit(nst12, nst34))

    ! Find symmetry operations that reduce the q-point to the irreducible
    ! part of the brillouin zone
    call findsymeqiv(input%xs%bse%fbzq, vq, vqr, nsc, sc, ivgsc)

    ! Find the map that rotates the g-vectors
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

    ! Calculate matrix elements of the plane wave
    input%xs%emattype = 2
    call ematbdcmbs(input%xs%emattype)
    call ematqalloc
    call ematqk1(iq, iknr)

    input%xs%emattype = 2
    call ematbdcmbs(input%xs%emattype)
    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task,sub,(k,kp)-pair; direct term of bse hamiltonian')

    ! Select screening level
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
        ! Full screening
        tm(:, :) = tmi(:, :)
    end select

    ! Combine indices for matrix elements of plane wave
    j1 = 0
    do ist2 = sta1, sto1
      do ist1 = sta1, sto1
        j1 = j1 + 1
        emat12(j1, :) = xiou(ist1-sta1+1, ist2-sta1+1, :)
      end do
    end do
    j2 = 0
    do ist4 = sta2, sto2
      do ist3 = sta2, sto2
        j2 = j2 + 1
        emat34(j2, :) = xiuo(ist3-sta2+1, ist4-sta2+1, :)
      end do
    end do

    ! Matrix elements of direct term (as in bse-code of peter and
    ! In the self-documentation of andrea marini)

    prefactor=1.0d0/(omega*dble(nkptnr))
    emat12(:,:)=conjg(emat12(:,:))
    allocate(zm(nst12,n))
    call zgemm('n', 'n', nst12, n, n, zone, emat12, nst12, tm, n, zzero, zm, nst12)
    call zgemm('n', 't', nst12, nst34, n, prefactor, zm,&
      & nst12, emat34, nst34, zzero, scclit, nst12)
    deallocate(zm)        

    ! Map back to individual band indices
    j2 = 0
    do ist4 = 1, rnst4
      do ist3 = 1, rnst3
        j2 = j2 + 1
        j1 = 0
        do ist2 = 1, rnst2
          do ist1 = 1, rnst1
            j1 = j1 + 1
            sccli(ist1, ist3, ist2, ist4) = scclit(j1, j2)
          end do
        end do
      end do
    end do

    ! Analyze bse diagonal
    if(iknr .eq. jknr) then
      do ist1 = 1, rnst1
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
    call putbsemat('SCCLI.OUT', sccli, ikkp, iknr, jknr,&
      & iq, iqr, rnst1, rnst3, rnst2, rnst4)

    deallocate(igqmap, emat12, emat34)
    deallocate(tm, tmi)
    deallocate(scclit)

  ! End loop over(k,kp)-pairs
  end do kkploop

  if(rank .eq. 0) write(un, '("# ikkp, iknr,ist1,ist3, jknr,ist2,ist4,&
    &    re(w),            im(w),             |w|^2,           ang/pi")')

  if(rank .eq. 0) close(un)

  call barrier

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(procs,3,zbuf=bsedt)

  ! Bse kernel diagonal parameters
  bsedl = minval(dble(bsedt(1, :)))
  bsedu = maxval(dble(bsedt(2, :)))
  bsedd = bsedu - bsedl
  bsed = sum(bsedt(3, :)) / nkptnr
  deallocate(bsedt, scclid)

  ! Write bse kernel diagonal parameters
  if(rank .eq. 0) call putbsediag('BSEDIAG.OUT')

  call findgntn0_clear

  write(unitout, '("Info(scrcoulint): Screened coulomb interaction&
    & finished")')

end subroutine scrcoulint
!eoc

