! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: exccoulint
! !INTERFACE:
subroutine exccoulint
! !USES:
  use mod_constants, only: zone, zzero
  use mod_misc, only: task
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: nqpt, iqmap
  use mod_kpoint, only: nkptnr, ivknr
  use mod_lattice, only: omega
  use modinput, only: input
  use modmpi, only: rank, barrier, mpi_allgatherv_ifc
  use modxs, only: sta1, sto1, sta2, sto2,&
                 & xsgnt, unitout, istocc0, istocc,&
                 & istunocc0, istunocc, isto0, isto,&
                 & istu0, istu, ksgap, ngq,&
                 & kpari, kparf, xiou, xiuo,&
                 & ppari, pparf, iqmapr, nst1,&
                 & nst2, qvkloff
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  use m_getunit
! !DESCRIPTION:
!   Calculates the exchange term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!   level for the treatment of core excitations(using local orbitals).
!   October 2010 (Weine Olovsson)
!EOP
!BOC      

  implicit none

  ! Local variables
  character(*), parameter :: thisnam = 'exccoulint'
  character(256) :: fnexcli
  integer, parameter :: iqmt = 1
  real(8), parameter :: epsortho = 1.d-12
  integer :: iknr, jknr, iqr, iq, igq1, n, un
  integer :: iv(3), j1, j2
  integer :: ist1, ist2, ist3, ist4
  integer :: nst12, nst34, nst13, nst24, ikkp, nkkp
  integer :: rnst1, rnst2, rnst3, rnst4
  real(8), allocatable :: potcl(:)
  complex(8), allocatable :: exclit(:, :), excli(:, :, :, :)
  complex(8), allocatable :: emat12(:, :), emat34(:, :)
  complex(8), allocatable :: emat12k(:, :, :, :)

  !---------------!
  !   main part   !
  !---------------!
  input%xs%emattype = 1
  call init0
  call init1
  call init2

  ! Set the range of valence/core and conduction states to use
  sta1 = input%xs%bse%nstlbse(1)
  sto1 = input%xs%bse%nstlbse(2)
  sta2 = input%xs%bse%nstlbse(3)
  sto2 = input%xs%bse%nstlbse(4)      
  rnst1 = sto1-sta1+1
  rnst2 = sto2-sta2+1
  rnst3 = sto2-sta2+1
  rnst4 = sto1-sta1+1

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

  write(unitout, '(a,3i8)') 'Info(' // thisnam // '):&
    & Gaunt coefficients generated within lmax values:',&
    & input%groundstate%lmaxapw, input%xs%lmaxemat, input%groundstate%lmaxapw

  write(unitout, '(a, i6)') 'Info(' // thisnam // '):&
    & Number of q-points: ', nqpt

  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
  call findocclims(0, istocc0, istocc, istunocc0, istunocc, isto0, isto, istu0, istu)

  ! Only for systems with a gap in energy
  if( .not. ksgap) then
    write(*,*)
    write(*, '("Warning(",a,"): there is no ks-gap present")') trim(thisnam)
    write(*,*)
  end if

  ! Check number of empty states
  if(input%xs%screening%nempty .lt. input%groundstate%nempty) then
    write(*,*)
    write(*, '("Error(",a,"): too few empty states in screening eigenvector file&
      & - the screening should include many empty states (bse/screening)", 2i8)')&
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

  n = ngq(iqmt)

  call ematrad(iqmt)

  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  allocate(potcl(n))
  allocate(excli(rnst1, rnst2, rnst1, rnst2))
  allocate(exclit(nst12, nst34))
  allocate(emat12k(nst1, nst2, n, nkptnr))
  allocate(emat12(nst12, n), emat34(nst34, n))

  potcl(:) = 0.d0
  excli(:, :, :, :) = zzero

  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
  call genparidxran('k', nkptnr)
  call init1offs(qvkloff(1, iqmt))
  call ematqalloc

  do iknr = kpari, kparf
    call chkpt(3, (/ task, 1, iknr /),&
      & 'task,sub,k-point; matrix elements of plane wave')

    ! Matrix elements for k and q=0
    call ematqk1(iqmt, iknr)
    emat12k(:, :, :, iknr) = xiou(:, :, :)
    deallocate(xiou, xiuo)

  end do

  ! Communicate array-parts wrt. k-points
  call mpi_allgatherv_ifc(nkptnr,nst1*nst2*n,zbuf=emat12k)

  input%xs%emattype = 1
  call ematbdcmbs(input%xs%emattype)

  !-------------------------------!
  !     loop over(k,kp) pairs     !
  !-------------------------------!
  nkkp = (nkptnr*(nkptnr+1)) / 2
  call genparidxran('p', nkkp)
  call genfilname(basename='EXCLI', asc=.true., filnam=fnexcli)
  call getunit(un)

  if(rank .eq. 0) open(un, file=trim(fnexcli),&
    & form='formatted', action='write', status='replace')

  kkp: do ikkp = ppari, pparf

    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task,sub,(k,kp)-pair; exchange term of bse hamiltonian')

    call kkpmap(ikkp, nkptnr, iknr, jknr)
    iv(:) = ivknr(:, jknr) - ivknr(:, iknr)
    iv(:) = modulo(iv(:), input%groundstate%ngridk(:))

    ! Q-point(reduced)
    iqr = iqmapr(iv(1), iv(2), iv(3))
    ! Q-point(non-reduced)
    iq = iqmap(iv(1), iv(2), iv(3))

    ! Set g=0 term of coulomb potential to zero [Ambegaokar-Kohn]
    potcl(1) = 0.d0

    ! Set up coulomb potential
    do igq1 = 2, n
      call genwiqggp(0, iqmt, igq1, igq1, potcl(igq1))
    end do

    call genfilname(dotext='_SCR.OUT', setfilext=.true.)

    j1 = 0
    do ist2 = sta2, sto2
      do ist1 = sta1, sto1
        j1 = j1 + 1
        emat12(j1, :) = emat12k(ist1-sta1+1, ist2-sta2+1, :, iknr)
      end do
    end do
    j2 = 0
    do ist4 = sta2, sto2
      do ist3 = sta1, sto1
        j2 = j2 + 1
        emat34(j2, :) = emat12k(ist3-sta1+1, ist4-sta2+1, :, jknr) * potcl(:)
      end do
    end do

    ! Calculate exchange matrix elements: v_{1234} = m_{12}^* m_{34}^t
    emat12 = conjg(emat12)
    call zgemm('n', 't', nst12, nst12, n, zone/omega/nkptnr,&
      & emat12, nst12, emat34, nst12, zzero, exclit, nst12)

    ! Map back to individual band indices
    j2 = 0
    do ist4 = 1, rnst2
      do ist3 = 1, rnst1
        j2 = j2 + 1
        j1 = 0
        do ist2 = 1, rnst2
          do ist1 = 1, rnst1
            j1 = j1 + 1
            excli(ist1, ist2, ist3, ist4) = exclit(j1, j2)
          end do
        end do
      end do
    end do

    ! Parallel write
    call putbsemat('EXCLI.OUT', excli, ikkp, iknr, jknr,&
      & iq, iqr, rnst1, rnst2, rnst4, rnst3)
    call genfilname(dotext='_SCI.OUT', setfilext=.true.)

  ! End loop over(k,kp) pairs
  end do kkp

  if(rank .eq. 0) then
    write(un, '("# ikkp, iknr,ist1,ist3, jknr,ist2,ist4,&
      &    re(v),            im(v),             |v|^2")')
    close(un)
  end if

  call barrier
  call findgntn0_clear
  deallocate(emat12k, exclit, emat12, emat34)
  deallocate(potcl, excli)

  write(unitout, '(a)') "Info(" // trim(thisnam) // "):&
    & exchange coulomb interaction finished"

end subroutine exccoulint
!EOC

