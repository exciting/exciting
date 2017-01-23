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
  use modmpi
  use modxs, only: sta1, sto1, sta2, sto2,&
                 & xsgnt, unitout, istocc0, istocc,&
                 & istunocc0, istunocc, isto0, isto,&
                 & istu0, istu, ksgap, ngq,&
                 & kpari, kparf, xiou, xiuo,&
                 & ppari, pparf, iqmapr, nst1,&
                 & nst2, qvkloff, ikmapikq
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
  integer, parameter :: iqmt = 1
  real(8), parameter :: epsortho = 1.d-12
  integer :: iknr, jknr, iqr, iq, igq1, n
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
  ! emattype=1 corresponds to 12=ou,34=uo combinations
  input%xs%emattype = 1

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
  rnst4 = sto1-sta1+1
  ! Number of unoccupied states
  rnst2 = sto2-sta2+1
  rnst3 = sto2-sta2+1

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
    & Gaunt coefficients generated within lmax values:',&
    & input%groundstate%lmaxapw, input%xs%lmaxemat, input%groundstate%lmaxapw

  write(unitout, '(a, i6)') 'Info(' // thisnam // '):&
    & Number of q-points: ', nqpt

  call flushifc(unitout)
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Find occupation limits for k and k+q, where q=0
  call findocclims(0, ikmapikq(:,1), istocc0, istunocc0, isto0, isto, istu0, istu)
  istunocc = istunocc0
  istocc = istocc0
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

  ! Set nstX, istlX, istuX variables for X: 1=o, 2=u, 3=o, 4=u
  call ematbdcmbs(input%xs%emattype)

  ! Number of ou-combinations
  nst12 = rnst1 * rnst2
  ! Number of uo-combinations
  nst34 = rnst3 * rnst4
  ! Number of ou-combinations
  nst13 = rnst1 * rnst3
  ! Number of uo-combinations
  nst24 = rnst2 * rnst4

  call genfilname(dotext='_SCI.OUT', setfilext=.true.)

  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if

  ! Number of G+q points (q=0)
  n = ngq(iqmt)

  ! Calculate radial integrals used in the construction 
  ! of the plane wave matrix elements (for q=0)
  call ematrad(iqmt)

  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  allocate(potcl(n))
  allocate(excli(rnst1, rnst2, rnst1, rnst2))
  allocate(exclit(nst12, nst34))
  allocate(emat12k(nst1, nst2, n, nkptnr))
  allocate(emat12(nst12, n), emat34(nst34, n))

  potcl(:) = 0.d0
  excli(:, :, :, :) = zzero

  !!<-- Generate M_ouk(G,0)
  !---------------------------!
  !     Loop over k-points    !
  !---------------------------!
  ! Parallelize over non reduced k-points
  call genparidxran('k', nkptnr)

  ! Call init1 with q as vkloff, so that
  ! k mesh goes over into k+q mesh.
  call init1offs(qvkloff(1, iqmt))

  ! Allocate eigenvalue/eigenvector related
  ! quantities for use in ematqk
  call ematqalloc

  do iknr = kpari, kparf
    call chkpt(3, (/ task, 1, iknr /),&
      & 'task,sub,k-point; matrix elements of plane wave')

    ! Matrix elements for k and q=0 (xiou and xiou)
    call ematqk1(iqmt, iknr)
    emat12k(:, :, :, iknr) = xiou(:, :, :)
    deallocate(xiou, xiuo)

  end do

  ! Communicate array-parts wrt. k-points
  call mpi_allgatherv_ifc(nkptnr,rlen=nst1*nst2*n,zbuf=emat12k)
  !!-->

!! emattype and bdcmbs should not have changed?
  ! Select 12=ou 34=uo
  input%xs%emattype = 1
  call ematbdcmbs(input%xs%emattype)

  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!
  nkkp = (nkptnr*(nkptnr+1)) / 2
  call genparidxran('p', nkkp)

  kkp: do ikkp = ppari, pparf

    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task,sub,(k,kp)-pair; exchange term of BSE Hamiltonian')

    ! Get individual k-point indices from combined kk' index.
    !   iknr runs from 1 to nkptnr, jknr from iknr to nkptnr
    call kkpmap(ikkp, nkptnr, iknr, jknr)
    !! Note: If the exchange matrix 
    !! has the elements v_{i,j} and the indices enumerate the
    !! states according to
    !! i = o1u1k1, o1u1k2, ..., o1u1kN, o2u1k1, ..., o2u1kN, ...,
    !!     oMu1kN, oMu2k1, ..., oMuMkN
    !! then because of v_{j,i} = v^*_{i,j} only kj = ki,..,kN is 
    !! needed.

    !!<-- The difference vector, i.e. q not needed in
    !! the calculations (yet?)
    ! K-point difference k_j-k_i on integer grid.
    iv(:) = ivknr(:, jknr) - ivknr(:, iknr)
    ! Map to reciprocal unit cell 
    iv(:) = modulo(iv(:), input%groundstate%ngridk(:))
    ! Q-point(reduced)
    iqr = iqmapr(iv(1), iv(2), iv(3))
    ! Q-point(non-reduced)
    iq = iqmap(iv(1), iv(2), iv(3))
    !!-->

    ! Set G=0 term of coulomb potential to zero [Ambegaokar-Kohn]
    potcl(1) = 0.d0

    ! Set up coulomb potential
    ! For G/=0 construct it via v^{1/2}(G,q)*v^{1/2}(G,q), where
    ! here iqmt=1 is fixed, which is q=0
    do igq1 = 2, n
      call genwiqggp(0, iqmt, igq1, igq1, potcl(igq1))
    end do

    call genfilname(dotext='_SCR.OUT', setfilext=.true.)

!! one could use getpemat as in dfq.F90
    j1 = 0
    ! Unoccupied
    do ist2 = sta2, sto2
      ! Occupied
      do ist1 = sta1, sto1
        j1 = j1 + 1
        ! emat12_j = M_o1u1ki, M_o2u1ki, ..., M_oNu1ki, M_o1u2ki, ..., M_oNuNki
        emat12(j1, :) = emat12k(ist1-sta1+1, ist2-sta2+1, :, iknr)
      end do
    end do
    j2 = 0
    ! Unoccupied
    do ist4 = sta2, sto2
      ! Occupied
      do ist3 = sta1, sto1
        j2 = j2 + 1
        ! emat34_j = (M_o1u1kj, M_o2u1kj, ..., M_oNu1kj, M_o1u2kj, ..., M_oNuNkj)*v(G,0)
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
    call putbsemat('EXCLI.OUT', 79, excli, ikkp, iknr, jknr,&
      & iq, iqr, rnst1, rnst2, rnst4, rnst3)
    call genfilname(dotext='_SCI.OUT', setfilext=.true.)

  ! End loop over(k,kp) pairs
  end do kkp

  call barrier
  call findgntn0_clear
  deallocate(emat12k, exclit, emat12, emat34)
  deallocate(potcl, excli)

  write(unitout, '(a)') "Info(" // trim(thisnam) // "):&
    & exchange coulomb interaction finished"

end subroutine exccoulint
!EOC

