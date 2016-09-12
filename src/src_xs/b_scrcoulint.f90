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
  use mod_misc, only: task
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
                 & vqlr, vgqc, dielten,&
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
  complex(8) :: zt1,prefactor
  complex(8), allocatable :: scclit(:, :), sccli(:, :, :, :), scclid(:, :)
  complex(8), allocatable :: scieffg(:, :, :), tm(:, :), tmi(:, :), bsedt(:, :),zm(:,:)
  complex(8), allocatable :: phf(:, :), cmoo(:, :), cmuu(:, :)
  real(8), parameter :: epsortho = 1.d-12
  real(8) :: vqr(3), vq(3), t1, ta
  integer :: ikkp, iknr, jknr, iqr, iq, iqrnr, jsym, jsymi, igq1, n, reclen
  integer :: nsc, iv(3), ivgsym(3), j1, j2, nkkp
  integer(4) :: nou, noo, nuu, no, nu
  integer(4) :: io1, io2, iu1, iu2
  integer :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer, allocatable :: igqmap(:)
  logical :: tq0, tphf

  ! Plane wave arrays
  type(bcbs) :: ematbc
  complex(8), allocatable :: muu(:, :, :)
  complex(8), allocatable :: moo(:, :, :)

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

  write(unitout, '(a,3i8)') 'info(' // thisnam // '):&
    & Gaunt coefficients generated within lmax values:', input%groundstate%lmaxapw,&
    & input%xs%lmaxemat, input%groundstate%lmaxapw

  write(unitout, '(a, i6)') 'info(' // thisnam // '): number of q-points: ', nqpt

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


  call flushifc(unitout)

  ! Set EVALSV_SCR.OUT as basis for the occupation limits search
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Set ist* variables in modxs using findocclims
  call setranges_modxs(0)

  ! Set band combinations (modbse:bcou & modbse:bcouabs)
  call setbcbs_bse

  ! Number of occupied states
  no = bcou%n1
  ! Number of unoccupied states
  nu = bcou%n2

  ! Number of oo-combinations
  noo = no * no
  ! Number of ou-combinations
  nou = no * nu

  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if

  ! Allocate local arrays for screened coulomb interaction
  ! Phased for transformations form reduced q points to non reduced ones.
  allocate(phf(ngqmax, ngqmax))
  ! W arrays for one k2-k1=q
  allocate(sccli(no, nu, no, nu), scclid(no, nu))
  ! W(G,G',q)
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

  do iqr = qpari, qparf ! Reduced q

    call chkpt(3, (/ task, 1, iqr /),&
      & 'task,sub,reduced q-point; generate effective screened coulomb potential')

    ! Locate reduced q-point in non-reduced set
    iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))

    ! Get number of G+q vectors for current q
    n = ngq(iqrnr)

    ! Calculate effective screened coulomb interaction
    ! by inverting the symmetrized RPA dielectric matrix for a given q and
    ! 0 frequency and then multiplying
    ! it with v^{1/2} from both sides.
    call genscclieff(iqr, ngqmax, n, scieffg(1,1,iqr))

    ! Generate radial integrals for matrix elements of plane wave
    call putematrad(iqr, iqrnr)

  end do

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(nqptr,ngqmax*ngqmax,zbuf=scieffg)
  call barrier


  ! Combinations of k and k', where ik'>=ik
  nkkp = (nkptnr*(nkptnr+1)) / 2

  ! Information on size of output file
  inquire(iolength=reclen) ikkp, iknr, jknr, iq, iqr, no, nu, no, nu, sccli
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
      & 'task, sub,(k,kp)-pair; direct term of BSE hamiltonian')

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
    n = ngq(iq)

    allocate(igqmap(n), cmoo(noo, n), cmuu(nuu, n))
    allocate(tm(n, n), tmi(n, n))
    allocate(scclit(noo, nuu))

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
    ! Get radial integrals (previously calculated for reduced q set)
    call getematrad(iqr, iq)
    ! Rotate radial integrals
    ! (i.e. apply symmetry transformation to the result for non reduced q point)
    call rotematrad(n, igqmap)
    ! Rotate inverse of screening, coulomb potential and radial integrals
    ! (i.e. apply symmetry transformation to the result for non reduced q point)
    tmi(:,:) = phf(:n, :n) * scieffg(igqmap, igqmap, iqr)
    !!-->

    ! Allocate arrays used in ematqk
    ! (those independent of the band selection) 
    call ematqalloc

    !! Calculate the oo plane wave elements
    ! Pack selected o-o combinations in struct
    ematbc%n1=no
    ematbc%il1=bcouabs%il1
    ematbc%iu1=bcouabs%iu1
    ematbc%n2=no
    ematbc%il2=bcouabs%il1
    ematbc%iu2=bcouabs%iu1
    ! Allocate space for M_{o1o2,G} at fixed (k, q)
    if(allocated(moo)) deallocate(moo)
    allocate(moo(no,no,n))
    ! Calculate M_{o1o2,G} at fixed (k, q)
    call ematqk(iq, iknr, moo, ematbc)

    !! Calculate the uu plane wave elements
    ! Pack selected o-o combinations in struct
    ematbc%n1=nu
    ematbc%il1=bcouabs%il2
    ematbc%iu1=bcouabs%iu2
    ematbc%n2=nu
    ematbc%il2=bcouabs%il2
    ematbc%iu2=bcouabs%iu2
    ! Allocate space for M_{u1u2,G} at fixed (k, q)
    if(allocated(muu)) deallocate(muu)
    allocate(muu(nu,nu,n))
    ! Calculate M_{o1o2,G} at fixed (k, q)
    call ematqk(iq, iknr, muu, ematbc)

    call chkpt(3, (/ task, 2, ikkp /),&
      & 'task,sub,(k,kp)-pair; direct term of BSE hamiltonian')

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
        ! Full screening tm is set to
        ! W_{GG'}(q,\omega=0)
        tm(:, :) = tmi(:, :)
    end select

    ! Combine indices for matrix elements of plane wave
    j1 = 0
    ! Occupied
    do io1 = 1, no
      ! Occupied
      do io2 = 1, no
        j1 = j1 + 1
        ! cmoo_j = M_o1o1, M_o2o1, ..., M_oNo1, M_o1o2, ..., M_oNoN
        cmoo(j1, :) = moo(io2, io1, :)
      end do
    end do
    j2 = 0
    ! Unoccupied
    do iu1 = 1, nu
      ! Unoccupied
      do iu2 = 1, nu
        j2 = j2 + 1
        ! cmuu_j = M_u1u1, M_u2u1, ..., M_uNu1, M_u1u2, ..., M_uNuN
        cmuu(j2, :) = muu(iu2, iu1, :)
      end do
    end do

    ! Matrix elements of direct term (as in BSE-code of Peter and
    ! In the self-documentation of Andrea Marini)

    prefactor=1.0d0/(omega*dble(nkptnr))

    ! M_oioj -> M^*_oioj
    cmoo(:,:)=conjg(cmoo(:,:))

    ! Allocate helper array of dimension (#o*#o,#G) (same as cmoo)
    allocate(zm(noo,n))

    ! Calculate matrix elements of screened coulomb interaction scclit_{o_j1 o'_j1, u_j2 u'_j2}(q)
    ! zm = cmoo * tm
    !   i.e. zm_{j,G} = \Sum_{G'} cmoo_{j,G'} tm_{G',G}
    !   i.e. zm_{o_j,o'_j}(G,q) = \Sum_{G'} M^*_{o_j,o'_j}(G',q) W(G',G, q)
    call zgemm('n', 'n', noo, n, n, zone, cmoo, noo, tm, n, zzero, zm, noo)
    ! scclit = prefactor * zm * cmuu^T
    !   i.e. scclit(j1, j2) = \Sum_{G} zm_{j1,G} (cmuu^T)_{G,j2}
    !   i.e. scclit_{o_j1 o'_j2, u_j2 u'_j2} = \Sum{G,G'} M^*_{o_j1,o'_j1}(G,q) W(G,G', q) M_{u_j2 u'_j2}(G',q)
    call zgemm('n', 't', noo, nuu, n, prefactor, zm,&
      & noo, cmuu, nuu, zzero, scclit, noo)
    deallocate(zm)        

    ! Map back to individual band indices
    j2 = 0
    ! Unoccupied
    do iu1 = 1, nu
      ! Unoccupied
      do iu2 = 1, nu
        j2 = j2 + 1
        j1 = 0
        ! Occupied
        do io1 = 1, no
          ! Occupied
          do io2 = 1, no
            j1 = j1 + 1
            ! scclit_{o_j1 o'_j1, u_j2 u'_j2} -> sccli_{o_j1 u_j2, o'_j1 u'_j2}
            sccli(io2, iu2, io1, iu1) = scclit(j1, j2)
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
          t1 = dble(zt1)
          bsedt(1, rank) = min(dble(bsedt(1, rank)), t1)
          bsedt(2, rank) = max(dble(bsedt(2, rank)), t1)
          bsedt(3, rank) = bsedt(3, rank) + zt1 / (no*nu)
        end do
      end do
    end if

    ! Parallel write
    call putbsemat('SCCLI.OUT', sccli, ikkp, iknr, jknr,&
      & iq, iqr, no, nu, no, nu)

    deallocate(moo, muu)
    deallocate(igqmap, cmoo, cmuu)
    deallocate(tm, tmi)
    deallocate(scclit)

  ! End loop over(k,kp)-pairs
  end do kkploop

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

  write(unitout, '("Info(scrcoulint): Screened coulomb interaction&
    & finished")')

end subroutine b_scrcoulint
!EOC

