! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_exccoulint
! !INTERFACE:
subroutine b_exccoulint(iqmt, fra, fti)
! !USES:
  use mod_misc, only: filext
  use mod_constants, only: zone, zzero
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: iqmap
  use mod_kpoint, only: nkptnr, ivknr
  use mod_lattice, only: omega
  use modinput, only: input
  use modmpi, only: rank, barrier, mpi_allgatherv_ifc
  use modxs, only: xsgnt, unitout,&
                 & ngq, nqptr,&
                 & kpari, kparf,&
                 & ppari, pparf, iqmapr,&
                 & qvkloff, bcbs
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  use m_getunit
  use modbse
  use m_b_ematqk
  use m_putgetbsemat
! !DESCRIPTION:
!   Calculates the exchange term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!   level for the treatment of core excitations(using local orbitals).
!   October 2010 (Weine Olovsson)
!   Forked from exccoulint.F90 and adapted for non TDA BSE. (Aurich)
!EOP
!BOC      

  implicit none

  ! I/O

  integer, intent(in) :: iqmt ! Index of momentum transfer Q
  logical, intent(in) :: fra  ! Construct RA coupling block
  logical, intent(in) :: fti  ! Use time inverted anti-resonant basis

  ! Local variables

  character(*), parameter :: thisnam = 'b_exccoulint'
  ! ik,jk block of V matrix (final product)
  complex(8), allocatable :: excli(:, :)
  ! Auxilliary arrays for the construction of excli
  complex(8), allocatable :: ematouk(:, :, :, :), ematuok(:, :, :, :)
  complex(8), allocatable :: mou(:, :, :), muo(:,:,:)
  ! Truncated coulomb potential
  real(8), allocatable :: potcl(:)
  ! ik jk q points 
  integer(4) :: ikkp
  integer(4) :: ik, jk, iknr, jknr, iqr, iq
  ! Number of occupied/unoccupied states at ik and jk
  integer(4) :: ino, inu, jno, jnu
  ! Number of transitions at ik and jk
  integer(4) :: inou, jnou
  ! State loop indices
  integer(4) :: io, jo, iu, ju
  ! Aux.
  integer(4) :: igq1, numgq
  integer(4) :: iv(3), j1, j2
  ! Timinig vars
  real(8) :: tpw1, tpw0

  !---------------!
  !   main part   !
  !---------------!

write(*,*) "Hello, this is b_exccoulint at rank:", mpiglobal%rank

  ! Sanity check 
  if(fra) then
    if(fti) then 
      if(mpiglobal%rank == 0) then
        write(unitout,*) "Info(b_exccoulint): V^{RA,ti} = V^{RR}, returning."
      end if
      return
    end if
  end if

  ! General setup
  call init0
  ! k-point setup
  call init1
  ! q-point and qmt-point setup
  call init2

  ! Check number of empty states
  if(input%xs%screening%nempty .lt. input%groundstate%nempty) then
    write(*,*)
    write(*, '("Error(",a,"): too few empty states in screening eigenvector file&
      & - the screening should include many empty states (bse/screening)", 2i8)')&
      & trim(thisnam), input%groundstate%nempty, input%xs%screening%nempty
    write(*,*)
    call terminate
  end if

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
    write(unitout, '(a,3i8)') 'Info(' // thisnam // '):&
      & Gaunt coefficients generated within lmax values:',&
      & input%groundstate%lmaxapw, input%xs%lmaxemat, input%groundstate%lmaxapw
    write(unitout, '(a, i6)') 'Info(' // thisnam // '):&
      & Number of reduced q-points: ', nqptr
    call flushifc(unitout)
  end if

  ! Read Fermi energy from file
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
      call genfilname(basename=trim(infofbasename)//'_'//trim(exclicfbasename),&
        & iqmt=iqmt, filnam=infofname)
      call b_putbseinfo(infofname, iqmt)
    else
      call genfilname(basename=trim(infofbasename)//'_'//trim(exclifbasename),&
        & iqmt=iqmt, filnam=infofname)
      call b_putbseinfo(infofname, iqmt)
    end if
  end if

  ! Set output file names
  if(fra) then
    call genfilname(basename=exclicfbasename, iqmt=iqmt, filnam=exclifname)
  else
    call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=exclifname)
  end if

  ! Change file extension and write out k an q points
  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(mpiglobal%rank == 0) then
    call writekpts
    call writeqpts
  end if
  ! Revert to previous file extension
  call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

  ! Number of G+qmt points 
  numgq = ngq(iqmt+1) ! ngq(1) is the qmt=0 case

  ! Calculate radial integrals used in the construction 
  ! of the plane wave matrix elements
  call ematrad(iqmt+1)

  ! Allocate \bar{v}_{G}(qmt)
  allocate(potcl(numgq))
  potcl = 0.d0

  ! Call init1 with a vkloff that is derived
  ! from the momentum transfer q vectors (qmt).
  ! Currently BSE only works for vqmt=0, so that
  ! at the moment this basically does nothing, it
  ! just calls ini1 again with the same offset as 
  ! the groundstate vkloff.
  call init1offs(qvkloff(1, iqmt+1))

  ! Allocate eigenvalue/eigenvector related
  ! quantities for use in ematqk
  call ematqalloc

  ! Work array to store result of plane wave matrix elements
  ! for each considered k point 
! Note: Potentially too large to keep in memory 
  if(fra) then
    allocate(ematuok(nu_bse_max, no_bse_max, numgq, nk_bse))
    allocate(ematouk(no_bse_max, nu_bse_max, numgq, nk_bse))
  else
    allocate(ematouk(no_bse_max, nu_bse_max, numgq, nk_bse))
  end if

  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_exccoulint): Generating plane wave matrix elements")')
    call timesec(tpw0)
  end if

  !! Plane wave matrix elements calculation.
  !! RR:  M_ouk(G,qmt) for all k.
  !! RA:  M_uok-qmt(G,qmt) for all k.
  !! Currently only works for qmt=0.
  !---------------------------!
  !     Loop over k-points    !
  !---------------------------!
  ! Parallelize over non reduced k-points
  ! participating in the BSE
  call genparidxran('k', nk_bse)

  ! MPI distributed loop
  do ik = kpari, kparf

    ! Get golbal non reduced k point index
    ! from the BSE k-point index set.
    iknr = kmap_bse_rg(ik)

    ! Get the number of participating occupied/unoccupied
    ! states at current k point
    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1

    ! Get plane wave elements for io iu combinations 
    ! for non reduced k and qmt=0 
    ! and save them for all k points for later use.
    if(fra) then 
      call getmuo(iqmt+1, iknr, muo)
      ematuok(1:inu, 1:ino, 1:numgq, ik) = muo
      call getmou(iqmt+1, iknr, mou)
      ematouk(1:ino, 1:inu, 1:numgq, ik) = mou
    else
      call getmou(iqmt+1, iknr, mou)
      ematouk(1:ino, 1:inu, 1:numgq, ik) = mou
    end if

    if(mpiglobal%rank == 0) then
      write(6, '(a,"Exccoulint - mou/muo progess:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(ik-kpari+1)/dble(kparf-kpari+1)
      flush(6)
    end if

  end do

  ! Helper no longer needed
  if(allocated(mou)) deallocate(mou)
  if(allocated(muo)) deallocate(muo)

  ! Communicate array-parts wrt. k-points
  if(fra) then
    call mpi_allgatherv_ifc(set=nk_bse, rlen=nu_bse_max*no_bse_max*numgq, zbuf=ematuok,&
      & inplace=.true., comm=mpiglobal)
    call mpi_allgatherv_ifc(set=nk_bse, rlen=no_bse_max*nu_bse_max*numgq, zbuf=ematouk,&
      & inplace=.true., comm=mpiglobal)
  else
    call mpi_allgatherv_ifc(set=nk_bse, rlen=no_bse_max*nu_bse_max*numgq, zbuf=ematouk,&
      & inplace=.true., comm=mpiglobal)
  end if

  if(mpiglobal%rank == 0) then
    write(*,*)
  end if
  if(mpiglobal%rank == 0) then
    call timesec(tpw1)
    write(unitout, '("  Timing (in seconds)	   :", f12.3)') tpw1 - tpw0
  end if

  !! Generation of V matrix.
  !! Currently only works for qmt=0.
  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!
  call genparidxran('p', nkkp_bse)

  allocate(excli(nou_bse_max, nou_bse_max))

  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_exccoulint): Generating V matrix elements")')
    call timesec(tpw0)
  end if

  kkp: do ikkp = ppari, pparf

    ! Get individual k-point indices from combined kk' index.
    !   ik runs from 1 to nk_bse, jk from ik to nk_bse
    call kkpmap(ikkp, nk_bse, ik, jk)
    !! Note: If the exchange matrix 
    !! has the elements v_{i,j} and the indices enumerate the
    !! states according to
    !! i = {u1o1k1, u2o1k1, ..., uMo1k1,
    !!      uMo2k1, ..., uMoNk1, o1u1k2, ..., oMuNkO} -> {1,...,M*N*O}
    !! then because of v_{j,i} = v^*_{i,j} only kj = ki,..,kN is 
    !! needed.

    ! Get total k point indices
    iknr = kmap_bse_rg(ik)
    jknr = kmap_bse_rg(jk) 

    ! Get number of transitions at ik,jk
    inou = kousize(iknr)
    jnou = kousize(jknr)

    !!<-- The difference vector, i.e. q not needed in
    !! the calculations (yet?)
    ! K-point difference k_j-k_i on integer grid.
    iv(:) = ivknr(:, jknr) - ivknr(:, iknr)
    ! Map to reciprocal unit cell 01
    iv(:) = modulo(iv(:), input%groundstate%ngridk(:))
    ! Q-point(reduced)
    iqr = iqmapr(iv(1), iv(2), iv(3))
    ! Q-point(non-reduced)
    iq = iqmap(iv(1), iv(2), iv(3))
    !!-->

    !! THIS IS NOT (YET) DEPENDENT ON THE LOOP!
    ! Set G=0 term of coulomb potential to zero [Ambegaokar-Kohn]
    potcl(1) = 0.d0
    ! Set up coulomb potential
    ! For G/=0 construct it via v^{1/2}(G,q)*v^{1/2}(G,q),
    ! which corresponds to the first passed flag=0.
    ! Here iqmt=1 is fixed, which correspons to qmt=0.
    do igq1 = 2, numgq
      call genwiqggp(0, iqmt+1, igq1, igq1, potcl(igq1))
    end do

    call makeexcli(excli(1:inou,1:jnou))

    ! Parallel write
    if(fra) then
      call b_putbsemat(exclifname, 78, ikkp, iqmt, excli)
    else
      call b_putbsemat(exclifname, 77, ikkp, iqmt, excli)
    end if

    if(mpiglobal%rank == 0) then
      write(6, '(a,"Exccoulint progess:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(ikkp-ppari+1)/dble(pparf-ppari+1)
      flush(6)
    end if

  ! End loop over(k,kp) pairs
  end do kkp

  if(mpiglobal%rank == 0) then
    write(*,*)
  end if
  
  call barrier

  if(mpiglobal%rank == 0) then
    call timesec(tpw1)
    write(unitout, '("  Timing (in seconds)	   :", f12.3)') tpw1 - tpw0
  end if

  call findgntn0_clear

  deallocate(potcl) 
  if(allocated(ematouk)) deallocate(ematouk)
  if(allocated(ematuok)) deallocate(ematuok)
  deallocate(excli)

  if(mpiglobal%rank .eq. 0) then
    write(unitout, '("Info(b_exccoulint): Exchange coulomb interaction&
      & finished for iqmt=",i8)') iqmt
  end if

  contains 

    subroutine getmou(iqnr, iknr, mou)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: mou(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the io iu plane wave elements
      ! Pack selected o-u combinations in struct
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=inu
      ematbc%il2=koulims(1,iknr)
      ematbc%iu2=koulims(2,iknr)

      ! Allocate space for M_{io iu,G} at fixed (k, q)
      if(allocated(mou)) deallocate(mou)
      allocate(mou(ino,inu,numgq))

      ! Calculate M_{io iu,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, mou, ematbc)

    end subroutine getmou

    subroutine getmuo(iqnr, iknr, muo)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: muo(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the iu io plane wave elements
      ! Pack selected uo combinations in struct
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=ino
      ematbc%il2=koulims(3,iknr)
      ematbc%iu2=koulims(4,iknr)

      !! As long as q=0 this reallocating should not be needed
      ! Allocate space for M_{iu io,G} at fixed (k, q)
      if(allocated(muo)) deallocate(muo)
      allocate(muo(inu,ino,numgq))

      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, muo, ematbc)

    end subroutine getmuo

    subroutine makeexcli(excli)

      complex(8), intent(out) :: excli(inou, jnou) 

      ! Work arrays
      complex(8) :: emat12(inou, numgq), emat34(jnou, numgq)
      integer(4) :: iaoff, jaoff, ia, ja

      ! Offset in combined index for ik and jk
      iaoff = sum(kousize(1:iknr-1))
      jaoff = sum(kousize(1:jknr-1))

      !! RR and RA part
      do ia = 1, inou
        ! Get iu index relative to iu range of ik
        iu = smap_rel(1, ia+iaoff) ! iu (ik+qmt)
        ! Get io index relative to io range of ik
        io = smap_rel(2, ia+iaoff) ! io (ik)
        ! emat12_ia = M_o1u1ki, M_o2u1ki, ..., M_oNu1ki, M_o1u2ki, ..., M_oNuMki
        emat12(ia, :) = ematouk(io, iu, :, ik)
      end do
      ! M_ou -> M^*_ou
      emat12 = conjg(emat12)

      if(.not. fra) then

        !! RR part
        do ja = 1, jnou
          ! Get ju index relative to iu range of jk
          ju = smap_rel(1, ja+jaoff) ! ju (jk+qmt)
          ! Get jo index relative to io range of jk
          jo = smap_rel(2, ja+jaoff) ! jo (jk)
          ! emat34_ja = (M_o1u1kj,M_o2u1kj,...,M_oNu1kj,M_o1u2kj,...,M_oNuMkj)*\bar{v}
          emat34(ja, :) = ematouk(jo, ju, :, jk) * potcl(:)
        end do
        ! Calculate exchange matrix elements: 
        ! excli_{ia, ja} = \Sum_{G} emat12_{ia,G} (emat34^T)_{G,ja}
        !   i.e. excli_{iu io ik, ju jo jk}(qmt) =
        !          \Sum_{G} M^*_{io iu ik}(G,qmt) M_{jo ju jk}(G,qmt) v(G,qmt)
        call zgemm('n', 't', inou, jnou, numgq, zone/omega/nk_bse,&
          & emat12, inou, emat34, jnou, zzero, excli, inou)

      else

        !! RA part
        do ja = 1, jnou
          ! Get ju index relative to iu range of jk
          ju = smap_rel(1, ja+jaoff) ! ju (jk-qmt)
          ! Get jo index relative to io range of jk
          jo = smap_rel(2, ja+jaoff) ! jo (jk)
          ! emat34_ja = (M_u1o1kj-qmt,M_u2o1kj-qmt,...,
          !               M_uMo1kj-qmt,M_u1o2kj-qmt,...,M_uMoNkj-qmt)*v(G,0)
          emat34(ja, :) = ematuok(ju, jo, :, jk) * potcl(:)
        end do

        ! Calculate coupling exchange matrix elements: v_{1234} = m_{12}^* m_{34}^t
        ! excli_{j1, j2} = \Sum_{G} emat12_{j1,G} (emat34^T)_{G,j2}
        !   i.e. excli_{iu io ik, ju jo jk}(qmt) =
        !          \Sum_{G} M^*_{io iu ik}(G, qmt) M_{ju jo jk-qmt}(G, qmt) v(G, qmt)
        call zgemm('n', 't', inou, jnou, numgq, zone/omega/nk_bse,&
          & emat12, inou, emat34, jnou, zzero, excli, inou)

      end if

    end subroutine makeexcli

end subroutine b_exccoulint
!EOC

