! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_exccoulint
! !INTERFACE:
subroutine b_exccoulint
! !USES:
  use mod_misc, only: filext
  use mod_constants, only: zone, zzero
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: nqpt, iqmap
  use mod_kpoint, only: nkptnr, ivknr
  use mod_lattice, only: omega
  use modinput, only: input
  use modmpi, only: rank, barrier, mpi_allgatherv_ifc
  use modxs, only: xsgnt, unitout,&
                 & ngq,&
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

  ! Local variables
  character(*), parameter :: thisnam = 'exccoulint'
  integer(4) :: iqmt

  logical :: fcoup
  integer(4) :: iknr, jknr, iqr, iq, igq1, numgq
  integer(4) :: iv(3), j1, j2
  integer(4) :: ikkp

  integer(4) :: ik, jk, inou, jnou, ino, inu, jno, jnu
  integer(4) :: io, jo, iu, ju

  real(8), allocatable :: potcl(:)
  complex(8), allocatable :: excli(:, :, :, :), exclic(:, :, :, :)
  complex(8), allocatable :: ematouk(:, :, :, :), ematuok(:, :, :, :)
  complex(8), allocatable :: mou(:, :, :), muo(:,:,:)

  !---------------!
  !   main part   !
  !---------------!

write(*,*) "Hello, this is b_exccoulint at rank:", rank

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

  ! Read Fermi energy from file
  ! Which one .OUT _QMTXXX.OUT or _SCR.OUT ??
  ! When called during whole BSE chain, this is set to _SCR.OUT ?
  ! When calling it using doonly is set to .OUT ?
  write(*,'("b_exccoulint@rank",i3," : Reading fermi energy form file: ",a)')&
    & mpiglobal%rank, 'EFERMI'//trim(filext)
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
    write(unitout, '(a,3i8)') 'Info(' // thisnam // '):&
      & Gaunt coefficients generated within lmax values:',&
      & input%groundstate%lmaxapw, input%xs%lmaxemat, input%groundstate%lmaxapw
    write(unitout, '(a, i6)') 'Info(' // thisnam // '):&
      & Number of q-points: ', nqpt
    call flushifc(unitout)
  end if

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
  if(read_eval_occ_qmt) then
    call setranges_modxs(iqmt)
  end if

  ! Select relevant transitions for the construction
  ! of the BSE hamiltonian
  ! Also sets nkkp_bse, nk_bse 
  if(seltrans) then
    call select_transitions(iqmt)
  end if

  ! Include coupling terms
  fcoup = input%xs%bse%coupling

  ! Set output file names
  call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=exclifname)
  if(fcoup) then
    call genfilname(basename=exclicfbasename, iqmt=iqmt, filnam=exclicfname)
  end if

  ! Change file extension and write out k an q points
  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  if(rank .eq. 0) then
    call writekpts
    call writeqpts
  end if
  ! Revert to previous file extension
  call genfilname(iqmt=max(0, iqmt), setfilext=.true.)

  ! Number of G+q points (q=0)
  numgq = ngq(iqmt+1) ! ngq(1) is the q=0 case

  ! Calculate radial integrals used in the construction 
  ! of the plane wave matrix elements (for q=0)
  call ematrad(iqmt+1)

  allocate(potcl(numgq))
  potcl = 0.d0

  ! Call init1 with a vkloff that is derived
  ! from the momentum transfer q vectors.
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
  allocate(ematouk(no_bse_max, nu_bse_max, numgq, nk_bse))
  if(fcoup == .true.) then
    allocate(ematuok(nu_bse_max, no_bse_max, numgq, nk_bse))
  end if

  !!<-- Generate M_ouk(G,qmt) for all k.
  !! If not TDA, then also generate M_uok-qmt(G,qmt) needs to be created.
  !! Currently only works for qmt=0.
  !---------------------------!
  !     Loop over k-points    !
  !---------------------------!
  ! Parallelize over non reduced k-points
  ! participating in the BSE
  call genparidxran('k', nk_bse)

  ! MPI distributed loop
  do ik = kpari, kparf

    iknr = kmap_bse_rg(ik)

    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1

    ! Get plane wave elements for io iu combinations 
    ! for non reduced k and q=0 
    call getmou(iqmt+1, iknr, mou)

    ! Save them for all ks
    ematouk(1:ino,1:inu, :, ik) = mou

    if(fcoup == .true.) then
      call getmuo(iqmt+1, iknr, muo)
      ematuok(1:inu,1:ino,:,iknr) = muo
    end if

  end do

  ! Communicate array-parts wrt. k-points
  call mpi_allgatherv_ifc(set=nk_bse, rlen=nou_bse_max*numgq, zbuf=ematouk,&
    & inplace=.true., comm=mpiglobal)
  if(fcoup == .true.) then
    call mpi_allgatherv_ifc(set=nk_bse, rlen=nou_bse_max*numgq, zbuf=ematuok,&
      & inplace=.true., comm=mpiglobal)
  end if
  !!-->

  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!
  call genparidxran('p', nkkp_bse)

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

    ! Get ranges
    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1
    inou = inu*ino
    jnu = koulims(2,jknr) - koulims(1,jknr) + 1
    jno = koulims(4,jknr) - koulims(3,jknr) + 1
    jnou = jnu*jno

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
    ! Here iqmt=1 is fixed, which is q=0.
    do igq1 = 2, numgq
      call genwiqggp(0, iqmt+1, igq1, igq1, potcl(igq1))
    end do

    call makeexcli()

    ! Parallel write
    call b_putbsemat(exclifname, 77, ikkp, iqmt, excli, nou_bse_max**2)
    if(fcoup == .true.) then
      call b_putbsemat(exclicfname, 78, ikkp, iqmt, exclic, nou_bse_max**2)
    end if

  ! End loop over(k,kp) pairs
  end do kkp
  
  call barrier

  call findgntn0_clear

  deallocate(potcl) 
  deallocate(ematouk)
  deallocate(excli)
  deallocate(mou)
  if(fcoup == .true.) then
    deallocate(ematuok)
    deallocate(exclic)
    deallocate(muo)
  end if

  if(rank .eq. 0) then
    write(unitout, '(a)') "Info(" // trim(thisnam) // "):&
      & exchange coulomb interaction finished"
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

    subroutine makeexcli

      ! Work arrays
      complex(8) :: exclit(inou, jnou) 
      complex(8) :: emat12(inou, numgq), emat34(jnou, numgq)

      ! Result arryas
      if(allocated(excli)) deallocate(excli)
      allocate(excli(inu, ino, jnu, jno))
      if(fcoup == .true.) then
        if(allocated(exclic)) deallocate(exclic)
        allocate(exclic(inu, ino, jnu, jno))
      end if

      !!<-- RR and RA part
     ! do iu = 1, inu   ! iu (ik+qmt)
     !   do io = 1, ino ! io (ik)
     !     j1 = io + (iu-1)*inu ! ioiu
     !     ! emat12_j = M_o1u1ki, M_o2u1ki, ..., M_oNu1ki, M_o1u2ki, ..., M_oNuMki
     !     emat12(j1, :) = ematouk(io1, iu1, :, ik)
     !   end do
     ! end do
      ! Does the same commented code above 
      emat12 = reshape(ematouk(1:ino,1:inu,:,ik),[inou, numgq])
      ! M_ou -> M^*_ou
      emat12 = conjg(emat12)
      !!-->

      !!<-- RR part
      do ju =1, jnu    ! ju (jk+qmt)
        do jo = 1, jno ! jo (jk)
          j2 = jo + (ju - 1) * jnu ! joju
          ! emat34_j = (M_o1u1kj,M_o2u1kj,...,M_oNu1kj,M_o1u2kj,...,M_oNuMkj)*v(G,0)
          emat34(j2, :) = ematouk(jo, ju, :, jk) * potcl(:)
        end do
      end do
      ! Calculate exchange matrix elements: v_{1234} = m_{12}^* m_{34}^t
      ! exclit_{j1, j2} = \Sum_{G} emat12_{j1,G} (emat34^T)_{G,j2}
      !   i.e. exclit_{io iu ik, jo ju jk}(qmt) =
      !          \Sum_{G} M^*_{io iu ik}(G,qmt) M_{jo ju jk}(G,qmt) v(G,qmt)
      call zgemm('n', 't', inou, jnou, numgq, zone/omega/nk_bse,&
        & emat12, inou, emat34, jnou, zzero, exclit, inou)

      ! Map back to individual band indices
      do ju = 1, jnu   ! ju
        do jo = 1, jno ! jo
          j2 = jo + (ju-1)*jnu ! joju
          do iu = 1, inu   ! iu
            do io = 1, ino ! io
              j1 = io + (iu-1)*inu ! ioiu
              ! exclit_{io_j1 iu_j1, jo_j2 ju_j2} -> excli_{iu_j1 io_j1, ju_j2 jo_j2}
              excli(iu, io, ju, jo) = exclit(j1, j2)
            end do
          end do
        end do
      end do
      !!-->

      if(fcoup == .true.) then
        !!<-- RA part
        do jo = 1, jno   ! jo (jk)
          do ju =1, jnu  ! ju (jk-qmt)
            j2 = ju + (jo-1)*jno
            ! emat34_j = (M_u1o1kj-qmt,M_u2o1kj-qmt,...,
            !               M_uMo1kj-qmt,M_u1o2kj-qmt,...,M_uMoNkj-qmt)*v(G,0)
            emat34(j2, :) = ematuok(ju, jo, :, jk) * potcl(:)
          end do
        end do
        ! Calculate coupling exchange matrix elements: v_{1234} = m_{12}^* m_{34}^t
        ! exclitc_{j1, j2} = \Sum_{G} emat12_{j1,G} (emat34^T)_{G,j2}
        !   i.e. exclitc_{io iu ik, ju jo jk}(qmt) =
        !          \Sum_{G} M^*_{io iu ik}(G, qmt) M_{ju jo jk-qmt}(G, qmt) v(G, qmt)
        call zgemm('n', 't', inou, jnou, numgq, zone/omega/nk_bse,&
          & emat12, inou, emat34, jnou, zzero, exclit, inou)
        ! Map back to individual band indices
        do jo = 1, jno   ! jo
          do ju = 1, jnu ! ju
            j2 = ju + (jo-1)*jno !jujo
            do iu = 1, inu   ! iu
              do io = 1, ino ! io
                j1 = io + (iu-1)*inu !ioiu
                ! exclit_{io_j1 iu_j1, ju_j2 jo_j2} -> exclic_{iu_j1 io_j1, ju_j2 jo_j2}
                exclic(iu, io, ju, jo) = exclit(j1, j2)
              end do
            end do
          end do
        end do
        !!-->
      end if

    end subroutine makeexcli

end subroutine b_exccoulint
!EOC

