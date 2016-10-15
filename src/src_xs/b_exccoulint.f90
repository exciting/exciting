! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_exccoulint
! !INTERFACE:
subroutine b_exccoulint
! !USES:
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
                 & qvkloff
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  use m_getunit
  use modbse
  use m_b_ematqk
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
  integer, parameter :: iqmt = 1

  logical :: fcoup
  integer :: iknr, jknr, iqr, iq, igq1, numgq
  integer :: iv(3), j1, j2
  integer :: ikkp

  integer(4) :: io1, io2, iu1, iu2

  real(8), allocatable :: potcl(:)
  complex(8), allocatable :: emat12(:, :), emat34(:, :)
  complex(8), allocatable :: excli(:, :, :, :), exclic(:, :, :, :)
  complex(8), allocatable :: ematouk(:, :, :, :), ematuok(:, :, :, :)
  complex(8), allocatable :: mou(:, :, :), muo(:,:,:)

  !---------------!
  !   main part   !
  !---------------!

write(*,*) "Hello, this is b_exccoulint at rank:", rank

  ! General setup
  call init0
  ! K-point setup
  call init1
  ! Q-point setup
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

  ! The extension is not used anymore...?
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Number of G+q points (q=0)
  numgq = ngq(iqmt)

  ! Calculate radial integrals used in the construction 
  ! of the plane wave matrix elements (for q=0)
  call ematrad(iqmt)

  allocate(potcl(numgq))
  potcl = 0.d0

  allocate(excli(nu, no, nu, no))
  excli = zzero
  allocate(ematouk(no, nu, numgq, nkptnr))
  ematouk = zzero

  ! Include coupling terms
  fcoup = input%xs%bse%coupling
  if(fcoup == .true.) then
    allocate(exclic(nu, no, nu, no))
    exclic=zzero
    allocate(ematuok(nu, no, numgq, nkptnr))
    ematuok=zzero
  end if

  allocate(emat12(nou, numgq), emat34(nou, numgq))
  emat12=zzero
  emat34=zzero

  !!<-- Generate M_ouk(G,0) for all k
  !!<-- If not TDA, then also generate M_uok(G,0) for all k
  !---------------------------!
  !     Loop over k-points    !
  !---------------------------!
  ! Parallelize over non reduced k-points
  call genparidxran('k', nkptnr)

  ! Call init1 with a vkloff that is derived
  ! from the momentum transfer q vectors.
  ! Currently BSE only works for vqmt=0, so that
  ! at the moment this basically dones nothing, it
  ! just calls ini1 again with the same offset as 
  ! the groundstate vkloff.
  call init1offs(qvkloff(1, iqmt))

  ! Allocate eigenvalue/eigenvector related
  ! quantities for use in ematqk
  call ematqalloc

  ! MPI distributed loop
  do iknr = kpari, kparf
    ! Get plane wave elements for ou combinations 
    ! for non reduced k and q=0 
    call getmou(iqmt, iknr, mou)
    ! Save them for all ks
    ematouk(:, :, :, iknr) = mou

    if(fcoup == .true.) then
      call getmuo(iqmt, iknr, muo)
      ematuok(:,:,:,iknr) = muo
    end if

  end do

  ! Communicate array-parts wrt. k-points
  call mpi_allgatherv_ifc(nkptnr,no*nu*numgq,zbuf=ematouk)
  if(fcoup == .true.) then
    call mpi_allgatherv_ifc(nkptnr,no*nu*numgq,zbuf=ematuok)
  end if
  !!-->

  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!
  nkkp = (nkptnr*(nkptnr+1)) / 2
  call genparidxran('p', nkkp)

  kkp: do ikkp = ppari, pparf

    ! Get individual k-point indices from combined kk' index.
    !   iknr runs from 1 to nkptnr, jknr from iknr to nkptnr
    call kkpmap(ikkp, nkptnr, iknr, jknr)

    !! Note: If the exchange matrix 
    !! has the elements v_{i,j} and the indices enumerate the
    !! states according to
    !! i = {u1o1k1, u2o1k1, ..., uMo1k1,
    !!      uMo2k1, ..., uMoNk1, o1u1k2, ..., oMuNkO} -> {1,...,M*N*O}
    !! then because of v_{j,i} = v^*_{i,j} only kj = ki,..,kN is 
    !! needed.

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
      call genwiqggp(0, iqmt, igq1, igq1, potcl(igq1))
    end do

    call makeexcli()

    ! Parallel write
    call putbsemat('EXCLI.OUT', 77, excli, ikkp, iknr, jknr,&
      & iq, iqr, nu, no, nu, no)

    if(fcoup == .true.) then
      call putbsemat('EXCLIC.OUT', 78, exclic, ikkp, iknr, jknr,&
        & iq, iqr, nu, no, nu, no)
    end if

  ! End loop over(k,kp) pairs
  end do kkp
  
  call barrier

  call findgntn0_clear

  deallocate(potcl) 
  deallocate(emat12, emat34)
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

      !! Calculate the ou plane wave elements
      ! Pack selected o-u combinations in struct
      ematbc%n1=no
      ematbc%il1=bcouabs%il1
      ematbc%iu1=bcouabs%iu1
      ematbc%n2=nu
      ematbc%il2=bcouabs%il2
      ematbc%iu2=bcouabs%iu2

      !! As long as q=0 this should not be needed
      ! Allocate space for M_{o1o2,G} at fixed (k, q)
      if(allocated(mou)) deallocate(mou)
      allocate(mou(no,nu,numgq))

      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, mou, ematbc)

    end subroutine getmou

    subroutine getmuo(iqnr, iknr, muo)
      integer(4), intent(in) :: iqnr , iknr
      complex(8), allocatable, intent(out) :: muo(:,:,:)

      type(bcbs) :: ematbc

      !! Calculate the uo plane wave elements
      ! Pack selected uo combinations in struct
      ematbc%n1=nu
      ematbc%il1=bcouabs%il2
      ematbc%iu1=bcouabs%iu2
      ematbc%n2=no
      ematbc%il2=bcouabs%il1
      ematbc%iu2=bcouabs%iu1

      !! As long as q=0 this should not be needed
      ! Allocate space for M_{o1o2,G} at fixed (k, q)
      if(allocated(muo)) deallocate(muo)
      allocate(muo(nu,no,numgq))

      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iqnr, iknr, muo, ematbc)

    end subroutine getmuo

    subroutine makeexcli

      complex(8) :: exclit(nou,nou) 

      !!<-- RR and RA part
      j1 = 0           ! ou
      do iu1 = 1, nu   ! u (k+q)
        do io1 = 1, no ! o (k)
          j1 = j1 + 1
          ! emat12_j = M_o1u1ki, M_o2u1ki, ..., M_oNu1ki, M_o1u2ki, ..., M_oNuMki
          emat12(j1, :) = ematouk(io1, iu1, :, iknr)
        end do
      end do
! Try:
! emat12 = reshape(ematouk(:,:,:,ikr),[nou, numgq])
      ! M_ou -> M^*_ou
      emat12 = conjg(emat12)
      !!-->

      !!<-- RR part
      j2 = 0           ! o'u'
      do iu2 =1, nu    ! u' (k'+q)
        do io2 = 1, no ! o' (k')
          j2 = j2 + 1
          ! emat34_j = (M_o1u1kj,M_o2u1kj,...,M_oNu1kj,M_o1u2kj,...,M_oNuMkj)*v(G,0)
          emat34(j2, :) = ematouk(io2, iu2, :, jknr) * potcl(:)
        end do
      end do
      ! Calculate exchange matrix elements: v_{1234} = m_{12}^* m_{34}^t
      ! exclit_{j1, j2} = \Sum_{G} emat12_{j1,G} (emat34^T)_{G,j2}
      !   i.e. exclit_{o u ki, o' u' kj} =
      !          \Sum_{G} M^*_{o u ki}(G,0) M_{o' u' kj}(G,0) v(G,0)
      call zgemm('n', 't', nou, nou, numgq, zone/omega/nkptnr,&
        & emat12, nou, emat34, nou, zzero, exclit, nou)

      ! Map back to individual band indices
      j2 = 0           ! o'u'
      do iu2 = 1, nu   ! u'
        do io2 = 1, no ! o'
          j2 = j2 + 1
          j1 = 0           ! ou
          do iu1 = 1, nu   ! u
            do io1 = 1, no ! o
              j1 = j1 + 1
              ! exclit_{o_j1 u_j1, o'_j2 u'_j2} -> excli_{u_j2 o_j2, u'_j1 o'_j1}
              !   i.e. excli_{u_j2 o_j2 k, u'_j1 o'_j1 k'}
              excli(iu1, io1, iu2, io2) = exclit(j1, j2)
            end do
          end do
        end do
      end do
      !!-->

      if(fcoup == .true.) then
        !!<-- RA part
        j2 = 0           ! u'o'
        do io2 = 1, no   ! o' (k'+q)
          do iu2 =1, nu  ! u' (k')
            j2 = j2 + 1
            ! emat34_j = (M_u1o1kj,M_u2o1kj,...,M_uMo1kj,M_u1o2kj,...,M_uMoNkj)*v(G,0)
            emat34(j2, :) = ematuok(iu2, io2, :, jknr) * potcl(:)
          end do
        end do
        ! Calculate coupling exchange matrix elements: v_{1234} = m_{12}^* m_{34}^t
        ! exclitc_{j1, j2} = \Sum_{G} emat12_{j1,G} (emat34^T)_{G,j2}
        !   i.e. exclitc_{o u ki, u' o' kj} =
        !          \Sum_{G} M^*_{o u ki}(G,0) M_{u' o' kj}(G,0) v(G,0)
        call zgemm('n', 't', nou, nou, numgq, zone/omega/nkptnr,&
          & emat12, nou, emat34, nou, zzero, exclit, nou)
        ! Map back to individual band indices
        j2 = 0           ! u'o'
        do io2 = 1, no   ! o'
          do iu2 = 1, nu ! u'
            j2 = j2 + 1
            j1 = 0           ! ou
            do iu1 = 1, nu   ! u
              do io1 = 1, no ! o
                j1 = j1 + 1
                ! exclit_{o_j1 u_j1, u'_j2 o'_j2} -> exclic_{u_j2 o_j2, u'_j1 o'_j1}
                !   i.e. exclic_{u_j2 o_j2 k, u'_j1 o'_j1 k'}
                exclic(iu1, io1, iu2, io2) = exclit(j1, j2)
              end do
            end do
          end do
        end do
        !!-->
      end if
    end subroutine makeexcli

end subroutine b_exccoulint
!EOC

