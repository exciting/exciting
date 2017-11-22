! Copyright(C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: exccoulint
! !INTERFACE:
subroutine exccoulint(iqmt)
! !USES:
  use mod_misc, only: filext
  use mod_constants, only: zone, zzero
  use mod_APW_LO, only: lolmax
  use mod_lattice, only: omega
  use modinput, only: input
  use modmpi, only: rank, barrier, mpi_allgatherv_ifc
  use modxs, only: xsgnt, unitout,&
                 & ngq, nqmt,&
                 & totalqlmt, ivgmt,&
                 & kpari, kparf,&
                 & ppari, pparf,& 
                 & bcbs, iqmtgamma,&
                 & filext0, usefilext0, iqmt0, iqmt1, ivgigq
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genfilname
  use m_getunit
  use modbse
  use m_ematqk
  use m_putgetbsemat
  use mod_xsgrids
  use mod_Gkvector, only: gkmax
! !DESCRIPTION:
!   Calculates the exchange term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Forked from exccoulint.F90 and adapted for non-TDA and Q-dependent BSE. (Aurich)
!EOP
!BOC      

  implicit none

  ! I/O

  integer, intent(in) :: iqmt ! Index of momentum transfer Q

  ! Local variables

  character(*), parameter :: thisname = 'exccoulint'
  ! ik,jk block of V matrix (final product)
  complex(8), allocatable :: excli(:, :)
  ! Auxilliary arrays for the construction of excli
  complex(8), allocatable :: ematuok(:, :, :, :)
  complex(8), allocatable :: muo(:,:,:)
  ! Truncated coulomb potential
  real(8), allocatable :: potcl(:)
  ! ik jk q points 
  integer(4) :: ikkp
  integer(4) :: ik, iknr, ikpnr, ikmnr
  integer(4) :: jk, jknr, jkpnr, jkmnr
  integer(4), allocatable, target :: ikm2ikp_dummy(:,:)
  ! Number of occupied/unoccupied states at ik and jk
  integer(4) :: ino, inu
  ! Number of transitions at ik and jk
  integer(4) :: inou, jnou
  ! State loop indices
  integer(4) :: io, jo, iu, ju
  ! Aux.
  integer(4) :: igq1, numgq
  ! Maximal l used in the APWs and LOs
  !   Influences quality of plane wave matrix elements
  integer(4) :: maxl_apwlo
  ! Maximal l used in the Reghley expansion of exponential
  !   Influences quality of plane wave matrix elements
  integer(4) :: maxl_e
  ! Maximal l used in the APWs and LOs in the groundstate calculations
  !   Influences quality of eigencoefficients.
  integer(4) :: maxl_mat
  ! Timinig vars
  real(8) :: tpw1, tpw0
  integer(4) :: flg_analytic

  integer(4) :: igqmt
  logical :: fcoup
  real(8), parameter :: epslat = 1.0d-8
  logical :: fsamekp, fsamekm
  logical :: fchibarq

  fcoup = input%xs%bse%coupling
  fchibarq = input%xs%bse%chibarq

  !---------------!
  !   main part   !
  !---------------!

  ! General setup
  call init0
  ! k-point setup
  ! Also allocated the radial functions (mod_APW_LO)
  call init1
  ! Save variables of the unshifted (apart from xs:vkloff) k grid 
  ! to modxs (vkl0, ngk0, ...)
  call xssave0
  ! q-point and qmt-point setup
  !   Init 2 sets up (task 441):
  !   * A list of momentum transfer vectors form the q-point list 
  !     (modxs::vqmtl and mod_qpoint::vql)
  !   * Offset of the k+qmt grid derived from k offset an qmt point (modxs::qvkloff)
  !   * non-reduced mapping between ik,qmt and ik' grids (modxs::ikmapikq)
  !   * G+qmt quantities (modxs)
  !   * The square root of the Coulomb potential for the G+qmt points
  !   * Reads STATE.OUT
  !   * Generates radial functions (mod_APW_LO)
  call init2

  ! Generate gaunt coefficients used in the construction of 
  ! the plane wave matrix elements in ematqk.

  ! gs%lmaxapw defaults to 8, BUT xs%lmaxapw defaults to 10 
  ! in xsinit gs%lmaxapw is overwritten with the xs%lmaxapw value.
  ! --> xs%lmaxapw = 10
  maxl_apwlo = max(input%groundstate%lmaxapw, lolmax)

  ! xs%lmaxapwwf defaults to -1, in which case it is set to gs%lmaxmat by init2
  ! which in trun was previously overwritten by xs%lmaxmat by xsinit. 
  ! xs%lmatmax defaults to 5.
  ! --> xs%lmaxapwwf = xs%lmaxmat = 5
  maxl_mat = max(input%xs%lmaxapwwf, lolmax)

  ! xs%lmaxemat defaults to 3
  maxl_e = input%xs%lmaxemat

  ! Allocates xsgnt array of shape 
  ! ((maxl_apwlo+1)**2, (maxl_e+1)**2, (maxl_apwlo+1)**2)
  ! and fills it with the corresponding gaunt coefficients
  ! Note: In the generation of the gaunt coefficients l1 and l3 correspond
  !       to the APW/LOs while l2 corresponds to the exponetial.
  call xsgauntgen(maxl_apwlo, maxl_e, maxl_apwlo)
  ! Find indices for non-zero gaunt coefficients in xsgnt,
  ! and creates index maps, e.g. given l1,m1,l2,m2 -> non zero l3,m3
  ! Up l1 up to maxl_mat, l2 up to maxl_mat, l3 up to maxl_e
  ! Note: In the maps and following plane wave matrix elements l1 and l2 correspond
  !       to the APW/LOs while l3 corresponds to the exponetial.
  call findgntn0(maxl_mat, maxl_mat, maxl_e, xsgnt)

  if(rank .eq. 0) then
    write(unitout, '(a,3i8)') 'Info(' // thisname // '):&
      & Gaunt coefficients generated within lmax values:',&
      & input%groundstate%lmaxapw, input%xs%lmaxemat, input%groundstate%lmaxapw
    call flushifc(unitout)
  end if

  ! Read Fermi energy from file
  ! Use EFERMI_QMT001.OUT
  call genfilname(iqmt=iqmtgamma, setfilext=.true.)
  call readfermi

  ! Set ist* variables and ksgap in modxs using findocclims
  ! This also reads in 
  ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
  ! modxs:evalsv0, modxs:occsv0
  call setranges_modxs(iqmt)

  ! Select relevant transitions for the construction
  ! of the BSE hamiltonian
  ! Also sets nkkp_bse, nk_bse 
  call select_transitions(iqmt, serial=.false.)

  ! Write support information to file
  if(mpiglobal%rank == 0) then
    if(fchibarq) then 
      call genfilname(basename=trim(infofbasename)//'_'//trim(exclifbasename),&
        &  bsetype='-BAR', iqmt=iqmt, filnam=infofname)
      call putbseinfo(infofname, iqmt)
    else
      call genfilname(basename=trim(infofbasename)//'_'//trim(exclifbasename),&
        & iqmt=iqmt, filnam=infofname)
      call putbseinfo(infofname, iqmt)
    end if
  end if

  ! Set output file names
  if(fchibarq) then
    call genfilname(basename=exclifbasename, bsetype='-BAR', iqmt=iqmt, filnam=exclifname)
  else
    call genfilname(basename=exclifbasename, iqmt=iqmt, filnam=exclifname)
  end if

  write(unitout, '("Info(",a,"): Size of file ",a," will be about ", f12.6, " GB" )')&
    & trim(thisname), trim(exclifbasename),&
    & int(nou_bse_max,8)**2*int(nkkp_bse,8)*16.0d0/1024.0d0**3
  call flushifc(unitout)

  call xsgrids_init(totalqlmt(1:3, iqmt), gkmax) 

  ! Change file extension and write out k points
  call genfilname(iqmt=iqmtgamma, dotext='_EXC.OUT', setfilext=.true.)
  if(mpiglobal%rank == 0) then
    call writekpts
  end if

  !--------------------------------------------------------------------------------!
  ! Setup k grids
  !--------------------------------------------------------------------------------!
  ! Save the k,G arrays of the k-qmt/2 grid to modxs::vkl0 etc
  call init1offs(k_kqmtm%kqmtset%vkloff)
  call xssave0
  ! Save the k,G arrays of the k+qmt/2 grid to the default locations
  call init1offs(k_kqmtp%kqmtset%vkloff)
  ! Check whether k+-qmt/2 grids are identical to k grid
  if(all(abs(k_kqmtp%kqmtset%vkloff-k_kqmtp%kset%vkloff) < epslat)) then 
    write(unitout, '("Info(exccoulint):&
      & k+qmt/2-grid is identical to iqmt=1 grid for, iqmt=",i3)') iqmt
    fsamekp=.true.
  else
    fsamekp=.false.
  end if
  if(all(abs(k_kqmtm%kqmtset%vkloff-k_kqmtm%kset%vkloff) < epslat)) then 
    write(unitout, '("Info(exccoulint):&
      & k-qmt/2-grid is identical to iqmt=1 grid for, iqmt=",i3)') iqmt
    fsamekm=.true.
  else
    fsamekm=.false.
  end if

  ! Make link map connecting ik-qmt/2 to ik+qmt/2
  ! ( only the iqmt component will be used ...)
  if(allocated(ikm2ikp_dummy)) deallocate(ikm2ikp_dummy)
  allocate(ikm2ikp_dummy(k_kqmtp%kset%nkptnr, nqmt))
  ikm2ikp_dummy=0
  ikm2ikp_dummy(:, iqmt) = ikm2ikp(:)

  !--------------------------------------------------------------------------------!

  ! Change file extension and write out k points
  call genfilname(iqmt=iqmt, dotext='_EXC.OUT', setfilext=.true.)
  if(mpiglobal%rank == 0) then
    call writekpts
  end if

  ! Number of G+qmt points 
  numgq = ngq(iqmt) 

  ! Calculate radial integrals used in the construction 
  ! of the plane wave matrix elements for exponent (qmt+G)
  call ematrad(iqmt)

  ! Allocate \bar{v}_{G}(qmt) (or v_{G}(qmt))
  allocate(potcl(numgq))
  potcl = 0.d0

  ! Allocate eigenvalue/eigenvector related
  ! quantities for use in ematqk
  call ematqalloc

  ! Work array to store result of plane wave matrix elements
  ! for each considered k point 
  allocate(ematuok(nu_bse_max, no_bse_max, numgq, nk_bse))
  allocate(muo(nu_bse_max, no_bse_max, numgq))
  muo = zzero

  ! Info out
  if(mpiglobal%rank == 0) then
    write(unitout, *)
    write(unitout, '("Info(exccoulint):&
     & Generating plane wave matrix elements for momentum transfer iqmt=",i4)') iqmt
    call timesec(tpw0)
  end if

  !! Plane wave matrix elements calculation.
  !---------------------------!
  !     Loop over k-points    !
  !---------------------------!
  ! Parallelize over non reduced k-points
  ! participating in the BSE
  call genparidxran('k', nk_bse)

  !! RR:  M_uok(G,qmt) = <iu ikm|e^{-i(G+qmt)r}|io ikp>
  !!      with ikp = ik+qmt
  do ik = kpari, kparf

    ! Get golbal non reduced k point index
    ! from the BSE k-point index set.
    iknr = kmap_bse_rg(ik)

    ! Get index of ik+qmt/2 and ik-qmt/2
    ikpnr = k_kqmtp%ik2ikqmt(iknr)
    ikmnr = k_kqmtm%ik2ikqmt(iknr)

    ! Get the number of participating occupied/unoccupied
    ! states at current k point
    ! Note: The saved ranges refer to the to k associated k_- points for the 
    !       unoccupied, and to the k_+ points for the occupied states
    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1

    ! Get M_uokm(G,qmt) = <iu ikm|e^{-i(G+qmt)r}|io ikp>
    ! with ikp = ik+qmt/2 and ikm = ik-qmt/2
    call getmuo(muo(1:inu, 1:ino, 1:numgq))

    ! and save it for all ik
    ematuok(1:inu, 1:ino, 1:numgq, ik) = muo(1:inu, 1:ino, 1:numgq)

    if(mpiglobal%rank == 0) then
      write(6, '(a,"Exccoulint - muo progress:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(ik-kpari+1)/dble(kparf-kpari+1)
      flush(6)
    end if

  end do

  if(mpiglobal%rank == 0) then
    write(6, *)
  end if

  ! Helper no longer needed
  if(allocated(muo)) deallocate(muo)

  ! Communicate array-parts wrt. k-points
  call mpi_allgatherv_ifc(set=nk_bse,&
    & rlen=nu_bse_max*no_bse_max*numgq, zbuf=ematuok,&
    & inplace=.true., comm=mpiglobal)

  if(mpiglobal%rank == 0) then
    call timesec(tpw1)
    write(unitout, '("  Timing (in seconds)	   :", f12.3)') tpw1 - tpw0
  end if

  !! Generation of V matrix.
  !-------------------------------!
  !     Loop over(k,kp) pairs     !
  !-------------------------------!
  call genparidxran('p', nkkp_bse)

  allocate(excli(nou_bse_max, nou_bse_max))

  ! Set up coulomb potential
  ! Construct it via v^{1/2}(G,qmt)*v^{1/2}(G,qmt),
  ! which corresponds to the first passed flag=0.
  ! Use all G for which |G+qmt|<gqmax
  ! Use cutoff if requested
  select case(trim(input%xs%bse%cuttype))
    case("none")
      flg_analytic = 0
    case("0d")
      flg_analytic = 4
      write(unitout, '("Info(exccoulint):&
        & Applying Coulomb cutoff for low dimensional systems. Type:",a)') &
        & trim(input%xs%bse%cuttype)
    case("2d")
      flg_analytic = 5
      write(unitout, '("Info(exccoulint):&
        & Applying Coulomb cutoff for low dimensional systems. Type:",a)') &
        & trim(input%xs%bse%cuttype)
    case default
      write(*,*) "Error(exccoulint): Invalid cuttype"
      call terminate
  end select
  do igq1 = 1, numgq
    call genwiqggp(flg_analytic, iqmt, igq1, igq1, potcl(igq1))
  end do

  ! If Q=0 compute \bar{P}, so use the truncated Coulomb potential.
  ! If Q/=0 and TDA compute \bar{P}, if chibarq = .true. (default)
  ! If Q/=0 and non-TDA compute \chi, i.e. no zeroing of the Coulomb potential.
  if(iqmt==1 .or. input%xs%bse%chibarq) then 
    igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
    potcl(igqmt) = 0.d0
  end if

  if(mpiglobal%rank == 0) then
    write(unitout, *)
    write(unitout, '("Info(exccoulint): Generating V matrix elements")')
    if(iqmt == 1 .or. input%xs%bse%chibarq) then 
      write(unitout, '("Info(exccoulint): Zeroing Coulomb potential at G+qmt index:", i3)') igqmt
    end if
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

    ! Get index of ik+qmt/2 and ik-qmt/2
    ikpnr = k_kqmtp%ik2ikqmt(iknr)
    ikmnr = k_kqmtm%ik2ikqmt(iknr)

    ! Get index of jk+qmt/2 and jk-qmt/2
    jkpnr = k_kqmtp%ik2ikqmt(jknr)
    jkmnr = k_kqmtm%ik2ikqmt(jknr)

    ! Get number of transitions at reference k-points
    ! ik, jk for current qmt 
    inou = kousize(iknr)
    jnou = kousize(jknr)

    call makeexcli(excli(1:inou,1:jnou))

    ! Parallel write
    call putbsemat(exclifname, 77, ikkp, iqmt, excli)

    if(mpiglobal%rank == 0) then
      write(6, '(a,"Exccoulint progess:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(ikkp-ppari+1)/dble(pparf-ppari+1)
      flush(6)
    end if

  ! End loop over(k,kp) pairs
  end do kkp

  if(mpiglobal%rank == 0) then
    write(6,*)
  end if
  
  call barrier(callername=trim(thisname))

  if(mpiglobal%rank == 0) then
    call timesec(tpw1)
    write(unitout, '("  Timing (in seconds)	   :", f12.3)') tpw1 - tpw0
  end if

  call findgntn0_clear
  call xsgrids_finalize()

  deallocate(potcl) 
  deallocate(ikm2ikp_dummy)
  if(allocated(ematuok)) deallocate(ematuok)
  deallocate(excli)

  contains 

    subroutine getmuo(muo)
      use mod_variation, only: ematqk_sv
      complex(8), intent(out) :: muo(:,:,:)

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save

      fileext0_save = filext0
      fileext_save = filext
      
      !------------------------------------------------------------------!
      ! Calculate M_{iu io ikm}(G, qmt) = <iu ikm|e^{-i(G+qmt)r}|io ikp> !
      ! with ikp = ik-qmt                                                !
      !------------------------------------------------------------------!

      ! Bands
      ! Note: The saved ranges refer to the to k associated k_- points for the 
      !       unoccupied, and to the k_+ points for the occupied states
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=ino
      ematbc%il2=koulims(3,iknr)
      ematbc%iu2=koulims(4,iknr)

      usefilext0 = .true.
      if(fsamekm) then 
        ! Set EVECFV_QMT001.OUT as bra state file
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true., fileext=filext0)
      else
        ! Set EVECFV_QMTXYZ_m.OUT as bra state file (k-qmt/2 grid)
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt0, auxtype="m", fileext=filext0)
      end if

      ! Set EVECFV_QMTXYZ.OUT as ket state file (k+qmt/2 grid)
      iqmt1 = iqmt
      ! Set to EVECFV_QMT001.OUT if it is the same grid
      if(fsamekp) iqmt1 = iqmtgamma
      call genfilname(iqmt=iqmt1, setfilext=.true.)

      ! Use normal plane wave routine that calculates M 
      emat_ccket=.false.

      ! Set up ikmapikq to link (jkm,iqmt) to (jkp)
      ikmapikq_ptr => ikm2ikp_dummy

      ! Set vkl0_ptr k-qmt/2-grid, vkl1_ptr, ... to k+qmt/2-grid
      call setptr01()

      ! Calculate M_{iu io,G}(ikm, qmt)
      if (input%xs%bse%xas) then
        call xasgauntgen (input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
        call ematqk_core(iqmt, ikmnr, muo, ematbc, 'uo')

      else
        if (.not. (input%groundstate%tevecsv)) then ! 1st variation
          call ematqk(iqmt, ikmnr, muo, ematbc)
        else                                        ! 2nd variation 
          call ematqk_sv(iqmt, ikmnr, muo, ematbc)  
        end if
        
      end if
      !------------------------------------------------------------------!

      filext0 = fileext0_save
      filext = fileext_save

    end subroutine getmuo

    subroutine makeexcli(excli)

      complex(8), intent(out) :: excli(inou, jnou) 

      ! Work arrays
      complex(8) :: emat12(inou, numgq), emat34(jnou, numgq)
      integer(4) :: iaoff, jaoff, ia, ja

      ! Offset in combined index for ik and jk
      iaoff = sum(kousize(1:iknr-1))
      jaoff = sum(kousize(1:jknr-1))

      !! RR 
      do ia = 1, inou
        ! Get iu index relative to iu range of ik
        iu = smap_rel(1, ia+iaoff) ! iu (ik-qmt/2)
        ! Get io index relative to io range of ik
        io = smap_rel(2, ia+iaoff) ! io (ik+qmt/2)
        ! emat12_ia = M_u1o1ki, M_u2o1ki, ..., M_uMo1ki, M_uMo2ki, ..., M_uMoNki
        emat12(ia, :) = ematuok(iu, io, :, ik)
      end do
      ! M_uo -> M^*_uo
      emat12 = conjg(emat12)

      do ja = 1, jnou
        ! Get ju index relative to iu range of jk
        ju = smap_rel(1, ja+jaoff) ! ju (jk-qmt/2)
        ! Get jo index relative to io range of jk
        jo = smap_rel(2, ja+jaoff) ! jo (jk+qmt/2)
        ! emat34_ja = (M_u1o1kj,M_u2o1kj,...,M_uMo1kj,M_u1o2kj,...,M_uMoMkj)*v
        emat34(ja, :) = ematuok(ju, jo, :, jk) * potcl(:)
      end do

      ! Calculate exchange matrix elements: 
      ! excli_{ia, ja} = \Sum_{G} emat12_{ia,G} (emat34^T)_{G,ja}
      !   i.e. excli_{iu io ik, ju jo jk}(qmt) =
      !          \Sum_{G} M^*_{iu io ikm}(G,qmt) M_{ju jo jkm}(G,qmt) v(G,qmt)
      call zgemm('n', 't', inou, jnou, numgq, zone/omega/nk_bse,&
        & emat12, inou, emat34, jnou, zzero, excli, inou)

    end subroutine makeexcli

end subroutine exccoulint
!EOC

