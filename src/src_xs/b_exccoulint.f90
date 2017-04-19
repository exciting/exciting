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
  use mod_lattice, only: omega
  use modinput, only: input
  use modmpi, only: rank, barrier, mpi_allgatherv_ifc
  use modxs, only: xsgnt, unitout,&
                 & ngq,&
                 & vqlmt, ivgmt,&
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
  use m_b_ematqk
  use m_putgetbsemat
  use mod_xsgrids
  use mod_Gkvector, only: gkmax

use m_writecmplxparts
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
  integer(4) :: ik, jk, iknr, ikpnr, jknr, jkpnr
  ! Number of occupied/unoccupied states at ik and jk
  integer(4) :: ino, inu, jno, jnu
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

  integer(4) :: igq, igqmt
  character(256) :: m_write, dirname
  logical :: fwp
  logical :: fcoup
  real(8), parameter :: epslat = 1.0d-8
  logical :: fsamek

  fcoup = input%xs%bse%coupling
  fwp = input%xs%bse%writeparts

  !---------------!
  !   main part   !
  !---------------!

  !write(*,*) "Hello, this is b_exccoulint at rank:", mpiglobal%rank

  if(iqmt /= 1 .and. .not. fti) then 
    write(*, '("Error(",a,"):&
     & Finite momentum tranfer currently only supported&
     & when using time reversal symmetry.")') trim(thisnam)
    call terminate
  end if

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
  call xasgauntgen (input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax))
  ! Find indices for non-zero gaunt coefficients in xsgnt,
  ! and creates index maps, e.g. given l1,m1,l2,m2 -> non zero l3,m3
  ! Up l1 up to maxl_mat, l2 up to maxl_mat, l3 up to maxl_e
  ! Note: In the maps and following plane wave matrix elements l1 and l2 correspond
  !       to the APW/LOs while l3 corresponds to the exponetial.
  call findgntn0(maxl_mat, maxl_mat, maxl_e, xsgnt)

  if(rank .eq. 0) then
    write(unitout, '(a,3i8)') 'Info(' // thisnam // '):&
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
  call setranges_modxs(iqmt, input%xs%bse%coupling, input%xs%bse%ti)

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

  call xsgrids_init(vqlmt(1:3, iqmt), gkmax) 

  !! Set vkl0 to k-grid
  !call init1offs(k_kqmtp%kset%vkloff)
  !call xssave0

  ! Change file extension and write out k points
  call genfilname(iqmt=iqmtgamma, dotext='_EXC.OUT', setfilext=.true.)
  if(mpiglobal%rank == 0) then
    call writekpts
  end if

  ! Set vkl to k+qmt-grid
  call init1offs(k_kqmtp%kqmtset%vkloff) ! = init1offs(qvkloff(iqmt))

  if(all(abs(k_kqmtp%kset%vkloff-k_kqmtp%kqmtset%vkloff) < epslat)) then 
    write(unitout, '("Info(b_exccoulint): Same k-grids for iqmt=1 and iqmt=",i3)') iqmt
    fsamek=.true.
  else
    write(unitout, '("Info(b_exccoulint): Different k-grids for iqmt=1 and iqmt=",i3)') iqmt
    fsamek=.false.
  end if

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

  ! Allocate \bar{v}_{G}(qmt)
  allocate(potcl(numgq))
  potcl = 0.d0

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

  allocate(mou(no_bse_max, nu_bse_max, numgq))
  mou = zzero

  if(fra) then 
    allocate(muo(nu_bse_max, no_bse_max, numgq))
  end if

  if(mpiglobal%rank == 0) then
    write(unitout, *)
    write(unitout, '("Info(b_exccoulint):&
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

  !! RR:  M_ouk(G,qmt) = <io ik|e^{-i(qmt+G)r}|iu ikp>
  !!      with ikp = ik+qmt
  do ik = kpari, kparf

    ! Get golbal non reduced k point index
    ! from the BSE k-point index set.
    iknr = kmap_bse_rg(ik)
    ! Get index of ik+qmt
    ikpnr = k_kqmtp%ik2ikqmt(iknr)

    !write(*,*) "iknr =", iknr
    !write(*,*) "ikpnr =", ikpnr

    ! Get the number of participating occupied/unoccupied
    ! states at current k point
    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1

    !write(*,*) "inu =", inu
    !write(*,*) "ino =", ino

    ! Get M_ouk(G,qmt) = <io ik|e^{-i(qmt+G)r}|iu ikp>
    ! with ikp = ik+qmt
    call getmou(mou(1:ino, 1:inu, 1:numgq))

    ! and save it for all ik
    ematouk(1:ino, 1:inu, 1:numgq, ik) = mou(1:ino, 1:inu, 1:numgq)

    if(mpiglobal%rank == 0) then
      write(6, '(a,"Exccoulint - mou progess:", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(ik-kpari+1)/dble(kparf-kpari+1)
      flush(6)
    end if

  end do

  if(mpiglobal%rank == 0) then
    write(6, *)
  end if

  !! RA:  M_uok-qmt(G,qmt) = < ju jkp|e^{-i(qmt+G)r}|jo jk>
  !!      with jkp = jk-qmt
  if(fra) then 

    ! Set vkl0 to k-qmt-grid
    call init1offs(k_kqmtm%kqmtset%vkloff)
    call xssave0

    ! Change file extension and write out k points
    call genfilname(iqmt=iqmt, auxtype='mqmt', dotext='_EXC.OUT', setfilext=.true.)
    if(mpiglobal%rank == 0) then
      call writekpts
    end if

    ! Set vkl to k-grid
    call init1offs(k_kqmtm%kset%vkloff)
    ! Override ikmapikq(1:nkpt,iqmt) to link indices of k-qmt with k 
    ikmapikq(1:nkpt, iqmt) = k_kqmtm%ikqmt2ik(1:nkpt)

    do jk = kpari, kparf

      ! Get golbal non reduced k point index
      ! from the BSE k-point index set.
      jknr = kmap_bse_rg(jk)
      jkpnr = k_kqmtm%ik2ikqmt(jknr) 

      !write(*,*)
      !write(*,*) "jk, jkp=", jknr, jkpnr

      ! Get the number of participating occupied/unoccupied
      ! states at current k point
      jnu = koulims(2,jknr) - koulims(1,jknr) + 1
      jno = koulims(4,jknr) - koulims(3,jknr) + 1

      ! Get M_uok-qmt(G,qmt) = <ju jkp|e^{-i(qmt+G)r}|jo jk>
      ! with jkp = jk-qmt
      call getmuo(muo(1:jnu, 1:jno, 1:numgq))
      ematuok(1:jnu, 1:jno, 1:numgq, jk) = muo(1:jnu, 1:jno, 1:numgq)
    
      !if(mpiglobal%rank == 0) then
      !  write(6, '(a,"Exccoulint - muo progess:", f10.3)', advance="no")&
      !    & achar( 13), 100.0d0*dble(jk-kpari+1)/dble(kparf-kpari+1)
      !  flush(6)
      !end if

    end do
  end if

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

!if(rank == 0) then 
!write(*,*) "ematouk"
!do ik = 1, nk_bse
!  write(*,'(2E13.5)') ematouk(:,:,:,ik)
!end do
!end if

  if(mpiglobal%rank == 0) then
    !write(*,*)
  end if
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
  do igq1 = 1, numgq
    call genwiqggp(0, iqmt, igq1, igq1, potcl(igq1))
  end do
  ! Set Gmt component term of coulomb potential to zero [Ambegaokar-Kohn]
  igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
  potcl(igqmt) = 0.d0
  if(fwp) then 
    call writecmplxparts('Vfc', revec=dble(potcl), dirname='Vfc')
  end if


  if(mpiglobal%rank == 0) then
    write(unitout, *)
    write(unitout, '("Info(b_exccoulint): Generating V matrix elements")')
    write(unitout, '("Info(b_exccoulint): Zeroing Coulomb potential at G+qmt index:", i3)') igqmt
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
    !write(*,*)
    !write(*,'(a, i4)') "ikkp =", ikkp

    iknr = kmap_bse_rg(ik)
    jknr = kmap_bse_rg(jk) 

    !write(*,'(a, i4)') "iknr =", iknr
    !write(*,'(a, i4)') "jknr =", jknr

    ! Get index of ik+qmt
    ikpnr = k_kqmtp%ik2ikqmt(iknr)
    ! Get index of jk+qmt
    jkpnr = k_kqmtp%ik2ikqmt(jknr)
    if(fra) then 
      jkpnr = k_kqmtm%ik2ikqmt(jknr)
    end if

    ! Get number of transitions at ik,jk
    inou = kousize(iknr)
    jnou = kousize(jknr)

    call makeexcli(excli(1:inou,1:jnou))

    ! Parallel write
    if(fra) then
      if(fwp) then
        call writecmplxparts('Vra', dble(excli(1:inou,1:jnou)),&
          & aimag(excli(1:inou,1:jnou)), ik1=iknr, ik2=jknr, dirname='Vra')
      end if
      call b_putbsemat(exclifname, 78, ikkp, iqmt, excli)
    else
      if(fwp) then
        call writecmplxparts('Vrr', dble(excli(1:inou,1:jnou)),&
          & aimag(excli(1:inou,1:jnou)), ik1=iknr, ik2=jknr, dirname='Vrr')
      end if
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
    write(6,*)
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

  contains 

    subroutine getmou(mou)
      complex(8), intent(out) :: mou(:,:,:)

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save

      !write(*,*)
      !write(*,*) "getmou:"
      !write(*,*) "  Mou"

      !write(*,*) "shape(mou)", shape(mou)

      fileext0_save = filext0
      fileext_save = filext

      !----------------------------------------------------------------!
      ! Calculate M_{io iu ik}(G, qmt) = <io ik|e^{-i(qmt+G)r}|iu ikp> !
      ! with ikp = ik+qmt                                              !
      !----------------------------------------------------------------!

      ! Bands
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=inu
      ematbc%il2=koulims(1,iknr)
      ematbc%iu2=koulims(2,iknr)

      usefilext0 = .true.
      ! Set EVECFV_QMT001.OUT as bra state file
      iqmt0 = iqmtgamma
      call genfilname(iqmt=iqmt0, setfilext=.true.)
      filext0 = filext
      !write(*,*) "filext0 =", trim(filext0)

      if(fsamek) then 
        ! Set EVECFV_QMT001.OUT as ket state file
        iqmt1 = iqmtgamma
      else
        ! Set EVECFV_QMTXYZ.OUT as ket state file
        iqmt1 = iqmt
      end if
      call genfilname(iqmt=iqmt1, setfilext=.true.)
      !write(*,*) "filext =", trim(filext)

      emat_ccket=.false.
      ! Set up ikmapikq to link (ik,iqmt) to (ikp)
      ikmapikq_ptr => ikmapikq
      ! Set vkl0_ptr,... to k-grid and vkl1_ptr, ... to k+qmt-grid
      call setptr01
      ! Calculate M_{io iu,G}(ik, qmt)
      if (input%xs%bse%xas) then
        call b_ematqk_core(iqmt, iknr, mou, ematbc, 'ou')
      else
        call b_ematqk(iqmt, iknr, mou, ematbc)
      end if
      !------------------------------------------------------------------!
      if(.false.) then 
        if(fwp) then
          do igq=1,numgq
            call genfilname(iqmt=iqmt, iq=igq, dotext='', fileext=m_write)
            m_write='Mou'//trim(adjustl(m_write))
            call genfilname(iqmt=iqmt, dotext='', fileext=dirname)
            dirname='Mou_exc'//trim(adjustl(dirname))
            call writecmplxparts(trim(adjustl(m_write)), remat=dble(mou(:,:,igq)),&
              & immat=aimag(mou(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
          end do
        end if
      end if

      filext0 = fileext0_save
      filext = fileext_save

      !write(*,*)
      !write(*,*)

    end subroutine getmou

    subroutine getmuo(muo)
      complex(8), intent(out) :: muo(:,:,:)

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save

      !!write(*,*)
      !!write(*,*) "getmuo:"
      !!write(*,*) "  Muo"

      fileext0_save = filext0
      fileext_save = filext
      
      !-----------------------------------------------------------------!
      ! Calculate M_{ju jo jkp}(G, qmt) = <ju jkp|e^{-i(qmt+G)r}|jo jk> !
      ! with jkp = jk-qmt                                               !
      !-----------------------------------------------------------------!

      ! Bands
      ematbc%n1=jnu
      ematbc%il1=koulims(1,jknr)
      ematbc%iu1=koulims(2,jknr)
      ematbc%n2=jno
      ematbc%il2=koulims(3,jknr)
      ematbc%iu2=koulims(4,jknr)

      usefilext0 = .true.
      if(iqmt /= 1) then 
        ! Set EVECFV_QMTXYZ_mqmt.OUT as bra state file
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt, auxtype='mqmt', setfilext=.true.)
      else
        ! Set EVECFV_QMT001.OUT as bra state file
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true.)
      end if
      filext0 = filext
      !write(*,*) "filext0 =", trim(filext0)

      ! Set EVECFV_QMT001.OUT as ket state file
      iqmt1 = iqmtgamma
      call genfilname(iqmt=iqmt1, setfilext=.true.)
      !write(*,*) "filext =", trim(filext)

      emat_ccket=.false.
      ! Note: Currently working only for qmt=0
      ! Set up ikmapikq to link (jkp,iqmt) to (jk)
      ikmapikq_ptr => ikmapikq
      ! Set vkl0_ptr k-qmt-grid, vkl1_ptr, ... to k-grid
      call setptr01
      ! Calculate M_{ju jo,G}(jkp, qmt)
      if (input%xs%bse%xas) then
        call b_ematqk_core(iqmt, jkpnr, muo, ematbc, 'uo')
      else
        call b_ematqk(iqmt, jkpnr, muo, ematbc)
      end if
      !------------------------------------------------------------------!
      if(.false.) then 
        if(fwp) then
          do igq=1,numgq
            call genfilname(iqmt=iqmt, iq=igq, dotext='', fileext=m_write)
            m_write='Muo'//trim(adjustl(m_write))
            call genfilname(iqmt=iqmt, dotext='', fileext=dirname)
            dirname='Muo_exc'//trim(adjustl(dirname))
            call writecmplxparts(trim(adjustl(m_write)), remat=dble(muo(:,:,igq)),&
              & immat=aimag(muo(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
          end do
        end if
      end if

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

