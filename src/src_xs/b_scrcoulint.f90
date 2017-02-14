! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: b_scrcoulint
! !INTERFACE:
subroutine b_scrcoulint(iqmt, fra, fti)
! !USES:
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi
  use mod_constants, only: zzero, zone, fourpi
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: iqmap, ngridq, vql, vqc, nqpt, ivq, wqpt
  use mod_kpoint, only: nkptnr, ivknr
  use mod_lattice, only: omega
  use mod_symmetry, only: maxsymcrys
  use modxs, only: xsgnt, unitout,&
                 & ngqmax,&
                 & nqptr, qpari, qparf, ivqr,&
                 & ngq, ppari, pparf, iqmapr,&
                 & vqlr, vqcr, wqptr, ngqr, vqlmt,&
                 & bsedl, bsedu, bsedd,&
                 & bsed, bcbs, ematraddir, eps0dirname,&
                 & filext0, usefilext0, iqmt0, iqmt1,&
                 & iqmtgamma
  use m_xsgauntgen
  use m_findgntn0
  use m_writevars
  use m_genfilname
  use m_getunit
  use m_b_ematqk
  use m_b_ematqk2
  use m_putgetbsemat
  use modbse
  use mod_xsgrids
  use mod_Gkvector, only: gkmax

use m_writecmplxparts
! !DESCRIPTION:
!   Calculates the direct term of the Bethe-Salpeter Hamiltonian.
!
! !REVISION HISTORY:
!   Created June 2008 (S. Sagmeister)
!   Addition of explicit energy ranges for states below and above the Fermi
!      level for the treatment of core excitations (using local orbitals).
!      October 2010 (Weine Olovsson)
!   Forked from scrcoulint.F90 and adapted for non TDA BSE. (Aurich)
!EOP
!BOC      

  implicit none

  ! I/O
  integer, intent(in) :: iqmt ! Index of momentum transfer Q
  logical, intent(in) :: fra  ! Construct RA coupling block
  logical, intent(in) :: fti  ! Use time inverted anti-resonant basis

  ! Local variables
  character(*), parameter :: thisnam = 'b_scrcoulint'

  ! ik,jk block of W matrix (final product)
  complex(8), allocatable :: sccli(:,:)

  ! Auxilliary arrays for the construction of sccli
  complex(8), allocatable :: sccli_t1(:, :), sccli_t2(:, :, :, :)
  complex(8), allocatable :: zm(:,:)
  ! Plane wave arrays
  complex(8), allocatable :: muu(:, :, :), cmuu(:, :)
  complex(8), allocatable :: moo(:, :, :), cmoo(:, :)
  complex(8), allocatable :: mou(:, :, :), cmou(:, :)
  complex(8), allocatable :: muo(:, :, :), cmuo(:, :)
  ! The fourier coefficients of the screend coulomb potential
  complex(8), allocatable :: wfc(:, :)
  ! Auxilliary arrays/vars for the construction of wfc
  complex(8), allocatable :: scieffg(:, :, :)
  complex(8), allocatable :: phf(:, :)
  ! Symmerty maps creation 
  integer(4) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
  integer(4), allocatable :: igqmap(:)
  integer(4) :: jsym, jsymi
  integer(4) :: nsc, ivgsym(3)
  logical :: tphf
  ! Mappings of jk ik combinations to q points
  integer(4) :: ikkp, iknr, jknr, jmknr, jmkpnr, jmkpnr2, ikpnr, jkpnr, ik, jk
  integer(4) :: iqrnr, iqr, iq
  real(8) :: vqr(3), vq(3)
  integer(4) :: numgq
  integer(4) :: iv(3)
  logical :: tq0
  ! Number of occupied/unoccupied states at ik and jk
  integer(4) :: ino, inu, jno, jnu
  ! Number of transitions at ik and jk
  integer(4) :: inou, jnou
  ! Number of (o/u)_i (o/u)_j combinations
  integer(4) :: noo, nuu, nou, nuo
  ! State loop indices
  integer(4) :: io, jo, iu, ju
  ! Combinded loop indices 
  integer(4) :: jaoff, iaoff, ia, ja
  ! Maximal l used in the APWs and LOs
  !   Influences quality of plane wave matrix elements
  integer(4) :: maxl_apwlo
  ! Maximal l used in the Reghley expansion of exponential
  !   Influences quality of plane wave matrix elements
  integer(4) :: maxl_e
  ! Maximal l used in the APWs and LOs in the groundstate calculations
  !   Influences quality of eigencoefficients.
  integer(4) :: maxl_mat
  ! Aux.
  integer(4) :: j1, j2, ii,i,jj,j,igq
  complex(8) :: pref 
  ! Timing vars
  real(8) :: tscc1, tscc0

  ! Auxilliary strings
  character(256) :: syscommand, fileext_scr_read, fileext_ematrad_write
  character(256) :: wfc_write, dirname

  ! External functions
  integer, external :: idxkkp
  logical, external :: tqgamma

  real(8) :: vqoff(3)

  logical :: fcoup
  logical :: fwp

  fcoup = input%xs%bse%coupling
  fwp = input%xs%bse%writeparts

  !---------------!
  !   main part   !
  !---------------!

write(*,*) "Hello, this is b_scrcoulint at rank:", rank

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

  ! General setup
  call init0
  ! k-point setup
  call init1
  ! Save variables of the unshifted (apart from xs:vkloff) k grid 
  ! to modxs (vkl0, ngk0, ...)
  call xssave0
  ! q-point and qmt-point setup
  !   Init 2 sets up (task 440):
  !   * A list of momentum transfer vectors form the q-point list (modxs::vqmtl)
  !   * The reduced unshifted q-grid (modxs::vqlr etc)
  !   * The non-reduced unshifted q-grid (mod_qpoint::vql etc)
  !   * Offset of the k+q grid derived from k offset an q point (modxs::qvkloff)
  !     which is just equal to the k-offset, since it is the unshifted q grid
  !   * non-reduced mapping between ik,q and ik' grids (modxs::ikmapikq)
  !   * G+q quantities for the non-reduced unshifted q-grid (modxs)
  !   * The square root of the Coulomb potential for the non-reduced q points
  !   * Reads STATE.OUT
  !   * Generates radial functions (mod_APW_LO)
  call init2

  ! Making folder for the radial integals pertaining to the plane wave matrix elements
  ematraddir = 'EMATRAD'
  if(rank == 0) then 
    syscommand = '[[ ! -e '//trim(adjustl(ematraddir))//' ]] && mkdir '//trim(adjustl(ematraddir))
    call system(trim(adjustl(syscommand)))
  end if

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

  if(mpiglobal%rank == 0) then
    write(unitout, '(a)') 'Info(' // thisnam // '):&
      & Gaunt coefficients generated within lmax values:'
    write(unitout, '(a, i8)') "lmax1 = lmaxapw =", input%groundstate%lmaxapw
    write(unitout, '(a, i8)') "lmax2 = lmaxemat=", input%xs%lmaxemat
    write(unitout, '(a, i8)') "lmax3 = lmaxapw =", input%groundstate%lmaxapw
    call flushifc(unitout)
  end if

  ! Read Fermi energy from file EFERMI
  ! Use EFERMI_QMT001.OUT (corresponding to the xs groundstate run for the unshifted k grid)
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
      if(fti) then 
        call genfilname(basename=trim(infofbasename)//'_'//trim(scclictifbasename),&
          & iqmt=iqmt, filnam=infofname)
        call b_putbseinfo(infofname, iqmt)
      else
        call genfilname(basename=trim(infofbasename)//'_'//trim(scclicfbasename),&
          & iqmt=iqmt, filnam=infofname)
        call b_putbseinfo(infofname, iqmt)
      end if
    else
      call genfilname(basename=trim(infofbasename)//'_'//trim(scclifbasename),&
        & iqmt=iqmt, filnam=infofname)
      call b_putbseinfo(infofname, iqmt)
    end if
  end if

  ! Set output file name
  if(fra) then
    if(fti) then
      call genfilname(basename=scclictifbasename, iqmt=iqmt, filnam=scclifname)
    else
      call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=scclifname)
    end if
  else
    call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=scclifname)
  end if

  !------------------------------------!
  ! GENERATE FOURIER COEFFICIENTS OF W !     
  ! (and radial integrals for emat)    !
  !------------------------------------!

  !! q-grid setup

  ! Setup reduced q-points in modxs

  !write(*,*)
  call xsgrids_init(vqlmt(1:3,iqmt), gkmax)
  if(fra) then 
    !write(*,*) "Generating reduced q-grid for W Fourier coefficients for RA coupling block."
    if(fti) then 
      !write(*,*) "  Using time inverted anti-resonant basis."
      vqoff = p_pqmtp%pset%vkloff
    else
      !write(*,*) "  Using standard anti-resonant basis."
      vqoff = q_qmtm%qset%vkloff
    end if
  else
    !write(*,*) "Generating reduced q-grid for W Fourier coefficients for RR block."
    vqoff =  q_q%qset%vkloff
  end if

  !write(*,*)
  !write(*,'(a,3E10.3)') "vqoff = ", vqoff

  ! Make reduced q-grid with possible offset
  ! This sets up also G+q quantities and the square root of the Coulomb potential
  ! but the second call below will override these
  call init2offs(vqoff, input%xs%reduceq)
  ! Copy results q-ponit results form mod_qpoint into modxs variables
  nqptr = nqpt
  ivqr = ivq
  iqmapr = iqmap
  vqlr = vql
  vqcr = vqc
  wqptr = wqpt

  if(mpiglobal%rank == 0) then
    write(unitout, '(a, i6)') 'Info(' // thisnam // '): Number of reduced q-points: ', nqptr
    write(unitout,*)
    call flushifc(unitout)
  end if

  !write(*,*)
  !write(*,*) "iqr vqlr"
  do iq = 1, nqptr
    !write(*,'(i3, 3E10.3)') iq, vqlr(1:3,iq)
  end do
    
  ! Make non-reduced q-grid with possible offset
  ! This sets up also G+q quantities and the square root of the Coulomb potential
  call init2offs(vqoff, .false.)

  !write(*,*)
  !write(*,*) "iq vql"
  do iq = 1, nqpt
    !write(*,'(i3, 3E10.3)') iq, vql(1:3,iq)
  end do

  ! Make also ngqr
  if(allocated(ngqr)) deallocate(ngqr)
  allocate(ngqr(nqptr))
  ngqr=0
  do iq= 1, nqpt
    iqr = iqmapr(ivq(1,iq), ivq(2,iq), ivq(3,iq))
    ngqr(iqr) = ngq(iq)
  end do

  ! Change file extension and write out k an q points
  call genfilname(dotext='_SCI.OUT', setfilext=.true.)
  !write(*,*)
  !write(*,*) "writing k-points to file"
  !write(*,*) "filext=",trim(filext)
  if(mpiglobal%rank == 0) then
    call writekpts
  end if
  if(fra) then 
    if(fti) then 
      call genfilname(iqmt=iqmt, auxtype="m", dotext='_SCI.OUT', setfilext=.true.)
    else
      call genfilname(iqmt=iqmt, auxtype="mqmt", dotext='_SCI.OUT', setfilext=.true.)
    end if
  else
    call genfilname(iqmt=iqmt, dotext='_SCI.OUT', setfilext=.true.)
  end if
  !write(*,*)
  !write(*,*) "writing q-points to file"
  !write(*,*) "filext=",trim(filext)
  if(mpiglobal%rank == 0) then
    call writeqpts
  end if

  ! Allocate local arrays for screened coulomb interaction and
  ! W(G,G',qr)
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  scieffg(:, :, :) = zzero
  ! Phases for transformations form reduced q points to non reduced ones.
  allocate(phf(ngqmax, ngqmax))

  ! Parallelize over reduced q-point set
  call genparidxran('q', nqptr)

  ! Set file to read RPA screening from 
  eps0dirname = 'EPS0'
  if(fra) then 
    if(.not. fti) then 
      if(iqmt == 1) then 
        call genfilname(iqmt=iqmtgamma, fileext=fileext_scr_read)
        call genfilname(iqmt=iqmt, fileext=fileext_ematrad_write)
      else
        call genfilname(iqmt=iqmt, auxtype='mqmt', fileext=fileext_scr_read)
        call genfilname(iqmt=iqmt, auxtype='mqmt', fileext=fileext_ematrad_write)
      end if
    else
      if(all(vqoff == 0.0d0)) then 
        call genfilname(iqmt=iqmtgamma, fileext=fileext_scr_read)
        call genfilname(iqmt=iqmt, fileext=fileext_ematrad_write)
      else
        call genfilname(iqmt=iqmt, auxtype='m', fileext=fileext_scr_read)
        call genfilname(iqmt=iqmt, auxtype='m', fileext=fileext_ematrad_write)
      end if
    end if
  else
    call genfilname(iqmt=iqmtgamma, fileext=fileext_scr_read)
    call genfilname(iqmt=iqmt, fileext=fileext_ematrad_write)
  end if

  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_scrcoulint):&
      & Calculating W(G1,G2,qr) fourier coefficients")')
    call timesec(tscc0)
  end if

  !write(*,*) "W Fourier coefficients"
  do iqr = qpari, qparf ! Reduced q

    !write(*,*) "iqr=", iqr
    ! Locate reduced q-point in non-reduced set
    iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))
    !write(*,*) "iqrnr=", iqrnr

    ! Get number of G+q vectors for current q
    numgq = ngq(iqrnr)

    ! Calculate effective screened coulomb interaction
    ! by inverting the symmetrized RPA dielectric matrix for a given q and
    ! 0 frequency and then multiplying
    ! it with v^{1/2} from both sides.
    filext = fileext_scr_read
    !write(*,*) "reading screening form =", trim(filext)
    call genscclieff(iqr, iqrnr, ngqmax, numgq, scieffg(:,:,iqr))

    ! Generate radial integrals for matrix elements of plane wave
    ! and save them to disk.
    filext = fileext_ematrad_write
    !write(*,*) "writing emat to =", trim(filext)
    call putematrad(iqr, iqrnr)

  end do

  ! Set file extesion for later read EMATRAD in getematrad
  ! (some ranks may not participate in the qr loop above)
  filext = fileext_ematrad_write

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(set=nqptr, rlen=ngqmax*ngqmax,&
    & zbuf=scieffg, inplace=.true., comm=mpiglobal)

  if(mpiglobal%rank == 0) then
    call timesec(tscc1)
    write(unitout, '("  Timing (in seconds):", f12.3)') tscc1 - tscc0
  end if

  !-------------------------------!
  ! CONSTRUCT W MATRIX ELEMENTS   !
  !-------------------------------!

  ! Normalization factor 1/V and per k point
  pref=1.0d0/(omega*dble(nk_bse))

  ! Allocate arrays used in ematqk (do not change in the following loop)
  call ematqalloc

  ! Work arrays (allocate for maximal size over all participating k points)
  allocate(sccli_t2(nu_bse_max, no_bse_max, nu_bse_max, no_bse_max))
  allocate(sccli_t1(no_bse_max**2, nu_bse_max**2))
  allocate(sccli(nou_bse_max, nou_bse_max))

  if(.not. fra) then 
    allocate(muu(nu_bse_max, nu_bse_max, ngqmax))
    allocate(moo(no_bse_max, no_bse_max, ngqmax))
  else
    allocate(mou(no_bse_max, nu_bse_max, ngqmax))
    allocate(muo(nu_bse_max, no_bse_max, ngqmax))
  end if


  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(b_scrcoulint): W matrix elements")')
    call timesec(tscc0)
  end if

  ! Distributed loop over combinations of non-reduced k-point
  ! that contribute to the desired energy range (or bands).
  call genparidxran('p', nkkp_bse)

  kkploop: do ikkp = ppari, pparf

    ! Get individual k-point indices from combined ik jk index.
    !   ik runs from 1 to nk_bse, jk from ik to nk_bse
    call kkpmap(ikkp, nk_bse, ik, jk)
    !! Note: If the screened coulomb interaction matrix 
    !! has the elements W_{i,j} and the indices enumerate the
    !! states according to
    !! i = {u1o1k1, u2o1k1, ..., uMo1k1,
    !!      u1o2k1, ..., uMoNk1, u1o1k2, ..., uMoNkO} -> {1,...,M*N*O}
    !! then because of W_{j,i} = W^*_{i,j} only kj = ki,..,kN is 
    !! needed (the diagonal blocks where k'=k will be fully computed 
    !! but only the upper triangle will be needed in the end).
    !! (The RA part in the standard basis is symmetric instead of 
    !!  hermitian, but one still just needs the upper triangle)

    ! Get total k point indices
    !write(*,*)
    !write(*,'(a, i4)') "ikkp =", ikkp

    ! Get k point indices  
    iknr = kmap_bse_rg(ik)
    jknr = kmap_bse_rg(jk) 

    !write(*,'(a, i4)') "iknr =", iknr
    !write(*,'(a, i4)') "jknr =", jknr

    ! Get corresponding k' (R case k'=k+qmt, A case k'=k-qmt, A^ti case k'=-(k+qmt))
    if(fra .and. .not. fti) then
      ! Get index of ik+qmt
      ikpnr = k_kqmtp%ik2ikqmt(iknr)
      ! Get index of jk-qmt
      jkpnr = k_kqmtm%ik2ikqmt(jknr)
      !write(*,'(a, i4)') "ik+qmt: ikpnr =", ikpnr
      !write(*,'(a, i4)') "jk-qmt: jkpnt =", jkpnr
    else
      ! Get index of ik+qmt
      ikpnr = k_kqmtp%ik2ikqmt(iknr)
      ! Get index of jk+qmt
      jkpnr = k_kqmtp%ik2ikqmt(jknr)
      !write(*,'(a, i4)') "ik+qmt: ikpnr =", ikpnr
      !write(*,'(a, i4)') "jk+qmt: jkpnt =", jkpnr
    end if
    !write(*,*)

    ! Get corresponding q-point
    ! (RR case: q = jk-ik, RA case: q = jk-ik-qmt, RA^ti case: q = -(jk+qmt)-ik)
    if(fra) then 
      if(fti) then 
        iq = p_pqmtp%ikikp2ip_nr(iknr, jkpnr)
        iqr = p_pqmtp%pset%ik2ikp(iq)
        ! Get corresponding vector in lattice coordinated
        vqr(:) = p_pqmtp%pset%vkl(:, iqr)
        vq(:) = p_pqmtp%pset%vklnr(:, iq)
        !write(*,'(a, i4, 3E10.3)') "-jk-ik-qmt: iq,vq =", iq, vq
        !write(*,'(a, i4, 3E10.3)') "-jk-ik-qmt: iqr,vqr =", iqr, vqr
      else
        iq = q_qmtm%ikikp2iq_nr(iknr, jkpnr)
        iqr = q_qmtm%qset%ik2ikp(iq)
        vqr(:) = q_qmtm%qset%vkl(:, iqr)
        vq(:) = q_qmtm%qset%vklnr(:, iq)
        !write(*,'(a, i4, 3E10.3)') "jk-ik-qmt: iq,vq =", iq, vq
        !write(*,'(a, i4, 3E10.3)') "jk-ik-qmt: iqr,vqr =", iqr, vqr
      end if
    else
      iq = q_q%ikikp2iq_nr(iknr, jknr)
      iqr = q_q%qset%ik2ikp(iq)
      vqr(:) = q_q%qset%vkl(:, iqr)
      vq(:) = q_q%qset%vklnr(:, iq)
      !write(*,'(a, i4, 3E10.3)') "jk-ik: iq,vq =", iq, vq
      !write(*,'(a, i4, 3E10.3)') "jk-ik: iqr,vqr =", iqr, vqr
    end if
    !write(*,*)

    !write(*,*) "check: vq =", vq
    !write(*,*) "check: vql =", vql(1:3,iq)

    ! Check if iq is Gamma point (mod_qpoint::vqc(:,iq) has length < 1d-12)
    tq0 = tqgamma(iq)
    !write(*,*) "q is gamma?=", tq0, (norm2(vq)<1.0d-6)
    !write(*,*)

    ! Local field effects size (Number of G+q vectors)
    numgq = ngq(iq)
    !write(*,*) "numgq=", numgq

    if(allocated(igqmap)) deallocate(igqmap)
    allocate(igqmap(numgq))
    if(allocated(wfc)) deallocate(wfc)
    allocate(wfc(numgq, numgq))

    ! Get radial integrals for q_r (previously calculated for reduced q set)
    call getematrad(iqr, iq)

    if(input%xs%reduceq) then 
      !! Find results  
      !! for a non reduced q-point with the help of the result 
      !! for the corresponding reduced q-point using symmetry operations.
      !! (Radial emat integrals and screened coulomb potential Fourier coefficients)
      !!<--
      !! RR & RA
      ! Find symmetry operations that map the reduced q-point to the non reduced one
      call findsymeqiv(input%xs%bse%fbzq, vq, vqr, nsc, sc, ivgsc)
      ! Find a crystal symmetry operation that rotates the G+q-vectors onto G'+q_r-vectors
      ! and generate a Map G' --> G
      call findgqmap(iq, iqr, nsc, sc, ivgsc, numgq, jsym, jsymi, ivgsym, igqmap)
      ! Rotate radial integrals calculated for the reduced q to get those for non-reduced q
      call rotematrad(numgq, igqmap)
      ! Generate phase factor for dielectric matrix due to non-primitive
      ! translations
      call genphasedm(iq, jsym, ngqmax, numgq, phf, tphf)
      ! W(G,G',q) <-- W(\tilde{G},\tilde{G}',qr)
      wfc(:,:) = phf(:numgq, :numgq) * scieffg(igqmap, igqmap, iqr)
    else
      wfc(:,:) = scieffg(1:numgq, 1:numgq, iqr)
    end if

    if(fwp) then 
      if(fti) then 
        call genfilname(iq=iq, auxtype='m', dotext='', fileext=wfc_write)
      else
        call genfilname(iq=iq, dotext='', fileext=wfc_write)
      end if
      wfc_write = 'Wfc'//trim(adjustl(wfc_write))
      call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(wfc),&
        & immat=aimag(wfc), dirname='Wfc')
    end if

    ! Get ik & jk dependent band ranges for 
    ! plane wave matrix calculation (is nstlbse was used in input all k points
    ! contribute the same number of transitions)
    inu = koulims(2,iknr) - koulims(1,iknr) + 1
    ino = koulims(4,iknr) - koulims(3,iknr) + 1
    jnu = koulims(2,jknr) - koulims(1,jknr) + 1
    jno = koulims(4,jknr) - koulims(3,jknr) + 1
    noo = ino*jno
    nuu = inu*jnu
    nou = ino*jnu
    nuo = inu*jno
    ! Number of transitions at ik & jk
    inou = kousize(iknr)
    jnou = kousize(jknr)

    ! W^RR
    if(.not. fra) then 
      !-------------------------------!
      ! Resonant-Resonant Part        !
      !-------------------------------!

      allocate(cmoo(noo, numgq), cmuu(nuu, numgq))

      ! Calculate M_{io jo ik}(G, q) = <io ik|e^{-i(q+G)r}|jo jk>
      ! and       M_{iu ju ikp}(G, q) = <iu ikp|e^{-i(q+G)r}|ju jkp>
      ! where ikp=ik+qmt, jkp=jk+qmt, q=jk-ik
      call getpwesrr(moo(1:ino,1:jno,1:numgq), muu(1:inu,1:jnu,1:numgq))

      ! Combine indices for matrix elements of plane wave.

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(io,jo,j1,iu,ju,j2)
      !$OMP DO COLLAPSE(2)
      do jo = 1, jno   ! jo
        do io = 1, ino ! io
          j1 = io + (jo-1)*ino ! iojo
          ! cmoo_j = M_o1o1, M_o2o1, ..., M_oNo1, M_o1o2, ..., M_oNoM
          cmoo(j1, :) = moo(io, jo, 1:numgq)
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP DO COLLAPSE(2)
      do ju = 1, jnu   ! ju
        do iu = 1, inu ! iu
          j2 = iu + (ju-1)*inu ! iuju
          ! cmuu_j = M_u1u1, M_u2u1, ..., M_uNu1, M_u1u2, ..., M_uNuM
          cmuu(j2, :) = muu(iu, ju, 1:numgq)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! M_{iojo} -> M^*_{iojo}
      cmoo = conjg(cmoo)

      ! Allocate helper array
      allocate(zm(noo,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclit_{j1, j2}(q)
      ! zm = cmoo * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmoo_{j1,G} wfc_{G,G'}
      !   i.e. zm_{io_j1,jo_j1}(G',q) = \Sum_{G} M^*_{io_j1,jo_j1}(G,q) W(G,G', q)
      call zgemm('n', 'n', noo, numgq, numgq, zone, cmoo, noo,&
        & wfc, numgq, zzero, zm, noo)
      ! scclit = pref * zm * cmuu^T
      !   i.e. scclit_{j1, j2} = \Sum_{G'} zm_{j1,G'} (cmuu^T)_{G',j2}
      !   i.e. scclit_{io_j1 jo_j1, iu_j2 ju_j2} = 
      !          \Sum{G,G'} M^*_{io_j1,jo_j1}(G,q) W(G,G',q) M_{iu_j2 ju_j2}(G',q)
      call zgemm('n', 't', noo, nuu, numgq, pref, zm, noo,&
        & cmuu, nuu, zzero, sccli_t1(1:noo,1:nuu), noo)

      deallocate(zm)        
      deallocate(cmoo, cmuu)

      ! Save only the selected transitions
      jaoff = sum(kousize(1:jknr-1))
      iaoff = sum(kousize(1:iknr-1))

      ! Map back to individual band indices
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iu,ju,io,jo,j1,j2,ia,ja)
      !$OMP DO COLLAPSE(4)
      do ju = 1, jnu    ! ju
        do iu = 1, inu  ! iu
          do jo = 1, jno   ! jo
            do io = 1, ino ! io
              j2 = iu + (ju-1)*inu
              j1 = io + (jo-1)*ino
              ! scclit_{io_j1 jo_j1, iu_j2 ju_j2} -> sccli_{iu_j2 io_j1, ju_j2 jo_j1}
              sccli_t2(iu, io, ju, jo) = sccli_t1(j1, j2)
            end do
          end do
        end do
      end do
      !$OMP END DO
      ! W^RR matrix element arrays for one jk-ik=q
      !$OMP DO COLLAPSE(2)
      do ja = 1, jnou
        do ia = 1, inou
          ju = smap_rel(1,ja+jaoff)
          jo = smap_rel(2,ja+jaoff)
          iu = smap_rel(1,ia+iaoff)
          io = smap_rel(2,ia+iaoff)
          ! scclit_{io_j1 jo_j1, iu_j2 ju_j2} -> sccliab_{iu_a io_a, ju_b jo_b}
          sccli(ia, ja) = sccli_t2(iu, io, ju, jo)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      if(fwp) then 
        call writecmplxparts('Wrr', remat=dble(sccli(1:inou,1:jnou)),&
          & immat=aimag(sccli(1:inou,1:jnou)), ik1=iknr, ik2=jknr, dirname='Wrr')
      end if
      ! Parallel write
      call b_putbsemat(scclifname, 77, ikkp, iqmt, sccli)

    ! W^{RA} or W^{RA,ti}
    else
      !-------------------------------!
      ! Resonant-Anti-Resonant Part   !
      !-------------------------------!

      allocate(cmou(nou, numgq), cmuo(nuo, numgq))

      ! For non-time-inverted anti-resonant basis:
      ! Calculate M_{io ju ik}(G, q) = <io ik|e^{-i(q+G)r}|ju jkp>
      ! and       M_{iu jo ikp}(G, q) = <iu ikp|e^{-i(q+G)r}|jo jk>
      ! where ikp=ik+qmt, jkp=jk-qmt, q=jk-ik-qmt
      !
      ! For time-inverted anti-resonant basis:
      ! Calculate \tilde{M}_{io ju ik}(G, q) = <io ik|e^{-i(q+G)r}|(ju jkp)^*>
      ! and       \tilde{M}_{iu jo ikp}(G, q) = <iu ikp|e^{-i(q+G)r}|(jo jk)^*>
      ! where ikp=ik+qmt, jkp=jk+qmt, q=-(jk+qmt)-ik
      !
      call getpwesra(mou(1:ino,1:jnu,1:numgq), muo(1:inu,1:jno,1:numgq))

      ! Combine indices for matrix elements of plane wave.
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,ju,j1)
      do ju = 1, jnu   ! ju
        do io = 1, ino ! io
          j1 = io+(ju-1)*ino
          ! cmou_j = M_o1u1, M_o2u1, ..., M_oNu1, M_o1u2, ..., M_oNuM
          cmou(j1, :) = mou(io, ju, 1:numgq)
        end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,ju,j1)
      do jo = 1, jno   ! jo
        do iu = 1, inu ! iu
          j2 = iu+(jo-1)*inu
          ! cmuo_j = M_u1o1, M_u2o1, ..., M_uMo1, M_u1o2, ..., M_uMoN
          cmuo(j2, :) = muo(iu, jo, 1:numgq)
        end do
      end do
      !$OMP END PARALLEL DO

      ! M_{ioju} -> M^*_{ioju}
      cmou = conjg(cmou)

      ! Allocate helper array
      allocate(zm(nou,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclict_{j1, j2}(q)
      ! zm = cmou * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmou_{j1,G} wfc_{G,G'}
      !   i.e. zm_{io_j1,ju_j1}(G',q) =
      !          \Sum_{G} M^*_{io_j1,ju_j1,ik}(G,q-qmt) W(G,G',q-qmt)
      call zgemm('n', 'n', nou, numgq, numgq, zone, cmou, nou,&
        & wfc, numgq, zzero, zm, nou)
      ! scclit = pref * zm * cmuo^T
      !   i.e. scclict(j1, j2) = \Sum_{G'} zm_{j1,G'} (cmuo^T)_{G',j2}
      !   i.e. scclict_{io_j1 ju_j1, iu_j2 jo_j2} =
      !   \Sum{G,G'} 
      !   M^*_{io_j1,ju_j1,ik}(G,q-qmt) W(G, G',q-qmt) M_{iu_j2,jo_j2,ik+qmt}(G',q-qmt)
      call zgemm('n', 't', nou, nuo, numgq, pref, zm,&
        & nou, cmuo, nuo, zzero, sccli_t1(1:nou,1:nuo), nou)
      deallocate(zm)        
      deallocate(cmou, cmuo)

      ! Map back to individual band indices
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(4),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,jo,iu,ju,j1,j2)
      do ju = 1, jnu   ! ju
        do io = 1, ino ! io
          do jo = 1, jno   ! jo
            do iu = 1, inu ! iu
              j1 = io+(ju-1)*ino !ioju
              j2 = iu+(jo-1)*inu ! iujo
              ! scclict_{io_j1 ju_j1, iu_j2 jo_j2}(ik,jk)(qmt)
              ! -> scclic(iu_j2, io_j1, ju_j1, jo_j2)(ik,jk)(qmt)
              sccli_t2(iu, io, ju, jo) = sccli_t1(j1, j2)
            end do
          end do
        end do
      end do 
      !$OMP END PARALLEL DO

      ! W matrix elements

      ! Save only the selected transitions
      jaoff = sum(kousize(1:jknr-1))
      iaoff = sum(kousize(1:iknr-1))
      !$OMP PARALLEL DO &
      !$OMP& COLLAPSE(2),&
      !$OMP& DEFAULT(SHARED), PRIVATE(io,jo,iu,ju,ia,ja)
      do ja = 1, jnou
        do ia = 1, inou
          ju = smap_rel(1,ja+jaoff)
          jo = smap_rel(2,ja+jaoff)
          iu = smap_rel(1,ia+iaoff)
          io = smap_rel(2,ia+iaoff)
          ! scclit_{io_j1 jo_j1, iu_j2 ju_j2} -> sccliab_{iu_a io_a, ju_b jo_b}
          sccli(ia, ja) = sccli_t2(iu, io, ju, jo)
        end do
      end do
      !$OMP END PARALLEL DO

      if(fwp) then 
        if( fti ) then 
          call writecmplxparts('Wra_ti', dble(sccli(1:inou,1:jnou)),&
            & aimag(sccli(1:inou,1:jnou)), ik1=iknr, ik2=jknr, dirname='Wra_ti')
        else 
          call writecmplxparts('Wra', dble(sccli(1:inou,1:jnou)),&
            & aimag(sccli(1:inou,1:jnou)), ik1=iknr, ik2=jknr, dirname='Wra')
        end if
      end if

      ! Parallel write
      call b_putbsemat(scclifname, 78, ikkp, iqmt, sccli)

    end if

    ! Deallocate G+q dependent work arrays
    deallocate(igqmap)
    deallocate(wfc)

    !if(rank == 0) then
    !  write(6, '(a,"Scrcoulint progess:", f10.3)', advance="no")&
    !    & achar( 13), 100.0d0*dble(ikkp-ppari+1)/dble(pparf-ppari+1)
    !  flush(6)
    !end if

  ! End loop over(k,kp)-pairs
  end do kkploop

  !if(rank == 0) then
  !  write(*,*)
  !end if

  if(mpiglobal%rank == 0) then
    call timesec(tscc1)
    write(unitout, '("  Timing (in seconds):", f12.3)') tscc1 - tscc0
  end if

  ! Deallocate helper array 
  deallocate(sccli)
  deallocate(sccli_t1)
  deallocate(sccli_t2)
  if(.not. fra) then 
    deallocate(muu)
    deallocate(moo)
  else
    deallocate(mou)
    deallocate(muo)
  end if

  call xsgrids_finalize()

  call barrier

  call findgntn0_clear

  contains

    subroutine getpwesrr(moo, muu)
      complex(8), intent(out) :: moo(:,:,:), muu(:,:,:)

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save

      !write(*,*)
      !write(*,*) "getpwesrr:"
      !write(*,*) "  Moo"

      fileext0_save = filext0
      fileext_save = filext

      !-----------------------------------------------------------!
      ! Calculate M_{io jo ik}(G, q) = <io ik|e^{-i(q+G)r}|jo jk> !
      !-----------------------------------------------------------!
      ! Bands
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=jno
      ematbc%il2=koulims(3,jknr)
      ematbc%iu2=koulims(4,jknr)

      ! Set EVECFV_QMT001.OUT as bra state file
      usefilext0 = .true.
      iqmt0 = iqmtgamma
      call genfilname(iqmt=iqmt0, setfilext=.true.)
      filext0 = filext
      !write(*,*) "filext0 =", trim(filext0)

      ! Set EVECFV_QMT001.OUT as ket state file
      iqmt1 = iqmtgamma
      call genfilname(iqmt=iqmt1, setfilext=.true.)
      !write(*,*) "filext =", trim(filext)

      ! Set vkl and vkl0 to k-grid
      if(iqmt /= 1) then
        call init1offs(k_kqmtp%kset%vkloff)
        call xssave0
      end if
      ! Set up ikmapikq to link (ik,iq) to jk
      ! (all other q and G+q dependent variables need not be changed)
      ikmapikq(1:nkpt, 1:nqpt) = q_q%ikiq2ikp_nr(1:nkpt, 1:nqpt)

      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iq, iknr, moo, ematbc)
      !-----------------------------------------------------------!
      if(fwp) then
        do igq=1,numgq
          call genfilname(iqmt=iq, iq=igq, dotext='', fileext=wfc_write)
          wfc_write='Moo'//trim(adjustl(wfc_write))
          call genfilname(iqmt=iq, dotext='', fileext=dirname)
          dirname='Moo'//trim(adjustl(dirname))
          call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(moo(:,:,igq)),&
            & immat=aimag(moo(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
        end do
      end if

      !write(*,*) "  Muu"
      !------------------------------------------------------------------!
      ! Calculate M_{iu ju ik+qmt}(G, q) = <iu ikp|e^{-i(q+G)r}|ju jkp>  !
      !------------------------------------------------------------------!
      ! Bands
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=jnu
      ematbc%il2=koulims(1,jknr)
      ematbc%iu2=koulims(2,jknr)

      ! Set EVECFV_QMTXYZ.OUT as bra state file
      usefilext0 = .true.
      iqmt0 = iqmt
      call genfilname(iqmt=iqmt0, setfilext=.true.)
      filext0 = filext
      !write(*,*) "filext0 =", trim(filext0)

      ! Set EVECFV_QMTXYZ.OUT as ket state file
      iqmt1 = iqmt
      call genfilname(iqmt=iqmt1, setfilext=.true.)
      !write(*,*) "filext =", trim(filext)

      ! Set vkl and vkl0 to k+qmt-grid
      if(iqmt /= 1) then 
        call init1offs(k_kqmtp%kqmtset%vkloff)
        call xssave0
      end if
      ! Set up ikmapikq to link (ikp,iq) to (jkp)
      ! (all other q and G+q dependent variables need not be changed)
      ikmapikq(1:nkpt, 1:nqpt) = qmtp_qmtp%ikiq2ikp_nr(1:nkpt, 1:nqpt)

      ! Calculate M_{o1o2,G} at fixed (k, q)
      call b_ematqk(iq, ikpnr, muu, ematbc)
      !------------------------------------------------------------------!

      if(fwp) then
        do igq=1,numgq
          call genfilname(iqmt=iq, iq=igq, dotext='', fileext=wfc_write)
          wfc_write='Muu'//trim(adjustl(wfc_write))
          call genfilname(iqmt=iq, dotext='', fileext=dirname)
          dirname='Muu'//trim(adjustl(dirname))
          call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(muu(:,:,igq)),&
            & immat=aimag(muu(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
        end do
      end if

      filext0 = fileext0_save
      filext = fileext_save

    end subroutine getpwesrr

    subroutine getpwesra(mou, muo)
      complex(8), intent(out) :: mou(:,:,:), muo(:,:,:)

      character(256) :: fileext0_save, fileext_save
      type(bcbs) :: ematbc

      !write(*,*)
      !write(*,*) "getpwesra:"

      fileext0_save = filext0
      fileext_save = filext

      ! Standard anti-resonant basis
      if(.not. fti) then 

        !write(*,*) " Std"
        !write(*,*) "  Mou"
        !------------------------------------------------------------!
        ! Calculate M_{io ju ik}(G, q) = <io ik|e^{-i(q+G)r}|ju jkp> !
        ! where jkp=jk-qmt, q=jk-ik-qmt                              !
        !------------------------------------------------------------!
        
        ! Bands
        ematbc%n1=ino
        ematbc%il1=koulims(3,iknr)
        ematbc%iu1=koulims(4,iknr)
        ematbc%n2=jnu
        ematbc%il2=koulims(1,jknr)
        ematbc%iu2=koulims(2,jknr)

        ! Set EVECFV_QMT001.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true.)
        filext0 = filext
        !write(*,*) "filext0 =", trim(filext0)

        if(iqmt /= 1) then 
          ! Set EVECFV_QMTXYZ_mqmt.OUT as ket state file
          iqmt1 = iqmt
          call genfilname(iqmt=iqmt1, auxtype="mqmt", setfilext=.true.)
          !write(*,*) "filext =", trim(filext)
        else
          ! Set EVECFV_QMT001.OUT as ket state file
          iqmt1 = iqmtgamma
          call genfilname(iqmt=iqmt1, setfilext=.true.)
          !write(*,*) "filext =", trim(filext)
        end if

        if(iqmt /= 1) then 
          ! Set vkl0 to k-grid
          call init1offs(k_kqmtm%kset%vkloff)
          call xssave0
          ! Set vkl to k-qmt-grid
          call init1offs(k_kqmtm%kqmtset%vkloff)
        end if

        ! Set up ikmapikq to link (ik,iq) to jkp 
        ikmapikq(1:nkpt, 1:nqpt) = q_qmtm%ikiq2ikp_nr(1:nkpt, 1:nqpt)

        ! Calculate M_{ou,G} at fixed (k, q)
        call b_ematqk(iq, iknr, mou, ematbc)
        !------------------------------------------------------------!

        if(fwp) then
          do igq=1,numgq
            call genfilname(iqmt=iq, iq=igq, dotext='', fileext=wfc_write)
            wfc_write='Mou'//trim(adjustl(wfc_write))
            call genfilname(iqmt=iq, dotext='', fileext=dirname)
            dirname='Mou'//trim(adjustl(dirname))
            call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(mou(:,:,igq)),&
              & immat=aimag(mou(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
          end do
        end if

        !write(*,*) "  Mou"
        !-------------------------------------------------------------!
        ! Calculate M_{iu jo ikp}(G, q) = <iu ikp|e^{-i(q+G)r}|jo jk> !
        ! where ikp=ik+qmt, q=jk-ik-qmt
        !-------------------------------------------------------------!

        ! Bands
        ematbc%n1=inu
        ematbc%il1=koulims(1,iknr)
        ematbc%iu1=koulims(2,iknr)
        ematbc%n2=jno
        ematbc%il2=koulims(3,jknr)
        ematbc%iu2=koulims(4,jknr)

        ! Set EVECFV_QMTXYZ.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt0, setfilext=.true.)
        filext0 = filext
        !write(*,*) "filext0 =", trim(filext0)

        ! Set EVECFV_QMT001.OUT as ket state file
        iqmt1 = iqmtgamma
        call genfilname(iqmt=iqmt1, setfilext=.true.)
        !write(*,*) "filext =", trim(filext)

        if(iqmt /= 1) then 
          ! Set vkl0 to k+qmt-grid
          call init1offs(k_kqmtp%kqmtset%vkloff)
          call xssave0
          ! Set vkl to k-grid
          call init1offs(k_kqmtp%kset%vkloff)
        end if

        ! Set up ikmapikq to link (ikp,iq) to jk
        ikmapikq(1:nkpt, 1:nqpt) = qmtp_q%ikiq2ikp_nr(1:nkpt, 1:nqpt)

        ! Calculate M_{uo,G} at fixed (k, q)
        call b_ematqk(iq, ikpnr, muo, ematbc)
        !-------------------------------------------------------------!

        if(fwp) then
          do igq=1,numgq
            call genfilname(iqmt=iq, iq=igq, dotext='', fileext=wfc_write)
            wfc_write='Muo'//trim(adjustl(wfc_write))
            call genfilname(iqmt=iq, dotext='', fileext=dirname)
            dirname='Muo'//trim(adjustl(dirname))
            call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(muo(:,:,igq)),&
              & immat=aimag(muo(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
          end do
        end if

      else

        !write(*,*) " TI"
        !write(*,*) "  Mou"
        !------------------------------------------------------------------------!
        ! Calculate \tilde{M}_{io ju ik}(G, q) = <io ik|e^{-i(q+G)r}|(ju jkp)^*> !
        ! where jkp=jk+qmt, q=-(jk+qmt)-ik                                       !
        !------------------------------------------------------------------------!
        
        ! Bands
        ematbc%n1=ino
        ematbc%il1=koulims(3,iknr)
        ematbc%iu1=koulims(4,iknr)
        ematbc%n2=jnu
        ematbc%il2=koulims(1,jknr)
        ematbc%iu2=koulims(2,jknr)

        ! Set EVECFV_QMT001.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true.)
        filext0 = filext
        !write(*,*) "filext0 =", trim(filext0)

        iqmt1 = iqmt
        ! Set EVECFV_QMTYZ.OUT as ket state file
        call genfilname(iqmt=iqmt1, setfilext=.true.)
        !write(*,*) "filext =", trim(filext)

        ! Set vkl0 to k-grid
        if(iqmt /= 1) then 
          call init1offs(k_kqmtp%kset%vkloff)
          call xssave0
          ! Set vkl to k+qmt-grid
          call init1offs(k_kqmtp%kqmtset%vkloff)
        end if

        ! Set up non reduced ikmapikq to link (ik,iq) to jkp
        ikmapikq(1:nkpt, 1:nqpt) = p_pqmtp%ikip2ikp_nr(1:nkpt, 1:nqpt)

        ! Calculate M_{ou,G} at fixed (k, q)
        call b_ematqk2(iq, iknr, mou, ematbc)
        !------------------------------------------------------------!

        if(fwp) then
          do igq=1,numgq
            call genfilname(iqmt=iq, iq=igq, dotext='', fileext=wfc_write)
            wfc_write='Mou_ti'//trim(adjustl(wfc_write))
            call genfilname(iqmt=iq, dotext='', fileext=dirname)
            dirname='Mou_ti'//trim(adjustl(dirname))
            call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(mou(:,:,igq)),&
              & immat=aimag(mou(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
          end do
        end if

        !write(*,*) "  Muo"
        !-------------------------------------------------------------------------!
        ! Calculate \tilde{M}_{iu jo ikp}(G, q) = <iu ikp|e^{-i(q+G)r}|(jo jk)^*> !
        ! where ikp=ik+qmt, q=-jk-(ik+qmt)                                        !
        !-------------------------------------------------------------------------!

        ! Bands
        ematbc%n1=inu
        ematbc%il1=koulims(1,iknr)
        ematbc%iu1=koulims(2,iknr)
        ematbc%n2=jno
        ematbc%il2=koulims(3,jknr)
        ematbc%iu2=koulims(4,jknr)

        ! Set EVECFV_QMTXYZ.OUT as bra state file
        usefilext0 = .true.
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt0, setfilext=.true.)
        filext0 = filext
        !write(*,*) "filext0 =", trim(filext0)

        iqmt1=iqmtgamma
        ! Set EVECFV_QMT001.OUT as ket state file
        call genfilname(iqmt=iqmt1, setfilext=.true.)
        !write(*,*) "filext =", trim(filext)

        if(iqmt /= 1) then
          ! Set vkl0 to k+qmt-grid
          call init1offs(k_kqmtp%kqmtset%vkloff)
          call xssave0
          ! Set vkl to k-grid
          call init1offs(k_kqmtp%kset%vkloff)
        end if

        ! Set up ikmapikq to link ikp,iq to jk
        ikmapikq(1:nkpt, 1:nqpt) = pqmtp_p%ikip2ikp_nr(1:nkpt, 1:nqpt)

        ! Calculate M_{uo,G} at fixed (k, q)
        call b_ematqk2(iq, ikpnr, muo, ematbc)
        !-------------------------------------------------------------!
        if(fwp) then
          do igq=1,numgq
            call genfilname(iqmt=iq, iq=igq, dotext='', fileext=wfc_write)
            wfc_write='Muo_ti'//trim(adjustl(wfc_write))
            call genfilname(iqmt=iq, dotext='', fileext=dirname)
            dirname='Muo_ti'//trim(adjustl(dirname))
            call writecmplxparts(trim(adjustl(wfc_write)), remat=dble(muo(:,:,igq)),&
              & immat=aimag(muo(:,:,igq)), ik1=iknr, ik2=jknr, dirname=dirname)
          end do
        end if

      end if

      filext0 = fileext0_save
      filext = fileext_save
    end subroutine getpwesra

end subroutine b_scrcoulint
!EOC

