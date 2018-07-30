! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: scrcoulint
! !INTERFACE:
subroutine scrcoulint(iqmt, fra)
! !USES:
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi
  use mod_constants, only: zzero, zone, fourpi
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: iqmap, vql, vqc, nqpt, ivq, wqpt
  use mod_lattice, only: omega
  use mod_symmetry, only: maxsymcrys
  use modxs, only: xsgnt, unitout,&
                 & ngqmax,&
                 & nqptr, qpari, qparf, ivqr,&
                 & ngq, ppari, pparf, iqmapr,&
                 & vqlr, vqcr, wqptr, ngqr,&
                 & bcbs, ematraddir, eps0dirname,&
                 & filext0, usefilext0, iqmt0, iqmt1,&
                 & iqmtgamma,&
                 & bsedl, bsedu, bsedd, bsed
  use m_xsgauntgen
  use m_findgntn0
  use m_writevars
  use m_genfilname
  use m_getunit
  use m_ematqk
  use m_putgetbsemat
  use modbse
  use mod_xsgrids
  use mod_Gkvector, only: gkmax
! !DESCRIPTION:
!   Calculates the resonant-resonant or resonant-anit-resonant block of the
!   direct term of the Bethe-Salpeter Hamiltonian for a momentum transfer 
!   $\vec{Q}_{mt}$.
!
! !REVISION HISTORY:
!   Forked from scrcoulint.F90 and adapted for non-TDA BSE and finite Q. (Aurich)
!EOP
!BOC      

  implicit none

  ! I/O
  integer, intent(in) :: iqmt ! Index of momentum transfer Q
  logical, intent(in) :: fra  ! Construct RA coupling block

  ! Local variables
  character(*), parameter :: thisname = 'scrcoulint'

  ! ik,jk block of W matrix (final product)
  complex(8), allocatable :: sccli(:,:)

  ! Diagonal of W (needed for MB TDDFT kernels)
  complex(8), allocatable :: bsedt(:, :)

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
  integer(4) :: ikkp, iknr, jknr, ikpnr, ikmnr, jkpnr, jkmnr, ik, jk
  integer(4) :: iqrnr, iqr, iq
  real(8) :: vqr(3), vq(3)
  integer(4) :: numgq
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
  integer(4) :: j1, j2
  complex(8) :: pref, zt1
  real(8) :: t1
  ! Timing vars
  real(8) :: tscc1, tscc0

  ! Auxilliary strings
  character(256) :: syscommand, fileext_scr_read, fileext_ematrad_write

  ! External functions
  logical, external :: tqgamma

  real(8) :: vqoff(3)
  real(8), parameter :: epslat = 1.0d-8

  logical :: fsameq, fsamekp, fsamekm

  !---------------!
  !   main part   !
  !---------------!

  !write(*,*) "Hello, this is scrcoulint at rank:", rank

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

  syscommand = 'test ! -e '//trim(adjustl(ematraddir))&
    & //' && mkdir '//trim(adjustl(ematraddir))//' &> /dev/null'
  call system(trim(adjustl(syscommand)))

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

  if (input%xs%BSE%outputlevelnumber == 1) then
    write(unitout, '(a)') 'Info(' // thisname // '):&
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
  call printline(unitout, '-')
  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Inspecting occupations...'
  call flushifc(unitout)

  call setranges_modxs(iqmt)

  ! Select relevant transitions for the construction
  ! of the BSE hamiltonian
  ! Also sets nkkp_bse, nk_bse 
  call printline(unitout, '-')
  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Selecting transitions...'
  call flushifc(unitout)

  call select_transitions(iqmt, serial=.false.)

  ! Write support information to file
  if(mpiglobal%rank == 0) then
    if(fra) then
      call genfilname(basename=trim(infofbasename)//'_'//trim(scclicfbasename),&
        & iqmt=iqmt, filnam=infofname)
    else
      call genfilname(basename=trim(infofbasename)//'_'//trim(scclifbasename),&
        & iqmt=iqmt, filnam=infofname)
    end if
    call putbseinfo(infofname, iqmt)
  end if

  ! Set output file name
  if(fra) then
    call genfilname(basename=scclicfbasename, iqmt=iqmt, filnam=scclifname)
  else
    call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=scclifname)
  end if

  call printline(unitout, '-')
  write(unitout, '("Info(",a,"): Size of file ",a," will be about ", f12.6, " GB" )')&
    & trim(thisname), trim(scclifbasename),&
    & int(nou_bse_max,8)**2*int(nkkp_bse,8)*16.0d0/1024.0d0**3
  call flushifc(unitout)

  !------------------------------------!
  ! GENERATE FOURIER COEFFICIENTS OF W !     
  ! (and radial integrals for emat)    !
  !------------------------------------!

  call xsgrids_init(totalqlmt(1:3, iqmt), gkmax)

  !--------------------------------------------------------------------------------!
  ! Setup q points
  !--------------------------------------------------------------------------------!
  ! Get q grid offset
  if(fra) then 
    ! This is -2*vkloff for all qmt
    vqoff = pqmt%pset%vkloff
  else
    ! This is zero
    vqoff = q%qset%vkloff
  end if

  ! Check whether q and p grid are identical
  if(all(abs(vqoff) < epslat)) then 
    write(unitout, '("Info(scrcoulint): Using q-grid")')
    fsameq=.true.
  else
    write(unitout, '("Info(scrcoulint): Using shifted q-grid (p-grid)")')
    fsameq=.false.
  end if

  ! Make reduced q-grid with possible offset
  ! This sets up also G+q quantities and the square root of the Coulomb potential
  ! but the second call below will override these
  call init1offs(k_kqmtp%kset%vkloff)
  vql = 0.0d0
  call init2offs(vqoff, input%xs%reduceq)
  ! Copy results q-ponit results form mod_qpoint into modxs variables
  nqptr = nqpt
  ivqr = ivq
  iqmapr = iqmap
  vqlr = vql
  vqcr = vqc
  wqptr = wqpt
  ! Make non-reduced q-grid with possible offset (written into mod_qpoint variables)
  ! This sets up also G+q quantities and the square root of the Coulomb potential
  call init2offs(vqoff, .false.)
  write(unitout, '(a, i6)') 'Info(' // thisname // '):&
    & Number of reduced q-points: ', nqptr
  write(unitout, '(a, i6)') 'Info(' // thisname // '):&
    & Number of non-reduced q-points: ', nqpt
  write(unitout,*)
  call flushifc(unitout)
  ! Make also ngqr
  if(allocated(ngqr)) deallocate(ngqr)
  allocate(ngqr(nqptr))
  ngqr=0
  do iq= 1, nqpt
    iqr = iqmapr(ivq(1,iq), ivq(2,iq), ivq(3,iq))
    ngqr(iqr) = ngq(iq)
  end do
  !--------------------------------------------------------------------------------!

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
    if (iqmt .ne. 1) then
      write(unitout, '("Info(scrcoulint):&
        & k+qmt/2-grid is identical for to iqmt=1 grid, iqmt=",i3)') iqmt
    end if
    fsamekp=.true.
  else
    fsamekp=.false.
  end if
  if(all(abs(k_kqmtm%kqmtset%vkloff-k_kqmtm%kset%vkloff) < epslat)) then 
    if (iqmt .ne. 1) then
      write(unitout, '("Info(scrcoulint):&
        & k-qmt/2-grid is identical for to iqmt=1 grid, iqmt=",i3)') iqmt
    end if
    fsamekm=.true.
  else
    fsamekm=.false.
  end if
  !--------------------------------------------------------------------------------!

  !--------------------------------------------------------------------------------!
  ! Make W(G,G',q) Fourier coefficients and radial integrals for the PW elements
  !--------------------------------------------------------------------------------!
  ! Allocate local arrays for screened coulomb interaction and
  ! W(G,G',qr)
  allocate(scieffg(ngqmax, ngqmax, nqptr))
  scieffg(:, :, :) = zzero
  ! Phases for transformations form reduced q points to non-reduced ones.
  allocate(phf(ngqmax, ngqmax))

  ! Parallelize over reduced q-point set
  call genparidxran('q', nqptr)

  ! Set file to read the static, non-broadened RPA screening from 
  eps0dirname = 'EPS0'
  ! RA
  if(fra) then 
    ! p=-k'-k grid is identical to q=k'-k grid (zero offset)
    if(all(abs(vqoff) < epslat)) then 
      call genfilname(fileext=fileext_scr_read)
    ! p=-k'-k grid is shifted compared to q=k'-k grid (by -2*vkloff)
    else
      call genfilname(auxtype='m', fileext=fileext_scr_read)
    end if
    call genfilname(iqmt=iqmt, auxtype='m', fileext=fileext_ematrad_write)
  ! RR
  else
    ! q-grid (zero offset)
    call genfilname(fileext=fileext_scr_read)
    call genfilname(iqmt=iqmt, fileext=fileext_ematrad_write)
  end if

  write(unitout, '("Info(scrcoulint):&
    & Calculating W(G1,G2,qr) fourier coefficients")')
  call timesec(tscc0)

  do iqr = qpari, qparf ! Reduced q

    ! Locate reduced q-point in non-reduced set
    iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))

    ! Get number of G+q vectors for current q
    numgq = ngq(iqrnr)

    ! Calculate effective screened coulomb interaction
    ! by inverting the symmetrized IP dielectric matrix for a given q and
    ! 0 frequency and then multiplying
    ! it with v^{1/2} from both sides.
    filext = fileext_scr_read
    call genscclieff(iqr, iqrnr, ngqmax, numgq, scieffg(:,:,iqr))

    ! Generate radial integrals for matrix elements of plane wave
    ! and save them to disk.
    filext = fileext_ematrad_write
    call putematrad(iqr, iqrnr)
#ifndef MPI
    if(mpiglobal%rank == 0) then
      write(6, '(a,"Calculating Screened Coulomb Potential:          ", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(iqr-qpari+1)/dble(qparf-qpari+1)
      flush(6)
    end if
#endif
  end do

#ifndef MPI
  if(mpiglobal%rank == 0) then
    write(6, *)
  end if
#endif

  ! Set file extesion for later read EMATRAD in getematrad
  ! (some ranks may not participate in the qr loop above)
  filext = fileext_ematrad_write

  ! Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(set=nqptr, rlen=ngqmax*ngqmax,&
    & zbuf=scieffg, inplace=.true., comm=mpiglobal)

  if(mpiglobal%rank == 0) then
    call timesec(tscc1)
    if (input%xs%BSE%outputlevelnumber == 1) &
      & write(unitout, '("  Timing (in seconds):", f12.3)') tscc1 - tscc0
  end if
  !--------------------------------------------------------------------------------!

  !--------------------------------------------------------------------------------!
  ! CONSTRUCT W MATRIX ELEMENTS   
  !--------------------------------------------------------------------------------!

  ! Normalization factor 1/V and per k point
  pref=1.0d0/(omega*dble(nk_bse))

  ! Allocate arrays used in ematqk (does not change in the following loop)
  call ematqalloc

  ! Work arrays (allocate for maximal size over all participating k points)
  allocate(sccli(nou_bse_max, nou_bse_max))
  allocate(sccli_t2(nu_bse_max, no_bse_max, nu_bse_max, no_bse_max))
  if(fra) then 
    allocate(sccli_t1(nu_bse_max*no_bse_max, no_bse_max*nu_bse_max))
    allocate(mou(no_bse_max, nu_bse_max, ngqmax))
    allocate(muo(nu_bse_max, no_bse_max, ngqmax))
  else
    allocate(sccli_t1(nu_bse_max**2,no_bse_max**2))
    allocate(muu(nu_bse_max, nu_bse_max, ngqmax))
    allocate(moo(no_bse_max, no_bse_max, ngqmax))
  end if

  if(mpiglobal%rank == 0) then
    write(unitout, '("Info(scrcoulint): Calculating W matrix elements")')
    call timesec(tscc0)
  end if

  ! Distributed loop over combinations of non-reduced k-point
  ! that contribute to the desired energy range (or bands).
  call genparidxran('p', nkkp_bse)

  ! Gather W diagonal (only for Q=0 used)
  allocate(bsedt(3, 0:mpiglobal%procs-1))
  bsedt(1, :) = 1.d8
  bsedt(2, :) = -1.d8
  bsedt(3, :) = zzero

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

    ! Get global k point indices from
    ! the k point index of the selected set of k-points
    iknr = kmap_bse_rg(ik)
    jknr = kmap_bse_rg(jk) 

    ! Get index of ik+qmt/2 and ik-qmt/2
    ikpnr = k_kqmtp%ik2ikqmt(iknr)
    ikmnr = k_kqmtm%ik2ikqmt(iknr)

    ! Get index of jk+qmt/2 and jk-qmt/2
    jkpnr = k_kqmtp%ik2ikqmt(jknr)
    jkmnr = k_kqmtm%ik2ikqmt(jknr)

    ! Get corresponding none-reduced and reduced q-point
    ! (RR case: q = jk-ik,
    !  RA case: q = -jk-ik)
    if(fra) then 
      iq = pqmt%ikikp2ip_nr(iknr, jkpnr)
      iqr = pqmt%pset%ik2ikp(iq)
      ! Get corresponding vectors in lattice coordinated
      vq(:) = pqmt%pset%vklnr(:, iq)
      vqr(:) = pqmt%pset%vkl(:, iqr)
    else
      iq = q%ikikp2iq_nr(iknr, jknr)
      iqr = q%qset%ik2ikp(iq)
      vq(:) = q%qset%vklnr(:, iq)
      vqr(:) = q%qset%vkl(:, iqr)
    end if

    ! Check if iq is Gamma point (mod_qpoint::vqc(:,iq) has length < 1d-12)
    tq0 = tqgamma(iq)

    ! Local field effects size (Number of G+q vectors)
    numgq = ngq(iq)

    if(allocated(igqmap)) deallocate(igqmap)
    allocate(igqmap(numgq))
    if(allocated(wfc)) deallocate(wfc)
    allocate(wfc(numgq, numgq))

    ! Get radial integrals for q_r (previously calculated for reduced q set)
    call getematrad(iqr, iq)

    if(input%xs%reduceq) then 
      !! Find results  
      !! for a non reduced q-point with the help of the results
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

    ! Get ik & jk dependent band ranges for 
    ! plane wave matrix calculation (if nstlbse was used in the input all k points
    ! contribute the same number of transitions)
    ! Note: The saved ranges refer to the to k associated k_- points for the 
    !       unoccupied, and to the k_+ points for the occupied states
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

      ! Calculate M_{io jo ikp}(G, q) = <io ikp|e^{-i(q+G)r}|jo jkp>
      ! and       M_{iu ju ikm}(G, q) = <iu ikm|e^{-i(q+G)r}|ju jkm>
      ! where it is as above:
      ! q=jk-ik
      ! ikp=ik+qmt/2, ikm=ik-qmt/2
      ! jkp=jk+qmt/2, jkm=jk-qmt/2
      call getpwesrr(moo(1:ino,1:jno,1:numgq), muu(1:inu,1:jnu,1:numgq))

      ! Combine state indices for plane-wave matrix elements.

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(io,jo,j1,iu,ju,j2)
      !$OMP DO COLLAPSE(2)
      do ju = 1, jnu   ! ju
        do iu = 1, inu ! iu
          j1 = iu + (ju-1)*inu ! iuju
          ! cmuu_j = M_u1u1, M_u2u1, ..., M_uNu1, M_u1u2, ..., M_uNuM
          cmuu(j1, :) = muu(iu, ju, 1:numgq)
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP DO COLLAPSE(2)
      do jo = 1, jno   ! jo
        do io = 1, ino ! io
          j2 = io + (jo-1)*ino ! iojo
          ! cmoo_j = M_o1o1, M_o2o1, ..., M_oNo1, M_o1o2, ..., M_oNoM
          cmoo(j2, :) = moo(io, jo, 1:numgq)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! M_{iuju} -> M^*_{iuju}
      cmuu = conjg(cmuu)

      ! Allocate helper array
      allocate(zm(nuu,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclit_{j1, j2}(q)
      ! zm = cmuu * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmuu_{j1,G} wfc_{G,G'}
      !   i.e. zm_{iu_j1,ju_j1}(G',q) = \Sum_{G} M^*_{iu_j1,ju_j1}(G,q) W(G,G', q)
      call zgemm('n', 'n', nuu, numgq, numgq, zone, cmuu, nuu,&
        & wfc, numgq, zzero, zm, nuu)
      ! scclit = pref * zm * cmoo^T
      !   i.e. scclit_{j1, j2} = \Sum_{G'} zm_{j1, G'} (cmuu^T)_{G',j2}
      !   i.e. scclit_{iu_j1 ju_j1, io_j2 jo_j2} = 
      !          \Sum{G,G'} M^*_{iu_j1,ju_j1}(G,q) W(G,G',q) M_{io_j2 jo_j2}(G',q)
      call zgemm('n', 't', nuu, noo, numgq, pref, zm, nuu,&
        & cmoo, noo, zzero, sccli_t1(1:nuu,1:noo), nuu)

      deallocate(zm)
      deallocate(cmoo, cmuu)

      ! Offset in combined index for ik and jk
      jaoff = sum(kousize(1:jknr-1))
      iaoff = sum(kousize(1:iknr-1))

      ! Map back to individual band indices
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iu,ju,io,jo,j1,j2,ia,ja)
      !$OMP DO COLLAPSE(4)
      do ju = 1, jnu    ! ju
        do iu = 1, inu  ! iu
          do jo = 1, jno   ! jo
            do io = 1, ino ! io
              j1 = iu + (ju-1)*inu
              j2 = io + (jo-1)*ino
              ! scclit_{iu_j1 ju_j1, io_j2 jo_j2} -> sccli_{iu_j1 io_j2, ju_j1 jo_j2}
              sccli_t2(iu, io, ju, jo) = sccli_t1(j1, j2)
            end do
          end do
        end do
      end do
      !$OMP END DO

      ! W^RR matrix element arrays for one jk-ik=q
      ! Save only the selected transitions, i.e. W^RR_{alpha,alpha'}
      !$OMP DO COLLAPSE(2)
      do ja = 1, jnou
        do ia = 1, inou
          ju = smap_rel(1,ja+jaoff)
          jo = smap_rel(2,ja+jaoff)
          iu = smap_rel(1,ia+iaoff)
          io = smap_rel(2,ia+iaoff)
          ! scclit_{iu_j1 io_j2, ju_j1 jo_j2} -> sccliab_{iu_a io_a, ju_b jo_b}
          sccli(ia, ja) = sccli_t2(iu, io, ju, jo)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! Parallel write
      call putbsemat(scclifname, 77, ikkp, iqmt, sccli)

      ! Analyze BSE diagonal
      if(iknr .eq. jknr) then
        do ia = 1, inou
          zt1 = sccli(ia, ia)
          t1 = dble(zt1)
          ! Find minimal real part of diagonal elements
          bsedt(1, mpiglobal%rank) = min(dble(bsedt(1, mpiglobal%rank)), t1)
          ! Find maximal real part of diagonal elements
          bsedt(2, mpiglobal%rank) = max(dble(bsedt(2, mpiglobal%rank)), t1)
          ! Average diagonal element 
          bsedt(3, mpiglobal%rank) = bsedt(3, mpiglobal%rank) + zt1 / inou
        end do
      end if

    ! W^{RA}
    else
      !-------------------------------!
      ! Resonant-Anti-Resonant Part   !
      !-------------------------------!

      allocate(cmou(nou, numgq), cmuo(nuo, numgq))

      ! Using time-reversal symmetry:
      ! Calculate: N_{io ju ikp}(G, q) = <io ikp|e^{-i(G+q)r}|(ju jkm)^*>
      ! and        N_{iu jo ikm}(G, q) = <iu ikm|e^{-i(G+q)r}|(jo jkp)^*>
      ! where 
      ! q=-jk-ik
      ! ikp=ik+qmt/2, jkp=jk+qmt/2
      ! ikm=ik-qmt/2, jkm=jk-qmt/2
      call getpwesra(mou(1:ino,1:jnu,1:numgq), muo(1:inu,1:jno,1:numgq))

      ! Combine state indices for plane-wave matrix elements.

      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iu,io,ju,jo,j1,j2)
      !$OMP DO COLLAPSE(2)
      do jo = 1, jno   ! jo
        do iu = 1, inu ! iu
          j1 = iu+(jo-1)*inu
          ! cmuo_j = M_u1o1, M_u2o1, ..., M_uMo1, M_u1o2, ..., M_uMoN
          cmuo(j1, :) = muo(iu, jo, 1:numgq)
        end do
      end do
      !$OMP END DO NOWAIT
      !$OMP DO COLLAPSE(2)
      do ju = 1, jnu   ! ju
        do io = 1, ino ! io
          j2 = io+(ju-1)*ino
          ! cmou_j = M_o1u1, M_o2u1, ..., M_oNu1, M_o1u2, ..., M_oNuM
          cmou(j2, :) = mou(io, ju, 1:numgq)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! N_{iujo} -> N^*_{iujo}
      cmuo = conjg(cmuo)

      allocate(zm(nuo,numgq))

      ! Calculate matrix elements of screened coulomb interaction scclict_{j1, j2}(q)
      ! zm = cmuo * wfc
      !   i.e. zm_{j1,G'} = \Sum_{G} cmuo_{j1,G} wfc_{G,G'}
      !   i.e. zm_{iu_j1,jo_j1}(G',q) =
      !          \Sum_{G} N^*_{iu_j1,jo_j1,ikm}(G,q) W(G,G',q)
      call zgemm('n', 'n', nuo, numgq, numgq, zone, cmuo, nuo,&
        & wfc, numgq, zzero, zm, nuo)
      ! scclit = pref * zm * cmou^T
      !   i.e. scclict(j1, j2) = \Sum_{G'} zm_{j1,G'} (cmou^T)_{G',j2}
      !   i.e. scclict_{ju_j1 io_j1, jo_j2 iu_j2} =
      !   \Sum{G,G'} 
      !   N^*_{iu_j1,jo_j1,ikm}(G,q) W(G, G',q) N_{io_j2,ju_j2,ikp}(G',q)
      call zgemm('n', 't', nuo, nou, numgq, pref, zm,&
        & nuo, cmou, nou, zzero, sccli_t1(1:nuo,1:nou), nuo)

      deallocate(zm)        

      deallocate(cmou, cmuo)

      ! Save only the selected transitions
      jaoff = sum(kousize(1:jknr-1))
      iaoff = sum(kousize(1:iknr-1))

      ! Map back to individual band indices
      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iu,ju,io,jo,j1,j2,ia,ja)
      !$OMP DO COLLAPSE(4)
      do ju = 1, jnu   ! ju
        do io = 1, ino ! io
          do jo = 1, jno   ! jo
            do iu = 1, inu ! iu
              j2 = io+(ju-1)*ino ! ioju
              j1 = iu+(jo-1)*inu ! iujo
              ! scclict_{iu_j1 jo_j1 ,io_j2 ju_j2}(ik,jk)(qmt)
              ! -> scclic(iu_j1, io_j2, ju_j2, jo_j1)(ik,jk)(qmt)
              sccli_t2(iu, io, ju, jo) = sccli_t1(j1, j2)
            end do
          end do
        end do
      end do 
      !$OMP END DO

      ! W^RA matrix elements
      ! Save W^RA_{alpha,alpha'}
      !$OMP DO COLLAPSE(2)
      do ja = 1, jnou
        do ia = 1, inou
          ju = smap_rel(1,ja+jaoff)
          jo = smap_rel(2,ja+jaoff)
          iu = smap_rel(1,ia+iaoff)
          io = smap_rel(2,ia+iaoff)
          ! scclit_{iu_j1 io_j2, ju_j2 jo_j1} -> sccliab_{iu_a io_a, ju_b jo_b}
          sccli(ia, ja) = sccli_t2(iu, io, ju, jo)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      ! Parallel write
      call putbsemat(scclifname, 78, ikkp, iqmt, sccli)

    end if

    ! Deallocate G+q dependent work arrays

    deallocate(igqmap)
    deallocate(wfc)
#ifndef MPI
    if(mpiglobal%rank == 0) then
      write(6, '(a,"Calculating Screened Coulomb Matrix Elements:    ", f10.3)', advance="no")&
        & achar( 13), 100.0d0*dble(ikkp-ppari+1)/dble(pparf-ppari+1)
      flush(6)
    end if
#endif
  ! End loop over(k,kp)-pairs
  end do kkploop

  ! Gather info about resonant diagonal elements 
  !   Communicate array-parts wrt. q-points
  call mpi_allgatherv_ifc(mpiglobal%procs,rlen=3,zbuf=bsedt)
  !   BSE kernel diagonal parameters
  bsedl = minval(dble(bsedt(1, :)))
  bsedu = maxval(dble(bsedt(2, :)))
  bsedd = bsedu - bsedl
  bsed = sum(bsedt(3, :)) / nk_bse
  deallocate(bsedt)
  !   Write BSE kernel diagonal parameters
  if(mpiglobal%rank .eq. 0) call putbsediag('BSEDIAG.OUT')
#ifndef MPI 
  if(mpiglobal%rank == 0) then
    write(6, *)
  end if
#endif
  if(allocated(igqmap)) deallocate(igqmap)
  if(allocated(wfc)) deallocate(wfc)

  if(mpiglobal%rank == 0) then
    call timesec(tscc1)
    if (input%xs%BSE%outputlevelnumber == 1) &
      & write(unitout, '("  Timing (in seconds):", f12.3)') tscc1 - tscc0
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
  call ematqdealloc
  call findgntn0_clear

  call barrier(callername=trim(thisname))

  contains

    subroutine getpwesrr(moo, muu)
      use mod_variation, only: ematqk_sv
      complex(8), intent(out) :: moo(:,:,:), muu(:,:,:)

      type(bcbs) :: ematbc
      character(256) :: fileext0_save, fileext_save

      fileext0_save = filext0
      fileext_save = filext

      !--------------------------------------------------------------!
      ! Calculate M_{io jo ikp}(G, q) = <io ikp|e^{-i(q+G)r}|jo jkp> !
      !--------------------------------------------------------------!
      ! Bands
      ! Note: The saved ranges refer to the to k associated k_- points for the 
      !       unoccupied, and to the k_+ points for the occupied states
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=jno
      ematbc%il2=koulims(3,jknr)
      ematbc%iu2=koulims(4,jknr)

      ! Set EVECFV_QMTXYZ.OUT as bra state file (k+qmt/2 grid)
      usefilext0 = .true.
      iqmt0 = iqmt
      ! Set to EVECFV_QMT001.OUT if it is the same grid
      if(fsamekp) iqmt0 = iqmtgamma
      call genfilname(iqmt=iqmt0, setfilext=.true.)
      filext0 = filext

      ! Set EVECFV_QMTXYZ.OUT as ket state file (k+qmt/2 grid)
      iqmt1 = iqmt
      ! Set to EVECFV_QMT001.OUT if it is the same grid
      if(fsamekp) iqmt1 = iqmtgamma
      call genfilname(iqmt=iqmt1, setfilext=.true.)
      ! Use normal computation of the plane wave matrix elements M
      emat_ccket=.false.
      ! Set up ikmapikq to link (ikp,iq) to jkp
      ikmapikq_ptr => q%ikiq2ikp_nr
      ! Set vkl0_ptr, vkl1_ptr, ...  to k+qmt/2-grid stored in default locations 
      call setptr11()
      ! Calculate M_{o1o2,G} at fixed (k, q)
      if (input%xs%bse%xas) then
        call xasgauntgen (input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax)) 
        call ematqk_core(iq, ikpnr, moo, ematbc, 'oo')
      else
        if (.not. (input%groundstate%tevecsv)) then 
          call ematqk(iq, ikpnr, moo, ematbc)
        else
          call ematqk_sv(iq, ikpnr, moo, ematbc)
        end if
      end if
      !-----------------------------------------------------------!

      !--------------------------------------------------------------!
      ! Calculate M_{iu ju ikm}(G, q) = <iu ikm|e^{-i(q+G)r}|ju jkm> !
      !--------------------------------------------------------------!
      ! Bands
      ! Note: The saved ranges refer to the to k associated k_- points for the 
      !       unoccupied, and to the k_+ points for the occupied states
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=jnu
      ematbc%il2=koulims(1,jknr)
      ematbc%iu2=koulims(2,jknr)

      usefilext0 = .true.
      if(fsamekm) then
        ! Set EVECFV_QMT001.OUT as bra state file
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true.)
      else
        ! Set EVECFV_QMTXYZ_m.OUT as bra state file
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt0, auxtype="m", setfilext=.true.)
      end if
      filext0 = filext

      if(fsamekm) then 
        ! Set EVECFV_QMT001.OUT as ket state file
        iqmt1 = iqmtgamma
        call genfilname(iqmt=iqmt1, setfilext=.true.)
      else
        ! Set EVECFV_QMTXYZ_m.OUT as ket state file
        iqmt1 = iqmt
        call genfilname(iqmt=iqmt1, auxtype="m", setfilext=.true.)
      end if

      ! Use normal M matrix elements
      emat_ccket=.false.

      ! Set up ikmapikq to link (ikm,iq) to (jkm)
      ikmapikq_ptr => q%ikiq2ikp_nr
      ! Set vkl0_ptr, vkl1_ptr, ... to k-qmt-grid
      call setptr00()
      ! Calculate M_{u1u2,G} at fixed (k, q)
      if (.not. (input%groundstate%tevecsv)) then
        call ematqk(iq, ikmnr, muu, ematbc)
      else
        call ematqk_sv(iq, ikmnr, muu, ematbc)
      end if
      !------------------------------------------------------------------!

      filext0 = fileext0_save
      filext = fileext_save

    end subroutine getpwesrr

    subroutine getpwesra(mou, muo)
      use mod_variation, only: ematqk_sv 
      
      complex(8), intent(out) :: mou(:,:,:), muo(:,:,:)

      character(256) :: fileext0_save, fileext_save
      type(bcbs) :: ematbc

      fileext0_save = filext0
      fileext_save = filext

      !------------------------------------------------------------------------!
      ! Calculate N_{io ju ikp}(G, q) = <io ikp|e^{-i(q+G)r}|(ju jkm)^*>       !
      ! where q=-jk-ik                                                         !
      !------------------------------------------------------------------------!
      
      ! Bands
      ! Note: The saved ranges refer to the to k associated k_- points for the 
      !       unoccupied, and to the k_+ points for the occupied states
      ematbc%n1=ino
      ematbc%il1=koulims(3,iknr)
      ematbc%iu1=koulims(4,iknr)
      ematbc%n2=jnu
      ematbc%il2=koulims(1,jknr)
      ematbc%iu2=koulims(2,jknr)

      ! Set EVECFV_QMTXYZ.OUT as bra state file (k+qmt/2 grid)
      usefilext0 = .true.
      iqmt0 = iqmt
      ! Set to EVECFV_QMT001.OUT if it is the same grid
      if(fsamekp) iqmt0 = iqmtgamma
      call genfilname(iqmt=iqmt0, setfilext=.true.)
      filext0 = filext

      if(fsamekm) then 
        ! Set EVECFV_QMT001.OUT as ket state file
        iqmt1 = iqmtgamma
        call genfilname(iqmt=iqmt1, setfilext=.true.)
      else
        ! Set EVECFV_QMTXYZ_m.OUT as ket state file
        iqmt1 = iqmt
        call genfilname(iqmt=iqmt1, auxtype="m", setfilext=.true.)
      end if

      ! Use variant of the plane wave routine that calculates N 
      emat_ccket=.true.

      ! Set up non reduced ikmapikq to link (ikp,iq) to jkm
      ikmapikq_ptr => pqmt%ikip2ikp_nr

      ! Set vkl0_ptr,... to k+qmt/2-grid and vkl1_ptr,... to k-qmt/2-grid
      call setptr10()

      ! Calculate N_{ou,G} at fixed (k, q)
      if (.not. (input%xs%bse%xas)) then
        if (.not. (input%groundstate%tevecsv)) then
          call ematqk(iq, ikpnr, mou, ematbc)
        else
          call ematqk_sv(iq, ikpnr, mou, ematbc)

        end if
      else
        call xasgauntgen (input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax)) 
        call ematqk_core(iq, ikpnr, mou, ematbc, 'ou')
      end if
      !------------------------------------------------------------!

      !------------------------------------------------------------------!
      ! Calculate N_{iu jo ikm}(G, q) = <iu ikm|e^{-i(G+q)r}|(jo jkp)^*> !
      ! where q=-jk-ik                                                   !
      !------------------------------------------------------------------!

      ! Bands
      ! Note: The saved ranges refer to the to k associated k_- points for the 
      !       unoccupied, and to the k_+ points for the occupied states
      ematbc%n1=inu
      ematbc%il1=koulims(1,iknr)
      ematbc%iu1=koulims(2,iknr)
      ematbc%n2=jno
      ematbc%il2=koulims(3,jknr)
      ematbc%iu2=koulims(4,jknr)

      usefilext0 = .true.
      if(fsamekm) then 
        ! Set EVECFV_QMT001.OUT as bra state file
        iqmt0 = iqmtgamma
        call genfilname(iqmt=iqmt0, setfilext=.true.)
      else
        ! Set EVECFV_QMTXYZ_m.OUT as bra state file (k-qmt/2 grid)
        iqmt0 = iqmt
        call genfilname(iqmt=iqmt0, auxtype="m", setfilext=.true.)
      end if
      filext0 = filext

      ! Set EVECFV_QMTXYZ.OUT as ket state file (k+qmt/2 grid)
      iqmt1 = iqmt
      ! Set to EVECFV_QMT001.OUT if it is the same grid
      if(fsamekp) iqmt1 = iqmtgamma
      call genfilname(iqmt=iqmt1, setfilext=.true.)

      ! Use variant of the plane wave routine that calculates N 
      emat_ccket=.true.

      ! Set up ikmapikq to link (ikm,ip) to jkp
      ikmapikq_ptr => pqmt%ikip2ikp_nr

      ! Set vkl0_ptr,... to k-qmt/2-grid and vkl1_ptr,... to k+qmt/2-grid
      call setptr01()

      ! Calculate N_{uo,G} at fixed (k, q)
      if (.not. (input%xs%bse%xas)) then
        if (.not. (input%groundstate%tevecsv)) then
          call ematqk(iq, ikmnr, muo, ematbc)
        else
          call ematqk_sv(iq,ikmnr, muo,ematbc)
        end if
      else
        call xasgauntgen (input%xs%lmaxemat, Max(input%groundstate%lmaxapw, lolmax)) 
        call ematqk_core(iq, ikmnr, muo, ematbc, 'uo')
      end if
      !-------------------------------------------------------------!

      filext0 = fileext0_save
      filext = fileext_save
    end subroutine getpwesra

end subroutine scrcoulint
!EOC

