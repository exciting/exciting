! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: kernxc_bse
! !INTERFACE:
!
subroutine kernxc_bse
! !USES:
  use modmpi
  use modinput
  use mod_APW_LO, only: lolmax
  use mod_qpoint, only: nqpt
  use mod_kpoint, only: nkptnr
  use mod_misc, only: task
  use mod_lattice, only: omega
  use mod_constants, only: zi, zone, zzero
#ifdef TETRA      
  use modtetra
#endif
  use modxs, only: xsgnt, unitout, ikmapikq,&
    & istocc0, istunocc0, istocc, istunocc,&
    & isto0, istu0, isto, istu,&
    & nst1, nst2, nst3, nst4,&
    & istl1, istl2, istl3, istl4,&
    & istu1, istu2, istu3, istu4,&
    & ksgap, wpari, wparf, nwdf,&
    & bzsampl, ngq, xiou, xiuo,&
    & pmou, pmuo, deou, deuo,&
    & docc12, docc21, bsed, torfxc,&
    & iqmt1, iqmt0, iqmtgamma, bcbs,&
    & sta1, sto1, sta2, sto2,&
    & usefilext0, filext0
  use m_xsgauntgen
  use m_findgntn0
  use m_writegqpts
  use m_genwgrid
  use m_xszoutpr3
  use m_getpemat
  use m_getunit
  use m_genfilname
  use m_putgetbsemat
  use m_b_ematqk
  use modbse
! !INPUT/OUTPUT PARAMETERS:
!   oct   : optical diagonal tensor component (in,integer)
! !DESCRIPTION:
!   BSE-kernel of A. Marini, Phys. Rev. Lett. 91, 256402 (2003)
!
! !REVISION HISTORY:
!   Created March 2008 (Sagmeister)
!EOP
!BOC

  implicit none

  ! local variables
  character(*), parameter :: thisname = 'kernxc_bse'
  integer(4), parameter :: iqmt = 1, noptc = 3

  character(256) :: filnam2, filnam3, filnam4
  integer(4) :: iw, wi, wf, nwdfp, n, reclen, un, un2, un3, j1, j2, oct
  integer(4) :: ikkp, iknr, jknr, iknrq, jknrq, igq1, igq2
  integer(4) :: ist1, ist2, ist3, ist4, nst12, nst34, nst13, nst24
  real(8) :: t1, brd
  real(8) :: cpu_init1offs, cpu_ematrad, cpu_ematqalloc, cpu_ematqk1
  real(8) :: cpu_ematqdealloc, cpu_clph, cpu_suma, cpu_write
  complex(8) :: zt1

  ! allocatable arrays
  real(8), allocatable :: dek(:, :), dok(:, :), scisk(:, :)
  real(8), allocatable :: dekp(:, :), dokp(:, :), sciskp(:, :)
  real(8), allocatable :: deval(:, :, :), docc(:, :, :), scis(:, :, :)
  real(8), allocatable :: dde(:, :)
  complex(8), allocatable :: zmr(:, :), zmq(:, :), zmra(:, :), zmqa(:, :)
  complex(8), allocatable :: scclit(:, :), scclith(:, :), sccli(:, :, :, :)
  complex(8), allocatable :: emat(:, :, :, :), emata(:, :, :, :)
  complex(8), allocatable :: den1(:), den2(:), den1a(:), den2a(:)
  complex(8), allocatable :: emat12p(:, :), emat12pa(:, :)
  complex(8), allocatable :: emat12k(:, :, :), emat12kp(:, :, :)
  complex(8), allocatable :: emat12ka(:, :, :), emat12kpa(:, :, :)
  complex(8), allocatable :: residr(:, :), residq(:, :), osca(:, :), oscb(:, :)
  complex(8), allocatable :: residra(:, :), residqa(:, :), oscaa(:, :), oscba(:, :)
  complex(8), allocatable :: fxc(:, :, :), w(:), bsedg(:, :), bufou(:, :, :),&
    & bufuo(:, :, :), pufou(:, :, :), pufuo(:, :, :)

  character(256) :: sfname, sinfofname
  type(bcbs) :: bc
  logical :: sfcmpt, sfid

  ! external functions
  integer(4), external :: idxkkp, l2int
  logical, external :: tqgamma

  ! check that if kohn-sham response is time-ordered, so is the setting for the kernel
  if(input%xs%tddft%tordfxc .neqv. input%xs%tddft%torddf) then
    write(*,*)
    write(*,'("Error(kernxc_bse): both, the kohn-sham response function")')
    write(*,'(" and the bse-derived xc kernel have to be either causal or time-ordered.")')
    write(*,*)
    call terminate
  end if

  brd = input%xs%broad

  call init0
  call init1
  call init2

  ! save variables for the gamma q-point
  call xssave0

  ! Read Fermi energy from file EFERMI
  ! Use EFERMI_QMT001.OUT (corresponding to the xs groundstate run for the unshifted k grid)
  call genfilname(iqmt=iqmtgamma, setfilext=.true.)
  call readfermi

  ! Set ist* variables and ksgap in modxs using findocclims
  ! This also reads in 
  ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
  ! modxs:evalsv0, modxs:occsv0
  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Inspecting occupations...'
  call flushifc(unitout)
  call setranges_modxs(iqmt)

  ! Select relevant transitions for the construction
  ! of the BSE hamiltonian
  ! Also sets nkkp_bse, nk_bse 
  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Selecting transitions...'
  call flushifc(unitout)
  call select_transitions(iqmt, serial=.false.)

  ! generate gaunt coefficients
  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
    & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))

  ! find indices for non-zero gaunt coefficients
  call findgntn0(max(input%xs%lmaxapwwf, lolmax), &
    & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)

  write(unitout, '(a)') 'Info(' // thisname // '):&
    & Gaunt coefficients generated within lmax values:'
  write(unitout, '(a, i8)') "lmax1 = lmaxapw =", input%groundstate%lmaxapw
  write(unitout, '(a, i8)') "lmax2 = lmaxemat=", input%xs%lmaxemat
  write(unitout, '(a, i8)') "lmax3 = lmaxapw =", input%groundstate%lmaxapw
  write(unitout, '(a, i6)') 'Info(' // thisname // '): Number of q-points: ', nqpt
  call flushifc(unitout)

  call findocclims(0, ikmapikq(:,1), istocc0, istunocc0, isto0, isto, istu0, istu)
  istunocc = istunocc0
  istocc = istocc0

  ! only for systems with a gap in energy
  if( .not. ksgap) then
    write(*,*)
    write(*, '("Error(", a, "):&
      & Screened coulomb interaction works only for systems with ks-gap.")')&
      & trim(thisname)
    write(*,*)
    call terminate
  end if

  ! Set BSE ranges (o absolute, u relative) (yes, it is indeed madness)
  sta1 = input%xs%bse%nstlbse(1)
  sto1 = input%xs%bse%nstlbse(2)
  sta2 = input%xs%bse%nstlbse(3)
  sto2 = input%xs%bse%nstlbse(4)

  ! Selects 12 = ou , 34 = uo ranges. Using BSE range.
  input%xs%emattype = 1
  call ematbdcmbs(input%xs%emattype)

  nst12 = nst1 * nst2
  nst34 = nst3 * nst4
  nst13 = nst1 * nst3
  nst24 = nst2 * nst4

  call genparidxran('w', nwdf)

  ! sampling type for brillouin zone sampling
  bzsampl = l2int(input%xs%tetra%tetradf) ! defaults to false

  ! limits for w-points
  wi = wpari
  wf = wparf
  nwdfp = wparf - wpari + 1

  ! matrix size for local field effects (first q-point is gamma-point)
  n = ngq(iqmt)

  allocate(bufou(nst1, nst2, n))
  allocate(bufuo(nst2, nst1, n))
  allocate(pufou(3, nst1, nst2))
  allocate(pufuo(3, nst2, nst1))

  ! allocate global arrays
  if(allocated(xiou)) deallocate(xiou)
  allocate(xiou(nst1, nst2, n))
  if(allocated(xiuo)) deallocate(xiuo)
  allocate(xiuo(nst2, nst1, n))
  if(allocated(pmou)) deallocate(pmou)
  allocate(pmou(3, nst1, nst2))
  if(allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3, nst2, nst1))
  if(allocated(deou)) deallocate(deou)
  allocate(deou(nst1, nst2))
  if(allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst2, nst1))
  if(allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1, nst2))
  if(allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst2, nst1))

  ! allocate local arrays
  allocate(emat12p(nst12,-3:n), zmr(nst12, nst12), zmq(nst12, nst12))
  allocate(emat12pa(nst12,-3:n), zmra(nst12, nst12), zmqa(nst12, nst12))
  allocate(dek(nst1, nst2), dekp(nst1, nst2), dde(nst1, nst2))
  allocate(dok(nst1, nst2), dokp(nst1, nst2))
  allocate(scisk(nst1, nst2), sciskp(nst1, nst2))
  allocate(fxc(-3:n,-3:n, nwdf))
  allocate(scclit(nst12, nst12))
  allocate(scclith(nst12, nst12))
  allocate(emat12k(-3:n, nst1, nst2), emat12kp(nst1, nst2,-3:n))
  allocate(residr(nst12,-3:n), residq(nst12,-3:n))
  allocate(emat12ka(-3:n, nst2, nst1), emat12kpa(nst2, nst1,-3:n))
  allocate(residra(nst12,-3:n), residqa(nst12,-3:n))
  allocate(w(nwdf))
  allocate(osca(-3:n,-3:n), oscb(-3:n,-3:n))
  allocate(oscaa(-3:n,-3:n), oscba(-3:n,-3:n))
  allocate(den1(nwdf), den2(nwdf), den1a(nwdf), den2a(nwdf))
  fxc(:, :, :) = zzero
  allocate(emat(nst1, nst2, n, nkptnr))
  allocate(emata(nst2, nst1, n, nkptnr))
  allocate(deval(nst1, nst2, nkptnr))
  allocate(docc(nst1, nst2, nkptnr))
  allocate(scis(nst1, nst2, nkptnr))
  allocate(bsedg(nst1, nst2))
  allocate(sccli(nst1, nst2, nst1, nst2))
  sccli(:, :, :, :) = zzero

  if((input%xs%tddft%fxctypenumber .eq. 7) .or. &
    &(input%xs%tddft%fxctypenumber .eq. 8)) then

    call getbsediag
    write(unitout, '("Info(", a, "): Read diagonal of bse kernel")') trim(thisname)
    write(unitout, '(" mean value : ", 2g18.10)') bsed
  end if

  ! generate energy grid
  call genwgrid(nwdf, input%xs%energywindow%intv, input%xs%tddft%acont,&
    & 0.d0, w_cmplx=w)

  ! precalculate matrix elements
  !   calculate radial integrals
  call ematrad(iqmt)
  !   allocate work arrays
  call ematqalloc
  call setptr01
  ikmapikq_ptr => ikmapikq

  !---------------------------!
  !     loop over k-points    !
  !---------------------------!
  usefilext0 = .true.
  iqmt0 = iqmtgamma
  call genfilname(iqmt=iqmt0, fileext=filext0)

  iqmt1 = iqmtgamma
  call genfilname(iqmt=iqmt1, setfilext=.true.)

  do iknr = 1, nkptnr

    call chkpt(3,(/ task, 1, iknr /),&
      & 'task, sub, k - point; generate matrix elements of plane wave')

    iknrq = ikmapikq(iknr, iqmt)

    ! Get ou
    if(allocated(xiou)) deallocate(xiou)
    allocate(xiou(nst1, nst2, ngq(iqmt)))
    bc%n1 = nst1
    bc%n2 = nst2
    bc%il1 = istl1
    bc%il2 = istl2
    bc%iu1 = istu1
    bc%iu2 = istu2
    call b_ematqk(iqmt, iknr, xiou, bc)

    ! Get uo
    if(allocated(xiuo)) deallocate(xiuo)
    allocate(xiuo(nst3, nst4, ngq(iqmt)))
    bc%n1 = nst3
    bc%n2 = nst4
    bc%il1 = istl3
    bc%il2 = istl4
    bc%iu1 = istu3
    bc%iu2 = istu4
    call b_ematqk(iqmt, iknr, xiuo, bc)

    emat(:, :, :, iknr) = xiou(:, :, :)
    emata(:, :, :, iknr) = xiuo(:, :, :)

    deallocate(xiou, xiuo)

    call getdevaldoccsv(iqmt, iknr, iknrq, istl1, istu1,&
      & istl2, istu2, deou, docc12, scisk)

    call getdevaldoccsv(iqmt, iknr, iknrq, istl2, istu2,&
      & istl1, istu1, deuo, docc21, sciskp)

    deval(:, :, iknr) = deou(:, :)
    docc(:, :, iknr) = docc12(:, :)
    scis(:, :, iknr) = scisk(:, :)

  end do

  ! Get info about saved W matrix
  call genfilname(basename=scclifbasename, iqmt=iqmt, filnam=sfname)
  sinfofname = trim(infofbasename)//'_'//trim(sfname)
  call b_getbseinfo(trim(sinfofname), iqmt, fcmpt=sfcmpt, fid=sfid)
  write(unitout, '("  Reading info from ", a)') trim(sinfofname)
  write(unitout, '("  Reading W from ", a)') trim(sfname)
  write(unitout, '("    W compatible:",l," W identical:",l)') sfcmpt, sfid
  if(.not. sfcmpt) call terminate

  !-------------------------------!
  !     loop over(k,kp) pairs     !
  !-------------------------------!
  if(allocated(xiou)) deallocate(xiou)
  if(allocated(xiuo)) deallocate(xiuo)

  ikkp = 0

  ! first k-point
  do iknr = 1, nkptnr

    call chkpt(3,(/ task, 3, iknr /),&
      & 'task, sub, k-point; bse-fxc-kernel')

    iknrq = ikmapikq(iknr, iqmt)

    bsedg(:, :) = bsed

    allocate(xiou(nst1, nst2, n))
    allocate(xiuo(nst2, nst1, n))
    xiou(:, :, :) = emat(:, :, :, iknr)
    xiuo(:, :, :) = emata(:, :, :, iknr)
    deou(:, :) = deval(:, :, iknr)
    docc12(:, :) = docc(:, :, iknr)

    ! apply gauge wrt. symmetrized coulomb potential
    call getpemat(iqmt, iknr, 'PMAT_XS.OUT', '',&
      & m12=bufou, p12=pufou, m34=bufuo, p34=pufuo)

    dek(:, :) = deou(:, :)
    dok(:, :) = docc12(:, :)

    ! add bse diagonal
    scisk(:, :) = scis(:, :, iknr) + bsedg(:, :)

    ! assign optical components
    do oct = 1, noptc
      emat12k(-oct, :, :) = pufou(oct, :, :)
      emat12ka(-oct, :, :) = pufuo(oct, :, :)
    end do

    do igq1 = 1, n
      emat12k(igq1, :, :) = bufou(:, :, igq1)
      emat12ka(igq1, :, :) = bufuo(:, :, igq1)
    end do

    deallocate(xiou, xiuo)

    residr(:, :) = zzero
    residq(:, :) = zzero
    residra(:, :) = zzero
    residqa(:, :) = zzero

    ! second k-point
    do jknr = 1, nkptnr

      jknrq = ikmapikq(jknr, iqmt)

      cpu_init1offs = 0.d0
      cpu_ematrad = 0.d0
      cpu_ematqalloc = 0.d0
      cpu_ematqk1 = 0.d0
      cpu_ematqdealloc = 0.d0
      cpu_clph = 0.d0
      cpu_suma = 0.d0
      cpu_write = 0.d0

      if(iknr .le. jknr) then
        ! index for upper triangle
        ikkp = idxkkp(iknr, jknr, nkptnr)
      else
        ! swapped index for lower triangle
        ikkp = idxkkp(jknr, iknr, nkptnr)
      end if

      allocate(xiou(nst1, nst2, n))
      allocate(xiuo(nst2, nst1, n))
      xiou(:, :, :) = emat(:, :, :, jknr)
      xiuo(:, :, :) = emata(:, :, :, jknr)
      deou(:, :) = deval(:, :, jknr)
      docc12(:, :) = docc(:, :, jknr)

      ! apply gauge wrt. symmetrized coulomb potential
      call getpemat(iqmt, jknr, 'PMAT_XS.OUT', '',&
        & m12=bufou, p12=pufou, m34=bufuo, p34=pufuo)

      dekp(:, :) = deou(:, :)
      dokp(:, :) = docc12(:, :)
      sciskp(:, :) = scis(:, :, jknr)

      ! assign optical component
      do oct = 1, noptc
        emat12kp(:, :,-oct) = pufou(oct, :, :)
        emat12kpa(:, :,-oct) = pufuo(oct, :, :)
      end do

      emat12kp(:, :, 1:) = bufou(:, :, :)
      emat12kpa(:, :, 1:) = bufuo(:, :, :)

      deallocate(xiou, xiuo)


      if(iknr .le. jknr) then

        ! Get ikkp block of W
        call b_getbsemat(trim(sfname), iqmt, ikkp,&
          & scclit(1:nst12,1:nst12), check=.false.)

      else

        ! Get ikkp block of W
        call b_getbsemat(trim(sfname), iqmt, ikkp,&
          & scclith(1:nst12,1:nst12), check=.false.)

        ! use hermitian property for lower triangle
        do ist1 = 1, nst12
          do ist2 = 1, nst12

            scclit(ist1, ist2) = conjg(scclith(ist2, ist1))

          end do
        end do

      end if

      ! proper sign of screened coulomb interaction
      scclit = -scclit

      ! set diagonal of bethe-salpeter kernel to zero
      ! (cf. a. marini, prl 2003)
      if(iknr .eq. jknr) then
        do ist1 = 1, nst12
          scclit(ist1, ist1) = zzero
        end do
      end if

      ! Expand to 4 indices
      ! (saved is a matrix in the contracted indces uo, where u is the fastest index)
      j1=0
      do ist1 = 1, nst1 ! o1
        do ist2 = 1, nst2 ! u1
          j1 = j1+1

          j2=0
          do ist3 = 1, nst1 ! o2
            do ist4 = 1, nst2 ! u2
              j2 = j2+1

              ! W_o1u1,o2u2(ikkp) <-- W_j1,j2(ikkp)
              sccli(ist1, ist2, ist3, ist4) = scclit(j1, j2)

            end do
          end do

        end do
      end do

      j1 = 0
      do ist2 = 1, nst2
        do ist1 = 1, nst1
          j1 = j1 + 1
          emat12p(j1, :) = conjg(emat12kp(ist1, ist2, :))
          emat12pa(j1, :) = conjg(emat12kpa(ist2, ist1, :))
        end do
      end do

      ! map
      j2 = 0
      do ist2 = 1, nst2 ! u1
        do ist1 = 1, nst1 ! o1
          j2 = j2 + 1
          j1 = 0
          do ist4 = 1, nst2 ! u2
            do ist3 = 1, nst1 ! o2
              j1 = j1 + 1
              zt1 = sccli(ist1, ist2, ist3, ist4)
              ! four point energy difference
              t1 = dekp(ist3, ist4) - dek(ist1, ist2)
              ! arrays for r- and q-residuals

              if(abs(t1) .ge. input%xs%tddft%fxcbsesplit) then
                zmr(j2, j1) = zt1 / t1
                zmq(j2, j1) = zzero
                zmra(j2, j1) = conjg(zt1) / t1
                zmqa(j2, j1) = zzero
              else
                zmr(j2, j1) = zzero
                zmq(j2, j1) = zt1
                zmra(j2, j1) = zzero
                zmqa(j2, j1) = conjg(zt1)
              end if

            end do
          end do

        end do
      end do

      ! calculate residual "r"; partial fraction decomposition without double poles
      ! (cf. a. marini, phys. rev. lett. 91, 256402(2003))
      residr = residr + matmul(zmr, emat12p)
      residra = residra + matmul(zmra, emat12pa)

      ! calculate residual "q"; double poles part
      ! (cf. a. marini, phys. rev. lett. 91, 256402(2003))
      residq = residq + matmul(zmq, emat12p)
      residqa = residqa + matmul(zmqa, emat12pa)

    ! end inner loop over k-points
    end do

    !--------------------------!
    !     set up bse-kernel    !
    !--------------------------!
    t1 = 1.d0 /(nkptnr*omega)

    do ist2 = 1, nst2
      do ist1 = 1, nst1

        osca(:, :) = zzero
        oscb(:, :) = zzero
        oscaa(:, :) = zzero
        oscba(:, :) = zzero
        j1 = ist1 +(ist2-1) * nst1

        ! set up inner part of kernel
        ! generate oscillators
        call xszoutpr3(n+noptc+1, n+noptc+1, zone,&
          & emat12k(:, ist1, ist2), residr(j1, :), osca)
        call xszoutpr3(n+noptc+1, n+noptc+1, zone,&
          & emat12ka(:, ist2, ist1), residra(j1, :), oscaa)

        ! add hermitian transpose
        forall(igq1=-3:n, igq2=-3:n)
          osca(igq1, igq2) = osca(igq1, igq2) + conjg(osca(igq2, igq1))
          oscaa(igq1, igq2) = oscaa(igq1, igq2) + conjg(oscaa(igq2, igq1))
        end forall

        call xszoutpr3(n+noptc+1, n+noptc+1, zone,&
          & emat12k(:, ist1, ist2), residq(j1, :), oscb)
        call xszoutpr3(n+noptc+1, n+noptc+1, zone,&
          & emat12ka(:, ist2, ist1), residqa(j1, :), oscba)

        ! set up energy denominators
        den1(:) = 2.d0 * t1 /(w(:)+scisk(ist1, ist2)+dek(ist1, ist2)+zi*brd)
        den2(:) = 2.d0 * t1 /(w(:)+scisk(ist1, ist2)+dek(ist1, ist2)+zi*brd)**2
        den1a(:) = 2.d0 * t1 /(w(:)+scisk(ist1, ist2)-dek(ist1, ist2)+torfxc*zi*brd)
        den2a(:) = -2.d0 * t1 /(w(:)+scisk(ist1, ist2)-dek(ist1, ist2)+torfxc*zi*brd)**2

        ! update kernel
        do iw = 1, nwdf

          ! resonant contributions
          fxc(:, :, iw) = fxc(:, :, iw) + osca(:, :) * den1(iw) + oscb(:, :) * den2(iw)

          ! antiresonant contributions
          if(input%xs%tddft%aresfxc) then
            fxc(:, :, iw) = fxc(:, :, iw) + oscaa(:, :) * den1a(iw) &
              &+ oscba(:, :) * den2a(iw)
          end if

        end do

      ! end loop over states #1
      end do
    ! end loop over states #2
    end do
  ! end outer loop over k-points
  end do

  ! filename for xc-kernel (ascii)
  call genfilname(basename='KERNXC_BSE', asc=.true., bzsampl=bzsampl,&
    & acont=input%xs%tddft%acont, nar= .not. input%xs%tddft%aresfxc,&
    & tord=input%xs%tddft%tordfxc, iqmt=iqmt, filnam=filnam2)
   
  call getunit(un)
  open(un, file=trim(filnam2), form='formatted', action='write', status='replace')

  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE', asc=.false., bzsampl=bzsampl,&
    & acont=input%xs%tddft%acont, nar= .not. input%xs%tddft%aresfxc,&
    & tord=input%xs%tddft%tordfxc, iqmt=iqmt, filnam=filnam3)

  inquire(iolength=reclen) n, fxc(-3:-1,-3:-1, 1), fxc(-3:-1, 1:, 1),&
    & fxc(1:,-3:-1, 1), fxc(1:, 1:, 1)
  call getunit(un2)
  open(un2, file=trim(filnam3), form='unformatted', action='write',&
    & status='replace', access='direct', recl=reclen)

  ! filename for xc-kernel
  call genfilname(basename='FXC_BSE_HEAD', asc=.false., &
    & bzsampl=bzsampl, acont=input%xs%tddft%acont, nar= .not. input%xs%tddft%aresfxc,&
    & tord=input%xs%tddft%tordfxc, iqmt=iqmt, filnam=filnam4)

  call getunit(un3)
  open(un3, file=trim(filnam4), form='formatted', action='write', status='replace')

  do iw = 1, nwdf
    write(un2, rec=iw) n, fxc(-3:-1,-3:-1, iw), fxc(-3:-1, 1:, iw),&
      & fxc(1:,-3:-1, iw), fxc(1:, 1:, iw)
    write(un3, '(i6, 2x, g18.10, 2x, 6g18.10)') iw, dble(w(iw)),&
      & (fxc(-oct,-oct, iw), oct=1, noptc)
  end do

  do iw = 1, nwdf, 10
    do igq1 = - noptc, n
       do igq2 = - noptc, n
          write(un, '(3i6, 3g18.10)') iw, igq1, igq2,&
            & fxc(igq1, igq2, iw), abs(fxc(igq1, igq2, iw))
       end do
    end do
  end do

  close(un)
  close(un2)
  close(un3)

  ! deallocate
  deallocate(den1, den2, den1a, den2a)
  deallocate(emat12p, zmr, zmq, dek, dekp, dde, dok, dokp, scisk, fxc)
  deallocate(sccli, scclit, scclith, emat12k, emat12kp, residr, residq, w, osca, oscb)
  deallocate(emat, deval, docc, scis)

  ! deallocate antiresonant parts
  deallocate(emata, emat12pa, emat12ka, emat12kpa, residra, residqa, zmra, zmqa)
  deallocate(oscaa, oscba)

  deallocate(bsedg)
  deallocate(bufou, bufuo, pufou, pufuo)

end subroutine kernxc_bse
!EOC
