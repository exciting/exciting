! Copyright (C) 2005-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dfq
! !INTERFACE:
subroutine dfq(iq)
! !USES:
  use modinput, only: input
  use modmpi, only: procs, rank, barrier
  use mod_misc, only: task
  use mod_constants, only: zzero, zone, zi, krondelta
  use mod_kpoint, only: nkpt, wkpt
  use mod_lattice, only: omega
  use modxs, only: tfxcbse, tscreen, bzsampl, wpari,&
                & wparf, ngq, fnwtet, fnpmat,&
                & fnetim, fnxtim, fnemat, fnchi0,&
                & fnchi0_t, qvkloff, istocc0, istocc,&
                & istunocc0, istunocc, isto0, isto,&
                & istu0, istu, unitout, nst1,&
                & nst2, nst3, nst4, istl1,&
                & istu1, istl2, istu2, istl3,&
                & istu3, istl4, istu4, deou,&
                & deuo, docc12, docc21, xiou,&
                & xiuo, pmou, pmuo, nwdf,&
                & bsed, ikmapikq, tordf, symt2
#ifdef TETRA      
  use mod_eigenvalue_occupancy, only: nstsv, evalsv, efermi 
  use mod_qpoint, only: vql
  use modtetra
#endif
  ! One subroutine modules
  use m_genwgrid
  use m_getpemat
  use m_dftim
  use m_gettetcw
  use m_putx0
  use m_getunit
  use m_writevars
  use m_filedel
  use m_genfilname
! !DESCRIPTION:
!   Calculates the symmetrized Kohn-Sham response function $\tilde{\chi}^0_{\bf{GG'}}
!   ({\bf q},\omega)$ for one ${\bf q}$-point according to
!   $$  \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega) = \sum_{ou{\bf k}} \left[
!      \tilde{M}^{\bf G}_{ou{\bf k}}({\bf q}) \tilde{M}^{\bf G'}_{ou{\bf k}}({\bf q})^*
!      w_{ou{\bf k}}({\bf q},\omega) +
!      \tilde{M}^{\bf G}_{uo{\bf k}}({\bf q}) \tilde{M}^{\bf G'}_{uo{\bf k}}({\bf q})^*
!      w_{uo{\bf k}}({\bf q},\omega) \right]
!   $$
!   It is related to the Kohn-Sham response function
!    $\chi^0_{\bf{GG'}}({\bf q},\omega)$ by
!   $$ \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega) = v^{\frac{1}{2}}_{\bf G}({\bf q})
!      \chi^0_{\bf{GG'}}({\bf q},\omega)
!      v^{\frac{1}{2}}_{\bf G'}({\bf q}) $$
!   and is well defined in the limit as ${\bf q}$ goes to zero.
!   The symmetrized matrix elements are defined as
!   $$   \tilde{M}^{\bf G}_{ou{\bf k}}({\bf q}) =
!         v^{\frac{1}{2}}_{\bf G}({\bf q})
!         M^{\bf G}_{ou{\bf k}}({\bf q}), $$
!   where
!   $$   M^{\bf G}_{ou{\bf k}}({\bf q}) =
!        \langle o{\bf k}|e^{-i({\bf{G+q}}){\bf r}}|u{\bf k+q} \rangle. $$
!   For ${\bf G}=0$ we have to consider three vectors stemming from the limits
!   as ${\bf q}\rightarrow 0$ along the Cartezian basis vectors ${\bf e_i}$,
!   i.e., we can think of ${\bf 0}_1,{\bf 0}_2,{\bf 0}_3$ in place of ${\bf 0}$.
!   The weights $w_{\rm ou{\bf k}}({\bf q},\omega)$ are defined as
!   $$   w_{nm{\bf k}}({\bf q},\omega) = \lambda_{\bf k}
!        \frac{f_{n{\bf k}}-f_{m{\bf k+q}}}
!        {\varepsilon_{n{\bf k}}-\varepsilon_{m{\bf k+q}}+
!        \Delta_{n{\bf k}}-\Delta_{m{\bf k+q}}+\omega+i\eta} $$
!   in the case where we use a Lorentzian broadening $\eta$.
!   In the above expression $\lambda_{\bf k}$ is the weight of the
!   ${\bf k}$-point, $\varepsilon_{n{\bf k}}$ and
!   $\varepsilon_{m{\bf k+q}}$ are the DFT Kohn-Sham energies,
!   $\Delta_{n{\bf k}}$ and $\Delta_{m{\bf k+q}}$ are the scissors corrections
!   that are non-zero in the case where $m{\bf k+q}$ corresponds to a
!   conduction state.
!   The indices $o$ and $u$ denote {\it at least partially occupied} and
!   {\it at least partially unoccupied} states, respectively.
!   The symmetrized Kohn-Sham response function can also be calculated
!   for imaginary frequencies $i\omega$ without broadening $\eta$. In this
!   case the replacement
!   $$ \omega+i\eta \mapsto i\omega $$
!   is applied to the expressions for the weights.
!   Optionally, the weights can be calculated with the help of the linear
!   tetrahedron method (including Bloechl's correction).
!   This routine can be run with MPI parallelization for energies.
!
! !REVISION HISTORY:
!   Created March 2005 (Sagmeister)
!   Added band and k-point analysis, 2007-2008 (Sagmeister)
!EOP
!BOC

  implicit none

  ! Arguments
  integer, intent(in) :: iq

  ! Local variables
  character(*), parameter :: thisnam = 'dfq'
  character(256) :: fnscreen
  complex(8) :: zt1, winv
  complex(8), allocatable :: w(:)
  complex(8), allocatable :: chi0 (:, :, :), hdg(:, :, :)
  complex(8), allocatable :: chi0w(:, :, :, :), chi0h(:, :, :), eps0(:, :, :)
  complex(8), allocatable :: chi0hahc(:, :)
  complex(8), allocatable :: wou(:,:,:), wuo(:,:,:), wouw(:), wuow(:), wouh(:), wuoh(:)
  complex(8), allocatable :: zvou(:), zvuo(:), chi0hs(:, :, :), bsedg(:, :), scis12c(:, :), scis21c(:, :),zm(:,:,:)
  real(8), parameter :: epstetra = 1.d-8
  real(8) :: ta,tb,tc,td
  real(8), allocatable :: wreal(:) 
  real(8), allocatable :: scis12 (:, :), scis21 (:, :)
  real(8) :: brd, cpu0, cpu1, cpuread, cpuosc, cpuupd, cputot, wintv(2), wplas, wrel
  integer :: n, j, i1, i2, ik, ikq, igq, iw, wi, wf, ist1, ist2, nwdfp
  integer :: oct1, oct2, un
  logical :: tq0

  ! External functions
  logical, external :: tqgamma, transik, transijst

  call timesec(ta)     

  ! Analytic continuation not allowed with screen
  if(input%xs%tddft%acont .and. tscreen) then
    write(*,*)
    write(*, '("Error(", a, "): analytic continuation&
     & does not work for screening")')
    write(*,*)
    call terminate
  end if

  tfxcbse = ((input%xs%tddft%fxctypenumber .eq. 7) .or. &
     & (input%xs%tddft%fxctypenumber .eq. 8)) .and. ( .not. tscreen)

  ! Sampling of Brillouin zone
  bzsampl = 0
  if(input%xs%tetra%tetradf) bzsampl = 1

  ! Initial and final w-point
  wi = wpari
  wf = wparf
  nwdfp = wf - wi + 1

  ! Matrix size for response function
  ! Get number of G+q vectors for current q
  n = ngq(iq)

  ! Zero broadening for analytic continuation
  brd = input%xs%broad
  if(input%xs%tddft%acont) brd = 0.d0

  ! Zero broadening for dielectric matrix (w=0) for band-gap systems
  ! Task 430 is 'screen'
  if(task .eq. 430) brd = 0.d0

  ! File extension for q-point
  if( .not. tscreen) call genfilname(iqmt=iq, setfilext=.true.)

  ! Filenames for output
  if(tscreen) then
    call genfilname(basename='TETW', iq=iq, appfilext=.true., filnam=fnwtet)
    call genfilname(basename='PMAT', appfilext=.true., filnam=fnpmat)
    call genfilname(basename='SCREEN', bzsampl=bzsampl, iq=iq, filnam=fnscreen)
    call genfilname(nodotpar=.true., basename='EMAT_TIMING', iq=iq,&
     & etype=input%xs%emattype, procs=procs, rank=rank, appfilext=.true., filnam=fnetim)
    call genfilname(nodotpar=.true., basename='X0_TIMING', iq=iq,&
     & procs=procs, rank=rank, appfilext=.true., filnam=fnxtim)
  else
    call genfilname(basename='TETW', iqmt=iq, filnam=fnwtet)
    call genfilname(basename='PMAT_XS', filnam=fnpmat)
    call genfilname(basename='EMAT', iqmt=iq, filnam=fnemat)
    call genfilname(nodotpar=.true., basename='X0_TIMING', bzsampl=bzsampl,&
     & iqmt=iq, procs=procs, rank=rank, filnam=fnxtim)
    call genfilname(basename='X0', bzsampl=bzsampl,&
     & acont=input%xs%tddft%acont, nar= .not. input%xs%tddft%aresdf,&
     & tord=input%xs%tddft%torddf, markfxcbse=tfxcbse, iqmt=iq, filnam=fnchi0)
    call genfilname(basename='X0', bzsampl=bzsampl,&
     & acont=input%xs%tddft%acont, nar= .not. input%xs%tddft%aresdf,&
     & tord=input%xs%tddft%torddf, markfxcbse=tfxcbse, iqmt=iq,&
     & procs=procs, rank=rank, filnam=fnchi0_t)
  end if

  ! Remove timing files from previous runs
  call filedel(trim(fnxtim))

  ! Calculate k+q and g+k+q related variables
  call init1offs(qvkloff(1, iq))

  ! Generate link array for tetrahedra
  if(input%xs%tetra%tetradf) then
#ifdef TETRA      
    call gentetlinkp(vql(1, iq),input%xs%tetra%qweights)
#else
    ! Added by DIN
    write(*,*) 'Tetrahedron method for XS is disabled!'
    write(*,*) 'Check -DTETRA option in make.inc' 
    stop
#endif            
  end if

  ! Find highest (partially) occupied and lowest (partially) unoccupied states
  ! for k and k+q points 
  call findocclims(iq, istocc0, istocc, istunocc0, istunocc,&
   & isto0, isto, istu0, istu)

  ! Find limits for band combinations
  ! When coming from the df routine, i.e. screen emattype is set to 1, so o-u,u-o
  call ematbdcmbs(input%xs%emattype)

  ! Check if q-point is gamma point
  tq0 = tqgamma(iq)
  if(tq0) then
    write(unitout, '(a)') 'Info(' // trim(thisnam) // '):&
     & Gamma q - point: using momentum matrix elements for dielectric function'
  end if

  ! Write out matrix size of response function
  write(unitout, '(a, i6)') 'Info(' // thisnam // '):&
    & number of G + q vectors (local field effects):', ngq(iq)
  write(unitout, '(a, 4i6)') 'Info(' // thisnam // '):&
    & lowest (partially)  unoccupied state: ', istunocc0
  write(unitout, '(a, 4i6)') 'Info(' // thisnam // '):&
    & highest (partially) occupied state  : ', istocc0
  write(unitout, '(a, 4i5)') 'Info(' // thisnam // '):&
    & band-combination limits  nst1,  nst2,  nst3,  nst4:',  nst1,  nst2,  nst3,  nst4 
  write(unitout, '(a, 4i5)') 'Info(' // thisnam // '):&
    & band-combination limits istl1, istu1, istl2, istu2:', istl1, istu1, istl2, istu2
  write(unitout, '(a, 4i5)') 'Info(' // thisnam // '):&
    & band-combination limits istl3, istu3, istl4, istu4:', istl3, istu3, istl4, istu4

  ! Allocate arrays for eigenvalue and occupation number differences
  if(allocated(deou)) deallocate(deou)
  allocate(deou(nst1, nst2))
  if(allocated(deuo)) deallocate(deuo)
  allocate(deuo(nst3, nst4))
  if(allocated(docc12)) deallocate(docc12)
  allocate(docc12(nst1, nst2))
  if(allocated(docc21)) deallocate(docc21)
  allocate(docc21(nst3, nst4))

  ! Allocate matrix elements arrays
  if(allocated(xiou)) deallocate(xiou)
  allocate(xiou(nst1, nst2, n))
  if(allocated(xiuo)) deallocate(xiuo)
  allocate(xiuo(nst3, nst4, n))
  if(allocated(pmou)) deallocate(pmou)
  allocate(pmou(3, nst1, nst2))
  if(allocated(pmuo)) deallocate(pmuo)
  allocate(pmuo(3, nst3, nst4))

  ! Allocate arrays
  allocate(hdg(nst1, nst2, nkpt))
  allocate(scis12(nst1, nst2),scis12c(nst1, nst2))
  allocate(scis21(nst2, nst1),scis21c(nst2, nst1))
  allocate(w(nwdf))
  allocate(wreal(nwdfp))
  ! Symmetrized KS response function arrays
  !   Head for each combination of Cartesian directions 
  allocate(chi0h(3, 3, nwdfp)) 
  allocate(chi0hahc(3, 3))
  !   Wings one of the G is zero
  allocate(chi0w(n, 2, 3, nwdfp))
  !   Full chi
  allocate(chi0(n, n, nwdfp))
  allocate(wouw(nwdf), wuow(nwdf), wouh(nwdf), wuoh(nwdf))
  allocate(zvou(n), zvuo(n))
  allocate(bsedg(nst1, nst2))
  
  scis12(:, :) = 0.d0
  scis21(:, :) = 0.d0
  bsedg(:,:)=zzero

  if(input%xs%tetra%tetradf) then
#ifdef TETRA      
    allocate(cw(nwdf), cwa(nwdf), cwsurf(nwdf))

    if(input%xs%tetra%cw1k) allocate(cwt(nstsv, nstsv), &
        & cw1k(nst1, nst2, nwdfp), cwa1k(nst1, nst2, nwdfp), &
        & cwsurf1k(nst1, nst2, nwdfp))
#else
    ! Added by DIN
    write(*,*) 'Tetrahedron method for XS is disabled!'
    write(*,*) 'Check -DTETRA option in make.inc' 
    stop
#endif            
  end if

  ! Generate complex energy grid
  wintv(1)=input%xs%energywindow%intv(1)
  wintv(2)=input%xs%energywindow%intv(2)

  ! For calculation of static screening the first frequency point should be zero
  if(task.eq.430) wintv(1)=0.d0
  call genwgrid(nwdf, wintv, input%xs%tddft%acont, 0.d0, w_cmplx=w)

  wreal(:) = dble(w(wi:wf))
  if(wreal(1) .lt. epstetra) wreal(1) = epstetra

  ! Initializations
  chi0(:, :, :) = zzero
  chi0w(:, :, :, :) = zzero
  chi0h(:, :, :) = zzero
  chi0hahc(:, :) = zzero

  if(tscreen) then
    ! Generate radial integrals wrt. sph. bessel functions
    call timesec(tc)
    call ematrad(iq)
    call timesec(td)

    ! Delete timing information of previous runs
    call filedel(trim(fnetim))

    ! Write information
    write(unitout, '(a, i6)') 'Info(' // thisnam // '):&
      & number of G + q vectors:', ngq(iq)
    call ematqalloc
  end if

  if(tfxcbse) then
    call getbsediag
    write(unitout, '("Info(", a, "): read diagonal of BSE kernel")') trim(thisnam)
    write(unitout, '(" mean value : ", 2g18.10)') bsed
    bsedg(:, :) = bsed
  end if

  call timesec(tb)

  ! Loop over k-points
  kloop: do ik = 1, nkpt

    ! K-point analysis
    ! Check whether the transitions at k-point ik should be included.
    if( .not. transik(ik)) cycle
    call chkpt(3, (/ task, iq, ik /), 'dfq: task, q-point index, k-point index')

    cpuosc = 0.d0
    cpuupd = 0.d0
    call timesec(cpu0)

    ikq = ikmapikq(ik, iq)

    ! Get second variational KS transition energies, scissor shifts 
    ! and occupancy differences for current k+q/k point combination 
    ! and the specified band ranges.
    call getdevaldoccsv(iq, ik, ikq, istl1, istu1, istl2, istu2,&
      & deou, docc12, scis12)
    call getdevaldoccsv(iq, ik, ikq, istl2, istu2, istl1, istu1,&
      & deuo, docc21, scis21)

    ! Create copies of scissor shifts
    scis12c(:,:)=scis12(:,:)
    scis21c(:,:)=scis21(:,:)

    if(tscreen) then
      ! Do not use scissors correction for screening
      if(task .eq. 430) then
        scis12c(:, :) = zzero
        scis21c(:, :) = zzero
      end if
      ! For screening calculate matrix elements of plane wave on the fly.
      ! The plane wave elements for occupied unoccupied transitions are 
      ! calculated and stored in xiou
!! xiuo is not calculated at the moment
      call ematqk1(iq, ik)
      ! Allocate anti-resonant plane wave matrix elements
      if( .not. allocated(xiuo)) allocate(xiuo(nst3, nst4, n))
      ! Allocate anti-resonant momentum matrix elements
      if( .not. allocated(pmuo)) allocate(pmuo(3, nst3, nst4))
    end if

    ! Add bse diagonal shift use with bse-kernel
    scis12c(:, :) = scis12c(:, :) + bsedg(:, :)
    scis21c(:, :) = scis21c(:, :) + transpose(bsedg(:, :))

    ! Get matrix elements (exp. expr. or momentum op.)
!! xiuo is used by getpemat
    ! Get m12=v^{1/2}*M_ou, m34=v^{1/2}*M_uo 
    !     p12=-Sqrt{4pi}P12/dE12, 
    !     p34=-Sqrt{4pi}P34/dE34
    call getpemat(iq, ik, trim(fnpmat), trim(fnemat),&
      & m12=xiou, m34=xiuo, p12=pmou, p34=pmuo)
      
    ! Set matrix elements to one for Lindhard function
    if(input%xs%tddft%lindhard) then
      ! Set g=0 components to one
      xiou(:, :, 1) = zone
      xiuo(:, :, 1) = zone
      ! Set g/=0 components to zero
      if(n .gt. 1) then
        xiou(:, :, 2:) = zzero
        xiuo(:, :, 2:) = zzero
      end if
      ! Set momentum matrix elements to one
      pmou(:, :, :) = zone
      pmuo(:, :, :) = zone
    end if

    if(input%xs%tetra%cw1k) then
#ifdef TETRA          
      do iw = 1, nwdfp
        call tetcwifc_1k(ik, nkpt, nstsv, evalsv, efermi,&
          & wreal(iw), 2, cwt)
        cw1k(:, :, iw) = cwt(istl1:istu1, istl2:istu2)
        call tetcwifc_1k(ik, nkpt, nstsv, evalsv,&
          & efermi,-wreal(iw), 2, cwt)
        cwa1k(:, :, iw) = cwt(istl1:istu1, istl2:istu2)
        call tetcwifc_1k(ik, nkpt, nstsv, evalsv, efermi,&
          & wreal(iw), 4, cwt)
        cwsurf1k(:, :, iw) = cwt(istl1:istu1, istl2:istu2)
      end do
#else
      ! Added by DIN
      write(*,*) 'Tetrahedron method for XS is disabled!'
      write(*,*) 'Check -DTETRA option in make.inc' 
      stop
#endif            
    end if

    if(tscreen) then
      ! We don't need anti-resonant parts here, assign them the same
      ! Value as for resonant parts, resulting in a factor of two.
      do igq = 1, n
        xiuo(:, :, igq) = transpose(xiou(:, :, igq))
      end do
      do j = 1, 3
        pmuo(j, :, :) = transpose(pmou(j, :, :))
      end do
      deuo(:, :) = transpose(deou(:, :))
      docc21 (:, :) = transpose(docc12(:, :))
      scis21c(:, :) = transpose(scis12c(:, :))
    end if
    ! Turn off anti-resonant terms (type 2-1 band combinations) for Kohn-Sham
    ! response function
    if(( .not. input%xs%tddft%aresdf) .and. ( .not. tscreen)) then
      xiuo(:, :, :) = zzero
      pmuo(:, :, :) = zzero
    end if

    do ist1 = 1, istocc0 - istunocc0 + 1
      do ist2 = 1, istocc0 - istunocc0 + 1
        j = ist1 + istunocc0 - 1
        ! Set lower triangle of first block to zero
        if(ist1 .gt. ist2) then
          xiou(j, ist2, :) = zzero
          pmou(:, j, ist2) = zzero
        end if
        ! Set diagonal to zero (project out intraband contributions)
        if(( .not. input%xs%tddft%intraband) .and. (ist1 .eq. ist2)) then
          xiou(j, ist2, :) = zzero
          pmou(:, j, ist2) = zzero
        end if
        ! Set upper triangle of second block to zero
        ! also set diagonal to zero to avoid double counting
        if(ist1 .ge. ist2) then
          xiuo(ist2, j, :) = zzero
          pmuo(:, ist2, j) = zzero
        end if
      end do
    end do

    call timesec(cpu1)
    cpuread = cpu1 - cpu0

    allocate(wou(nwdf,nst1,nst2))
    allocate(wuo(nwdf,nst1,nst2))

    ist1loop: do ist1 = 1, nst1
      ist2loop: do ist2 = 1, nst2
        !---------------------!
        !     denominator     !
        !---------------------!
        ! Absolute band indices
        i1 = ist1
        i2 = istunocc0 + ist2 - 1
        ! Band analysis
        ! Check if the considered transition is included.
        if( .not. transijst(ik, i1, i2)) cycle

        call timesec(cpu0)

        if(input%xs%tetra%tetradf) then
#ifdef TETRA               
          ! Mirror index pair on diagonal if necessary
          if(i1 .gt. i2) then
            j1 = ist2
            j2 = ist1 - istunocc0 + 1
          else
            j1 = ist1
            j2 = ist2
          end if
          ! Read weights for tetrahedron method
          if(input%xs%tetra%cw1k) then
            cw(wi:wf) = cw1k(ist1, ist2, :)
            cwa(wi:wf) = cwa1k(ist1, ist2, :)
            cwsurf(wi:wf) = cwsurf1k(ist1, ist2, :)
          else
            call gettetcw(iq, ik, j1, j2, nst1, nst2, nwdf,&
              & trim(fnwtet), cw, cwa, cwsurf)
          end if
          ! Include occupation number differences
          wou(wi:wf,ist1,ist2) = docc12 (ist1, ist2) * cmplx(cw(wi:wf),&
            & cwsurf(wi:wf), 8) / omega
          wuo(wi:wf,ist1,ist2) = - docc21 (ist2, ist1) *&
            & cmplx(cwa(wi:wf), 0.d0, 8) / omega

          if(tq0) then
            ! Rescale: use delta-function delta(e_nmk + scis_nmk - w)
            ! Take real part of bse diagonal(being contained in scis12c)
            ! Since tetrahedron method formalism implemented does not allow
            ! Otherwise
            wouw(wi:wf) = cmplx(dble(wou(wi:wf,ist1,ist2)),&
              & aimag(wou(wi:wf,ist1,ist2))*deou(ist1, ist2)/(-wreal(:)-dble(scis12c(ist1, ist2))))
            wuow(wi:wf) = cmplx(dble(wuo(wi:wf,ist1,ist2)),&
              & aimag(wuo(wi:wf,ist1,ist2))*deuo(ist2, ist1)/(-wreal(:)-dble(scis21c(ist2, ist1))))
            wouh(wi:wf) = cmplx(dble(wou(wi:wf,ist1,ist2)),&
              & aimag(wou(wi:wf,ist1,ist2))*deou(ist1, ist2)**2/(-wreal(:)-dble(scis12c(ist1, ist2)))**2)
            wuoh(wi:wf) = cmplx(dble(wuo(wi:wf,ist1,ist2)),&
              & aimag(wuo(wi:wf,ist1,ist2))*deuo(ist2, ist1)**2/(-wreal(:)-dble(scis21c(ist2, ist1)))**2)
          end if
#else
          ! Added by DIN
          write(*,*) 'Tetrahedron method for XS is disabled!'
          write(*,*) 'Check -DTETRA option in make.inc' 
          stop
#endif            
        else
          ! Include occupation number differences
          do iw=wi,wf
            ! Check for vanishing denominators in case of screening
            ! (no broadening)
            zt1=w(iw)+deou(ist1, ist2)+scis12c(ist1,ist2)+zi*brd
            if(abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
            wou(iw,ist1,ist2) = docc12(ist1, ist2) * wkpt(ik) / omega / zt1
            zt1=w(iw)+deuo(ist2, ist1)+scis21c(ist2,ist1)+tordf*zi*brd
            if(abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
            wuo(iw,ist1,ist2) = docc21(ist2, ist1) * wkpt(ik) / omega / zt1
          end do

          wouw(wi:wf) = wou(wi:wf,ist1,ist2)
          wuow(wi:wf) = wuo(wi:wf,ist1,ist2)
          wouh(wi:wf) = wou(wi:wf,ist1,ist2)
          wuoh(wi:wf) = wuo(wi:wf,ist1,ist2)
        end if

        call timesec(cpu1)
        cpuosc=cpuosc+cpu1-cpu0

        !----------------------------------!
        !     Update response function     !
        !----------------------------------!
        zvou(:) = xiou(ist1, ist2, :)
        zvuo(:) = xiuo(ist2, ist1, :)

        do iw = wi, wf

          G0: if(tq0) then

            do oct1 = 1, 3
              ! Wings
              chi0w(2:, 1, oct1, iw-wi+1) = chi0w(2:, 1, oct1, iw-wi+1)&
                &+ wouw(iw) * pmou(oct1, ist1, ist2) * conjg(zvou(2:))&
                &+ wuow(iw) * pmuo(oct1, ist2, ist1) * conjg(zvuo(2:))
              chi0w(2:, 2, oct1, iw-wi+1) = chi0w(2:, 2, oct1, iw-wi+1)&
                &+ wouw(iw) * zvou(2:) * conjg(pmou(oct1, ist1, ist2))&
                &+ wuow(iw) * zvuo(2:) * conjg(pmuo(oct1, ist2, ist1))

              do oct2 = 1, 3
                ! Head
                if(.not.input%xs%tddft%ahc) then
                  chi0h(oct1, oct2, iw-wi+1) = chi0h(oct1, oct2, iw-wi+1)&
                    &+ wouh(iw) * pmou(oct1, ist1, ist2) * conjg(pmou(oct2, ist1, ist2))&
                    &+ wuoh(iw) * pmuo(oct1, ist2, ist1) * conjg(pmuo(oct2, ist2, ist1))
                else
                  winv=1.0d0/(w(iw)+zi*brd)
                  if(abs(w(iw)).lt.1.d-8) winv=1.d0
                  chi0h(oct1, oct2, iw-wi+1) = chi0h(oct1, oct2, iw-wi+1)&
                    &+ wouh(iw) * pmou(oct1, ist1, ist2) * conjg(pmou(oct2, ist1, ist2))&
                    &* deou(ist1, ist2) * winv + wuoh(iw) * pmuo(oct1, ist2, ist1)&
                    &* conjg(pmuo(oct2, ist2, ist1)) * deuo(ist2, ist1) * winv 
                end if

              end do

            end do

          end if G0

        end do

        call timesec(cpu0)
        cpuupd = cpuupd + cpu0 - cpu1

      ! End loop over states combinations
      end do ist2loop
    end do ist1loop

    allocate(zm(n,nst1,nst2))
    do iw=wi,wf
    
      do ist2 = 1, nst2
        do ist1 = 1, nst1
          zm(:,ist1,ist2)=conjg(wou(iw,ist1,ist2)*xiou(ist1,ist2,:))
        end do
      end do

      call zgemm('n', 'n', n, n, nst1*nst2, zone,&
        & zm(1,1,1), n, xiou(1,1,1), nst1*nst2, zone, chi0(1,1,iw-wi+1), n)
    end do

    deallocate(zm)
    allocate(zm(n,nst2,nst1))

    do iw=wi,wf

      do ist2 = 1, nst2
        do ist1 = 1, nst1
          zm(:,ist2,ist1)=conjg(wuo(iw,ist1,ist2)*xiuo(ist2,ist1,:))
        end do
      end do
           
      call zgemm('n', 'n', n, n, nst1*nst2, zone, zm(1,1,1),&
        & n, xiuo(1,1,1), nst1*nst2, zone, chi0(1,1,iw-wi+1), n)

    end do

    call timesec(cpu1)
    cpuupd = cpuupd + cpu1 - cpu0

    deallocate(wou,wuo,zm)

!*****************************************************************************************************

    cputot = cpuread + cpuosc + cpuupd
    ! Timing information
    call dftim(iq, ik, trim(fnxtim), cpuread, cpuosc, cpuupd, cputot)

    ! Synchronize
    if( .not. tscreen) call barrier

  ! End loop over k-points
  end do kloop

  chi0(:,:,:)=conjg(chi0(:,:,:))

  wplas = input%xs%tddft%drude(1)
  wrel = input%xs%tddft%drude(2)
  if((wplas>1.d-8).and.(wrel>1.d-8)) then
    do iw = wi, wf
      winv=1.0d0/(w(iw)+zi*brd)
      if(abs(w(iw)).lt.1.d-8) winv=1.d0
      do oct1 = 1, 3
        chi0h(oct1, oct1, iw-wi+1) = chi0h(oct1, oct1, iw-wi+1)&
          &+ wplas**2/(w(iw)+zi*wrel)*winv
      end do
    end do
  end if

  if(tscreen) call ematqdealloc

  ! Symmetrize head
  head: if(tq0) then
    allocate(chi0hs(3, 3, nwdfp), eps0(3, 3, nwdf))

    ! Write dielectric tensor to file (unsymmetrized)
    forall(iw=1:nwdf)
      eps0(:, :, iw) = dble(krondelta) - chi0h(:, :, iw)
    end forall

    if(rank .eq. 0)&
      & call writedielt('DIELTENS0_NOSYM', 1, 0.d0, eps0(:, :, 1), 0)

    ! Symmetrize the macroscopic dielectric function tensor
    do oct1 = 1, 3
      do oct2 = 1, 3
        call symt2app(oct1, oct2, nwdfp, symt2, chi0h, chi0hs(oct1, oct2, :))
      end do
    end do

    ! Re-assign the symmetrized head
    chi0h(:, :, :) = chi0hs(:, :, :)

    ! Write dielectric tensor to file
    forall(iw=1:nwdf)
      eps0(:, :, iw) = dble(krondelta) - chi0hs(:, :, iw)
    end forall
    if(rank .eq. 0)&
     & call writedielt('DIELTENS0', 1, 0.d0, eps0(:, :, 1), 0)

    deallocate(chi0hs, eps0)

  end if head 

  ! Write response function to file
  if(tscreen) then
    ! Write out screening
    call getunit(un)
    open(un, file=trim(fnscreen), form='formatted', action='write', status='replace')
    call putscreen(un, tq0, n, chi0(:, :, 1), chi0h(:, :, 1), chi0w(:, :, :, 1))
    call writevars(un, iq, 0)
    close(un)
  else
    do j = 0, procs - 1
      if(rank .eq. j) then
        do iw = wi, wf
          call putx0(tq0, iq, iw-wi+1, trim(fnchi0_t), '',&
            & chi0(:, :, iw-wi+1), chi0w(:, :, :, iw-wi+1), chi0h(:, :, iw-wi+1))
        end do
      end if
      call barrier
    end do
  end if

  deallocate(chi0, chi0h, chi0w)
  deallocate(docc12, docc21, scis12, scis21, scis12c, scis21c)
  deallocate(deou, deuo, wouw, wuow, wouh, wuoh, zvou, zvuo)
  deallocate(xiou, xiuo, pmou, pmuo)
  deallocate(bsedg)
  deallocate(w, wreal)

  if(input%xs%tetra%tetradf) then
#ifdef TETRA      
    deallocate(cw, cwa, cwsurf)
    if(input%xs%tetra%cw1k) deallocate(cwt, cw1k, cwa1k, cwsurf1k)
#endif        
  end if

end subroutine dfq
!EOC
