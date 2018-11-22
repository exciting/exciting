! Copyright (C) 2005-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dfq
! !INTERFACE:
subroutine dfq(iq)
! !USES:
  use mod_misc, only: filext
  use modinput, only: input
  use modmpi, only: procs, rank, barrier
  use mod_misc, only: task
  use mod_constants, only: zzero, zone, zi, krondelta
  use mod_kpoint, only: nkpt, wkpt
  use mod_qpoint, only: nqpt, vql
  use mod_lattice, only: omega
  use modxs, only: tfxcbse, tscreen, bzsampl, wpari,&
                & wparf, ngq, fnpmat,&
                & fnetim, fnxtim, fnemat, fnchi0,&
                & fnchi0_t, qvkloff, istocc0, istocc,&
                & istunocc0, istunocc, isto0, isto,&
                & istu0, istu, unitout, nst1,&
                & nst2, nst3, nst4, istl1,&
                & istu1, istl2, istu2, istl3,&
                & istu3, istl4, istu4, deou,&
                & deuo, docc12, docc21, xiou,&
                & xiuo, pmou, pmuo, nwdf,&
                & bsed, ikmapikq, tordf, symt2,&
                & bcbs, filexteps,&
                & eps0dirname, scrdirname, timingdirname
#ifdef TETRA      
  use modxs, only: fnwtet
  use mod_eigenvalue_occupancy, only: nstsv, evalsv, efermi 
  use modtetra
#endif
  use m_genwgrid
  use m_getpemat
  use m_dftim
  use m_gettetcw
  use m_putx0
  use m_getunit
  use m_writevars
  use m_filedel
  use m_genfilname
  use m_ematqk
  use m_writecmplxparts
  use m_putgeteps0
  use mod_variation, only: ematqk_sv

! !INPUT/OUTPUT PARAMETERS:
! In:
! integer(4) :: iq  ! q-point index
!
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
!   as ${\bf q}\rightarrow 0$ along the Cartesian basis vectors ${\bf e_i}$,
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
!   Changed parts that are unique to execution of dfq with tscreen = true (Aurich)
!     Changed Plane wave matrix elements construction call
!     Changed write out of EPS
!     Added possibility to fully calculate anti-resonant part
!EOP
!BOC

  implicit none

  ! Arguments
  integer(4), intent(in) :: iq

  ! Local variables
  character(*), parameter :: thisname = 'dfq'
  character(256) :: fnscreen, fneps0, filex
  complex(8) :: zt1, winv
  complex(8), allocatable :: w(:)
  complex(8), allocatable :: chi0(:, :, :)
  complex(8), allocatable :: chi0w(:, :, :, :), chi0h(:, :, :), eps0(:, :, :)
  complex(8), allocatable :: chi0hahc(:, :)
  complex(8), allocatable :: wou(:,:,:), wuo(:,:,:), wouw(:), wuow(:), wouh(:), wuoh(:)
  complex(8), allocatable :: zvou(:), zvuo(:), chi0hs(:, :, :), bsedg(:, :)
  complex(8), allocatable :: scis12c(:, :), scis21c(:, :), zm(:,:,:)
  real(8), parameter :: epstetra = 1.d-8
  real(8) :: ta,tb,tc,td
  real(8), allocatable :: wreal(:) 
  real(8), allocatable :: scis12(:, :), scis21(:, :)
  real(8) :: brd, cpu0, cpu1, cpuread, cpuosc, cpuupd, cputot, wintv(2), wplas, wrel
  integer(4) :: n, j, i1, i2, ik, ikq, igq, iw, wi, wf, ist1, ist2, nwdfp
  integer(4) :: numpo
  integer(4) :: oct1, oct2, un
  logical :: tq0
  logical :: doares, fintraband
  type(bcbs) :: bc

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

  ! True if exchange-correlation kernel "MB1" or "BO" is use and
  ! it is not a BSE screening task.
  tfxcbse = ((input%xs%tddft%fxctypenumber .eq. 7) .or.&
    & (input%xs%tddft%fxctypenumber .eq. 8)) .and. ( .not. tscreen)

  ! Sampling of Brillouin zone (Tetra)
  bzsampl = 0
#ifdef TETRA
  if(input%xs%tetra%tetradf) bzsampl = 1
#endif

  ! Initial and final frequency-point (parallelization)
  wi = wpari
  wf = wparf
  nwdfp = wf - wi + 1

  ! Matrix size for response function
  ! Get number of G+q vectors for current q
  n = ngq(iq)

  ! Zero broadening for analytic continuation
  brd = input%xs%broad
  ! acont defaults to false
  if(input%xs%tddft%acont) brd = 0.d0

!! DISCUSS
  doares = .true.
  ! Task 430 is 'screen'
  if(task .eq. 430) then 
    ! Zero broadening for dielectric matrix (w=0) for band-gap systems
    brd = 0.d0
    ! Calculate anti-resonant part of chi explicitly ?
    doares = .not. input%xs%screening%tr
  end if

  ! Set whether intraband should be used
  fintraband=input%xs%tddft%intraband .or. input%xs%screening%intraband
  ! File extension for q-point (not in 'screen')
  if( .not. tscreen) call genfilname(iqmt=iq, setfilext=.true.)

  ! Filenames for output
  if(tscreen) then
#ifdef TETRA               
    call genfilname(basename='TETW', iq=iq, appfilext=.true., filnam=fnwtet)
#endif               
    call genfilname(basename='PMAT', appfilext=.true., filnam=fnpmat)
    filex=filext
    filext=filexteps
    call genfilname(basename=trim(adjustl(scrdirname))//'/'//'SCREEN', appfilext=.true., iq=iq, filnam=fnscreen)
    call genfilname(basename=trim(adjustl(eps0dirname))//'/'//'EPS0', appfilext=.true., iq=iq, filnam=fneps0)
    filext=filex

    filex=filext
    filext=filexteps
    call genfilname(nodotpar=.true., basename='EMAT_TIMING', iq=iq,&
     & etype=input%xs%emattype, procs=procs, rank=rank, appfilext=.true., filnam=fnetim)
    fnetim=trim(adjustl(timingdirname))//'/'//trim(adjustl(fnetim))
    call genfilname(nodotpar=.true., basename='X0_TIMING', iq=iq,&
     & procs=procs, rank=rank, appfilext=.true., filnam=fnxtim)
    fnxtim=trim(adjustl(timingdirname))//'/'//trim(adjustl(fnxtim))
    filext=filex

  else
#ifdef TETRA               
    call genfilname(basename='TETW', iqmt=iq, filnam=fnwtet)
#endif               
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

  ! Calculate k+q and G+k+q related variables
  ! by setting the offset generated from vkloff and the q point
  ! and then calling init1
  call init1offs(qvkloff(1:3, iq))

  !write(*,*) "df: qvkloff=", qvkloff(1:3,iq)

  ! TETRA: Generate link array for tetrahedra
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

  !write(*,*) "dfq: iq=", iq, " filext=", trim(filext)

  ! Find highest (partially) occupied and lowest (partially) unoccupied states
  ! for k and k+q points 
  !write(*,*) "dfq: ikmapikq(:,iq)", ikmapikq(:,iq)
  
  call findocclims(iq, ikmapikq(:,iq), istocc0, istunocc0, isto0, isto, istu0, istu)
  istunocc = istunocc0
  istocc = istocc0

  ! Find limits for band combinations
  !   When coming from the df routine, i.e. screen emattype is set to 1, so o-u,u-o
  call ematbdcmbs(input%xs%emattype)

  tq0 = tqgamma(iq)
  ! Check if q-point is gamma point (uses mod_qpoint::vqc)
  if (input%xs%BSE%outputlevelnumber == 1) then
    write(unitout, *)
    write(unitout, '(a,i4)') 'Info(' // trim(thisname) // '):&
    & Calculating screening for q-point: ', iq 
    if(tq0) then
      write(unitout, '(a)') 'Info(' // trim(thisname) // '):&
       & Gamma q-point: using momentum matrix elements for dielectric function'
    end if
    ! Write out matrix size of response function and contributing bands
    write(unitout, '(a, i6)') 'Info(' // thisname // '):&
      & number of G + q vectors (local field effects):', ngq(iq)
    write(unitout, '(a, 4i6)') 'Info(' // thisname // '):&
      & lowest (partially)  unoccupied state: ', istunocc0
    write(unitout, '(a, 4i6)') 'Info(' // thisname // '):&
      & highest (partially) occupied state  : ', istocc0
    write(unitout, '(a, 4i5)') 'Info(' // thisname // '):&
      & band-combination limits  nst1,  nst2,  nst3,  nst4:',  nst1,  nst2,  nst3,  nst4 
    write(unitout, '(a, 4i5)') 'Info(' // thisname // '):&
      & band-combination limits istl1, istu1, istl2, istu2:', istl1, istu1, istl2, istu2
    write(unitout, '(a, 4i5)') 'Info(' // thisname // '):&
      & band-combination limits istl3, istu3, istl4, istu4:', istl3, istu3, istl4, istu4
  end if
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
  allocate(scis12(nst1, nst2),scis12c(nst1, nst2))
  allocate(scis21(nst2, nst1),scis21c(nst2, nst1))
  allocate(w(nwdf))
  allocate(wreal(nwdfp))

  ! Symmetrized KS response function arrays
  !   Head for each combination of Cartesian directions 
  allocate(chi0h(3, 3, nwdfp)) 
  allocate(chi0hahc(3, 3))
  !   The 2 wings for 3 Cartesian directions
  allocate(chi0w(n, 2, 3, nwdfp))
  !   Full chi
  allocate(chi0(n, n, nwdfp))
  !   Weights
  allocate(wouw(nwdf), wuow(nwdf), wouh(nwdf), wuoh(nwdf))

  ! v^1/2*M_ou and v^1/2*M_uo
  allocate(zvou(n), zvuo(n))

  allocate(bsedg(nst1, nst2))
  
  ! Zeroing 
  scis12(:, :) = 0.d0
  scis21(:, :) = 0.d0
  bsedg(:,:)=zzero

  ! TETRA
  if(input%xs%tetra%tetradf) then
#ifdef TETRA      
    allocate(cw(nwdf), cwa(nwdf), cwsurf(nwdf))

    if(input%xs%tetra%cw1k) allocate(cwt(nstsv, nstsv),&
      & cw1k(nst1, nst2, nwdfp), cwa1k(nst1, nst2, nwdfp),&
      & cwsurf1k(nst1, nst2, nwdfp))
#else
    ! Added by DIN
    write(*,*) 'Tetrahedron method for XS is disabled!'
    write(*,*) 'Check -DTETRA option in make.inc' 
    stop
#endif            
  end if

  !! Generate complex energy grid
  ! Intervall limits
  wintv(1)=input%xs%energywindow%intv(1)
  wintv(2)=input%xs%energywindow%intv(2)
  ! For calculation of static screening the first frequency point should be zero
  if(task .eq. 430) wintv(1)=0.d0
  ! Make evenly spaced complex valued omega grid.
  !   acont defaults to false, broadening set to zero --> complex w grid = real w grid
  call genwgrid(nwdf, wintv, input%xs%tddft%acont, 0.d0, w_cmplx=w)

  ! Real frequency grid
  wreal(:) = dble(w(wi:wf))

#ifdef TETRA
  ! Set first real frequency to 10^{-8}
  if(wreal(1) .lt. epstetra) wreal(1) = epstetra
#endif

  ! Zeroing chi arrays
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
    if (input%xs%BSE%outputlevelnumber == 1) &
      & write(unitout, '(a, i6)') 'Info(' // thisname // '):&
      & number of G + q vectors:', ngq(iq)
    call ematqalloc
  end if

  ! Read BSE diagonal for BSE TDDFT kernel
  if(tfxcbse) then
    call getbsediag
    write(unitout, '("Info(", a, "): read diagonal of BSE kernel")') trim(thisname)
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

    ! Timing 
    cpuosc = 0.d0
    cpuupd = 0.d0
    call timesec(cpu0)

    ikq = ikmapikq(ik, iq)
    !write(*,*) "dfq: ik, iq, ikq", ik, iq, ikq

    ! Get second variational KS transition energies, scissor shifts 
    ! and occupancy differences for current k+q/k point combination 
    ! and the specified band ranges.
    ! ou
    !   deou = e_{ok} - e_{uk+q}
    !   docc12 = f_{ok} - f_{uk+q}
    call getdevaldoccsv(iq, ik, ikq, istl1, istu1, istl2, istu2,&
      & deou, docc12, scis12)

    if(doares) then
      ! uo
      !   deuo = e_{uk} - e_{ok+q}
      !   docc21 = f_{uk} - f_{ok+q}
      call getdevaldoccsv(iq, ik, ikq, istl2, istu2, istl1, istu1,&
        & deuo, docc21, scis21)
    end if

    ! Create copies of scissor shifts
    scis12c(:,:)=scis12(:,:)
    scis21c(:,:)=scis21(:,:)

    if(tscreen) then

      ! Do not use scissors correction for screening (why?)
      if(task .eq. 430) then
        scis12c(:, :) = zzero
        scis21c(:, :) = zzero
      end if

!*********** Plane wave calculation *****************************!
      ! For screening calculate matrix elements of plane wave on the fly.

      ! The plane wave elements for ou and uo transitions are 
      ! calculated and stored in xiou and xiuo
      ! Set 12=ou 34=uo
      call ematbdcmbs(1)
      ! Get ou
      if(allocated(xiou)) deallocate(xiou)
      allocate(xiou(nst1, nst2, n))
      bc%n1 = nst1
      bc%n2 = nst2
      bc%il1 = istl1
      bc%il2 = istl2
      bc%iu1 = istu1
      bc%iu2 = istu2
      ikmapikq_ptr => ikmapikq
      call setptr01
      if (.not. input%groundstate%tevecsv) then
        call ematqk(iq, ik, xiou, bc)
      else
        call ematqk_sv(iq, ik, xiou, bc)
      end if
      ! Get uo
      if(allocated(xiuo)) deallocate(xiuo)
      allocate(xiuo(nst3, nst4, n))

      ! Note:
      ! The uo plane wave matrix elements are not needed
      ! if the time reversal symmetry is applied to the 
      ! anti-resonant part of Chi0. 

      ! t.r. sym not used
      if( doares ) then
        bc%n1 = nst3
        bc%n2 = nst4
        bc%il1 = istl3
        bc%il2 = istl4
        bc%iu1 = istu3
        bc%iu2 = istu4
        ikmapikq_ptr => ikmapikq
        call setptr01
        if (.not. input%groundstate%tevecsv) then
          call ematqk(iq, ik, xiuo, bc)
        else
          call ematqk_sv(iq, ik, xiuo, bc)
        end if
      end if

!********************************************************!
    end if

    ! Add BSE diagonal shift use with bse-kernel
    ! For 'screen' this is all zero
    scis12c(:, :) = scis12c(:, :) + bsedg(:, :)
    scis21c(:, :) = scis21c(:, :) + transpose(bsedg(:, :))

    ! Get plane wave and momentum matrix elements 
    ! multiplied by the square root of the Coulomb potential.
    ! Get m12=v^{1/2}*M_12, m34=v^{1/2}*M_34 
    !     p12=-Sqrt{4pi}P12/dE12, 
    !     p34=-Sqrt{4pi}P34/dE34
    ! Momentum matrix elements are always read from file.
    ! NOTE: If tscreen, then xiou and xiou need to be already present 
    ! in modxs. And they are used implicitly as intent inout.

    if(doares) then
      call getpemat(iq, ik, trim(fnpmat), trim(fnemat),&
        & m12=xiou, m34=xiuo, p12=pmou, p34=pmuo)
    else
      call getpemat(iq, ik, trim(fnpmat), trim(fnemat),&
        & m12=xiou, p12=pmou)
    end if

    ! Set matrix elements to one for Lindhard function (default is false)
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

    ! TETRA
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
    ! Not screen: Turn off anti-resonant terms (type 2-1 band combinations) for Kohn-Sham
    ! response function (default skip if)
    if(( .not. input%xs%tddft%aresdf) .and. ( .not. tscreen)) then
      xiuo(:, :, :) = zzero
      pmuo(:, :, :) = zzero
    end if

    ! Treatment of partially occupied bands.
    !
    ! Calculate number of partially occupied bands
    !   istocc0: Highest at least partially occupied band
    !   istunocc0: Lowest at least partially occupied band
    numpo = istocc0 - istunocc0 + 1
    ! In the presence of partially occupied bands the plane wave matrix
    ! M_ou/M_uo has the following block form 
    ! where: (o)ccupied, (u)noccupied,
    !        (p)artially(o)ccupied, (p)artially(u)noccupied
    !          ______________            ______________
    !         | o/pu  | o/u  |          | pu/o | pu/po |
    ! M_ou =  |--------------|  M_uo =  |--------------|
    !         | po/pu | po/u |          | u/po |  u/po |
    !         |______________|          |______________|
    ! 
    !  The following loops are concerned with the po/pu (pu/po) part.
    !  Since we split M_{nm}M^*_{nm} into M_{ou}M^*_{ou}+M_{uo}M^*_{uo}
    !  we have to take care about double counting transitions if m AND n
    !  are referring to partially filled bands.
    !  
    do ist1 = 1, numpo
      do ist2 = 1, numpo
        ! Get band index of occupied state (counted from lowest energy state)
        j = ist1 + istunocc0 - 1
        ! Set lower triangle of po/pu block to zero, i.e. allow only
        ! transitions from energetically higher to lower bands.
        ! Exclude diagonal for momentum matrix, since weights should be 
        ! zero due to f_nk-f_nk (q is always 0 for momentum matrix).

        ! >= (p is used for q=0, and no intraband is possible)
        if(ist1 .ge. ist2) then
          pmou(:, j, ist2) = zzero
        end if
        if(ist1 .gt. ist2) then
          xiou(j, ist2, :) = zzero
        end if
        if(.not. doares) then 
          ! Account for double counting of interband contributions
          ! when Chi0 = 2*Chi0^Res is used.
          if(ist1 .eq. ist2) then
            xiou(j, ist2, :) = xiou(j, ist2, :)*0.5d0
          end if
        end if
        ! Set diagonal to zero (project out intra-band contributions)?
        ! Intraband default changed to true
        if( ist1 .eq. ist2) then 
          if(.not. fintraband) then
            xiou(j, ist2, :) = zzero
          end if
        end if
        if(doares) then
          ! Set upper triangle of pu/po block to zero, i.e. allow only
          ! transitions from energetically lower to higher bands.
          ! The diagonal needs to be excluded since those transitions
          ! are already included in M_ou.
          if(ist1 .ge. ist2) then
            xiuo(ist2, j, :) = zzero
            pmuo(:, ist2, j) = zzero
          end if
        end if
      ! End loop over partial occupied bands
      end do
    end do

    call timesec(cpu1)
    cpuread = cpu1 - cpu0

    ! Allocate resonant and anti-resonant weights
    allocate(wou(nwdf,nst1,nst2))
    allocate(wuo(nwdf,nst1,nst2))

    ! Occupied 
    ist1loop: do ist1 = 1, nst1
      ! Unoccupied
      ist2loop: do ist2 = 1, nst2
        !---------------------!
        !     Denominator     !
        !---------------------!
        ! Absolute band indices
        i1 = ist1
        i2 = istunocc0 + ist2 - 1
        ! Band analysis
        ! Check if the considered transition is included.
        if( .not. transijst(ik, i1, i2)) cycle

        call timesec(cpu0)

        ! TETRA
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
          wou(wi:wf,ist1,ist2) = docc12(ist1, ist2) * cmplx(cw(wi:wf),&
            & cwsurf(wi:wf), 8) / omega
          wuo(wi:wf,ist1,ist2) = - docc21(ist2, ist1) *&
            & cmplx(cwa(wi:wf), 0.d0, 8) / omega

          if(tq0) then
            ! Rescale: use delta-function delta(e_nmk + scis_nmk - w)
            ! Take real part of bse diagonal(being contained in scis12c)
            ! Since tetrahedron method formalism implemented does not allow
            ! Otherwise
            wouw(wi:wf) = cmplx(dble(wou(wi:wf,ist1,ist2)),&
              & aimag(wou(wi:wf,ist1,ist2))*deou(ist1, ist2)&
              &/(-wreal(:)-dble(scis12c(ist1, ist2))))
            wuow(wi:wf) = cmplx(dble(wuo(wi:wf,ist1,ist2)),&
              & aimag(wuo(wi:wf,ist1,ist2))*deuo(ist2, ist1)&
              &/(-wreal(:)-dble(scis21c(ist2, ist1))))
            wouh(wi:wf) = cmplx(dble(wou(wi:wf,ist1,ist2)),&
              & aimag(wou(wi:wf,ist1,ist2))*deou(ist1, ist2)**2&
              &/(-wreal(:)-dble(scis12c(ist1, ist2)))**2)
            wuoh(wi:wf) = cmplx(dble(wuo(wi:wf,ist1,ist2)),&
              & aimag(wuo(wi:wf,ist1,ist2))*deuo(ist2, ist1)**2&
              &/(-wreal(:)-dble(scis21c(ist2, ist1)))**2)
          end if
#else
          ! Added by DIN
          write(*,*) 'Tetrahedron method for XS is disabled!'
          write(*,*) 'Check -DTETRA option in make.inc' 
          stop
#endif            
        else

          ! Build weights for all frequencies including occupation number differences
          do iw=wi,wf
            ! Get resonant weight
            !   Get energy denominator 
            ! Note: For 'screen' scis12c=0, brd=0 and w only contains w(1) = 0
            !       --> zt1 = e_{ok} - e_{uk+q}
            zt1=w(iw)+deou(ist1, ist2)+scis12c(ist1,ist2)+zi*brd
            !   Check for vanishing denominators in case of screening (no broadening)
            if(abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
! Discuss: Should the weight be wkpt0?
            ! Note: For 'screen' : wou = (f_{ok} - f_{uk+q})/(e_{ok} - e_{uk+q})
            wou(iw,ist1,ist2) = docc12(ist1, ist2) * wkpt(ik) / omega / zt1
            if(doares) then
              ! Get anti-resonant weight 
              ! Note: tordf defaults to 1, yielding the retarded response. 
              !       For the time orderded resoponse set input%xs%tddft%torddf=true (tordf = -1)
              ! Note: For 'screen': zt1 = e_{uk} - e_{ok+q}
              !       or for t.r. sym. : zt1 = e_{ok} - e_{uk+q}
              zt1=w(iw)+deuo(ist2, ist1)+scis21c(ist2,ist1)+tordf*zi*brd
              if(abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
              ! Note: For 'screen' : wuo = (f_{uk} - f_{ok+q})/(e_{uk} - e_{ok+q})
              !       or for t.r. sym: wuo = wou (relies of w=0 and brd=0)
              wuo(iw,ist1,ist2) = docc21(ist2, ist1) * wkpt(ik) / omega / zt1
            end if
          end do

          ! Save weights for current ist1/ist2 combination
          ! to use in the treatment of head and wing components
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
        ! zv_ou = v^{1/2} M_ou 
        zvou(:) = xiou(ist1, ist2, :)
        if(doares) then
          ! zv_uo = v^{1/2} M_uo 
          ! Note: For 'screen' and t.r. sym xiou(1,2,:) = xiuo(2,1,:)
          !       i.e. zvou(:) = zvuo(:)
          zvuo(:) = xiuo(ist2, ist1, :)
        end if

        ! Treatments of head and wings of of chi
        do iw = wi, wf

          G0: if(tq0) then

            do oct1 = 1, 3

              if(doares) then
                ! Wings G=0,G'/=0,q=0
                chi0w(2:, 1, oct1, iw-wi+1) = chi0w(2:, 1, oct1, iw-wi+1)&
                  &+ wouw(iw) * pmou(oct1, ist1, ist2) * conjg(zvou(2:))&
                  &+ wuow(iw) * pmuo(oct1, ist2, ist1) * conjg(zvuo(2:))
                chi0w(2:, 2, oct1, iw-wi+1) = chi0w(2:, 2, oct1, iw-wi+1)&
                  &+ wouw(iw) * zvou(2:) * conjg(pmou(oct1, ist1, ist2))&
                  &+ wuow(iw) * zvuo(2:) * conjg(pmuo(oct1, ist2, ist1))
              else
                ! Wings G=0,G'/=0,q=0
                chi0w(2:, 1, oct1, iw-wi+1) = chi0w(2:, 1, oct1, iw-wi+1)&
                  &+ wouw(iw) * pmou(oct1, ist1, ist2) * conjg(zvou(2:))
                chi0w(2:, 2, oct1, iw-wi+1) = chi0w(2:, 2, oct1, iw-wi+1)&
                  &+ wouw(iw) * zvou(2:) * conjg(pmou(oct1, ist1, ist2))
              end if

              do oct2 = 1, 3

                if(doares) then
                  ! Head G=0,G'=0,q=0
                  if(.not.input%xs%tddft%ahc) then
                    chi0h(oct1, oct2, iw-wi+1) = chi0h(oct1, oct2, iw-wi+1)&
                      &+ wouh(iw) * pmou(oct1, ist1, ist2) * conjg(pmou(oct2, ist1, ist2))&
                      &+ wuoh(iw) * pmuo(oct1, ist2, ist1) * conjg(pmuo(oct2, ist2, ist1))
                  else
                    ! Treatment for anomalous hall conductivity
                    winv=1.0d0/(w(iw)+zi*brd)
                    if(abs(w(iw)).lt.1.d-8) winv=1.d0
                    chi0h(oct1, oct2, iw-wi+1) = chi0h(oct1, oct2, iw-wi+1)&
                      &+ wouh(iw) * pmou(oct1, ist1, ist2) * conjg(pmou(oct2, ist1, ist2))&
                      &* deou(ist1, ist2) * winv&
                      &+ wuoh(iw) * pmuo(oct1, ist2, ist1) * conjg(pmuo(oct2, ist2, ist1))&
                      &* deuo(ist2, ist1) * winv 
                  end if
                else
                  ! Head G=0,G'=0,q=0
                  chi0h(oct1, oct2, iw-wi+1) = chi0h(oct1, oct2, iw-wi+1)&
                    &+ wouh(iw) * pmou(oct1, ist1, ist2) * conjg(pmou(oct2, ist1, ist2))
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

    ! Calculate resonant contributions to chi (actually chi^*)
    ! Allocate helper array for matrix matrix multiplication 
    allocate(zm(n,nst1,nst2))
    do iw=wi,wf
      ! zm_ou(G) = ( w_ou(omega) * \tilde{M}_{ou}(G) )^*
      ! Note: For 'screen' only wou(0) is used and it is real-valued since no broadening is used
      do ist2 = 1, nst2
        do ist1 = 1, nst1
          zm(:,ist1,ist2)=conjg(wou(iw,ist1,ist2)*xiou(ist1,ist2,:))
        end do
      end do
      ! Chi0(omega)_{GG'} = Sum_{ou} zmou_{G ou} * xiou{ou G'} + Chi0(omega)_{GG'}
      !   The lapack routine considers zm to be of the shape zm(n,*)
      !   and xiou to be in the form xiou(*,n). 
      !   That means: zm(:,ist1,ist2) -> zm(:, ist1+(ist2-1)*nst1)
      !               xiou(ist1,ist2,:) -> xiou(ist1+(ist2-1)*nst1,:)
      !   So the sum over states is taken care of.
      call zgemm('n', 'n', n, n, nst1*nst2, zone,&
        & zm(1,1,1), n, xiou(1,1,1), nst1*nst2, zone, chi0(1,1,iw-wi+1), n)
    ! Energy loop
    end do

    if( doares ) then
      ! Calculate anti-resonant contributions to chi (actually chi^*)
      ! Allocate helper array for matrix matrix multiplication anti-resonant
      deallocate(zm)
      allocate(zm(n,nst2,nst1))
      do iw=wi,wf
        ! Note: For 'screen' only wuo(0) is used and it is real-valued since no broadening is used
        !       if t.r. sym. is used wuo = wou and xiuo = xiou. So the matrix multiplication 
        !       just adds the same as the resonant part.
        do ist2 = 1, nst2
          do ist1 = 1, nst1
            zm(:,ist2,ist1)=conjg(wuo(iw,ist1,ist2)*xiuo(ist2,ist1,:))
          end do
        end do
        call zgemm('n', 'n', n, n, nst1*nst2, zone, zm(1,1,1),&
          & n, xiuo(1,1,1), nst1*nst2, zone, chi0(1,1,iw-wi+1), n)
      end do
    end if

    ! Timing
    call timesec(cpu1)
    cpuupd = cpuupd + cpu1 - cpu0

    ! Deallocate weights and helper arrays
    deallocate(wou,wuo,zm)

!*****************************************************************************************************

    ! Timing information
    cputot = cpuread + cpuosc + cpuupd
    call dftim(iq, ik, trim(fnxtim), cpuread, cpuosc, cpuupd, cputot)

    ! Why is taht nescessary?
    ! Synchronize
    if( .not. tscreen) call barrier(callername=trim(thisname))

  ! End loop over k-points
  end do kloop

  ! The way the helper array zm was constructed, actually
  ! the complex conjugated of chi was computed. Now fix that:
  ! Also multiply by 2 to include antires contributions, if
  ! they were omitted
  if(.not. doares) then 
    chi0=2.0d0*conjg(chi0)
    chi0w = 2.0d0 * chi0w
    chi0h = 2.0d0 * chi0h
  else
    chi0=conjg(chi0)
  end if

  ! Semi-classical Drude approximation to intraband terms
  ! Default drude = [0.0 0.0]
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

  ! Screen - Deallocate eigenvalue and eigenvector related arrays
  if(tscreen) call ematqdealloc

  ! Symmetrize head, with respect to the crystal symmetry.
  head: if(tq0) then
    allocate(chi0hs(3, 3, nwdfp), eps0(3, 3, nwdf))

    ! Write coulomb-symmetrized dielectric tensor to file (lattice-unsymmetrized)
    !   \epsilon_{i,j}(\omega) = 1 - \tilde{\chi}_{i,j}(\omega)
    forall(iw=1:nwdf)
      eps0(:, :, iw) = dble(krondelta) - chi0h(:, :, iw)
    end forall

    if(rank .eq. 0)&
      & call writedielt('DIELTENS0_NOSYM', 1, [0.d0], eps0(:, :, 1), 0)

    ! Symmetrize the head of \chi, w.r.t. the lattice
    do oct1 = 1, 3
      do oct2 = 1, 3
        call symt2app(oct1, oct2, nwdfp, symt2, chi0h, chi0hs(oct1, oct2, :))
      end do
    end do
    ! Re-assign the symmetrized head
    chi0h(:, :, :) = chi0hs(:, :, :)

    ! Write symmetrized dielectric tensor to file
    forall(iw=1:nwdf)
      eps0(:, :, iw) = dble(krondelta) - chi0hs(:, :, iw)
    end forall
    if(rank .eq. 0)&
     & call writedielt('DIELTENS0', 1, [0.d0], eps0(:, :, 1), 0)

    deallocate(chi0hs, eps0)

  end if head 

  ! Write response function to file
  if(tscreen) then
    ! Write out symmetrized dielectric function/tensor 

    ! Write out static (\omega = 0) screening for current q to text file
    ! eps = 1 - chi0
    call getunit(un)
    open(un, file=trim(fnscreen), form='formatted', action='write', status='replace')
    call putscreen(un, tq0, n, chi0(:, :, 1), chi0h(:, :, 1), chi0w(:, :, :, 1))
    close(un)

    ! Parallel output of q- and w-dependent \epsilon=1-\chi to direct access file
    ! Make microscopic epsilon matrix
    ! Head 
    forall(iw=1:nwdf)
      chi0h(:, :, iw) = dble(krondelta) - chi0h(:, :, iw)
    end forall
    ! Wings
    chi0w = -chi0w
    ! Body
    chi0 = -chi0
    forall(igq=1:n)
      chi0(igq,igq,:) = zone + chi0(igq,igq,:)
    end forall
    ! Write to direct access file
    do iw = wi, wf
      ! It uses the reduced set, but the grid parameters are saved in 
      ! mod_qpoint, where in the case of task >430 the non-reduces parameters 
      ! are saved...
      call puteps0(reduced=.false.,&
        & iq=iq, iw=iw, w=wreal(iw-wi+1),&
        & eps0=chi0(:,:,iw-wi+1), eps0wg=chi0w(:,:,:,iw-wi+1),&
        & eps0hd=chi0h(:,:,iw-wi+1), fname=fneps0)
    end do
  ! Not tscreen
  else
    ! Parallel output of frequency dependent \chi to direct access file
    do j = 0, procs - 1
      if(rank .eq. j) then
        do iw = wi, wf
          call putx0(tq0, iq, iw-wi+1, trim(fnchi0_t), '',&
            & ch0=chi0(:, :, iw-wi+1), ch0wg=chi0w(:, :, :, iw-wi+1),&
            & ch0hd=chi0h(:, :, iw-wi+1))
        end do
      end if
      ! Why is that nescessary
      call barrier(callername=trim(thisname))
    end do
  end if

  ! Deallocations
  deallocate(chi0, chi0h, chi0w)
  deallocate(docc12, docc21, scis12, scis21, scis12c, scis21c)
  deallocate(deou, deuo, wouw, wuow, wouh, wuoh, zvou, zvuo)
  deallocate(xiou, xiuo, pmou, pmuo)
  deallocate(bsedg)
  deallocate(w, wreal)

  ! TETRA
  if(input%xs%tetra%tetradf) then
#ifdef TETRA      
    deallocate(cw, cwa, cwsurf)
    if(input%xs%tetra%cw1k) deallocate(cwt, cw1k, cwa1k, cwsurf1k)
#endif        
  end if

end subroutine dfq
!EOC
