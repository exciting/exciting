
! Copyright (C) 2005-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dfq
! !INTERFACE:
Subroutine dfq (iq)
! !USES:
use ioarray
      Use modinput
      Use mod_constants
      Use mod_eigenvalue_occupancy
      Use mod_misc
      Use mod_Gkvector
      Use mod_APW_LO
      Use mod_Gvector
      Use mod_kpoint
      Use mod_qpoint
      Use mod_lattice
      Use mod_DOS_optics_response
      Use modxs
      Use modtetra
      Use modmpi
      Use m_genwgrid
      Use m_getpemat
      Use m_dftim
      Use m_gettetcw
      Use m_putx0
      Use m_getunit
      Use m_writevars
      Use m_filedel
      Use m_genfilname
! !DESCRIPTION:
!   Calculates the symmetrized Kohn-Sham response function $\chi^0_{\bf{GG'}}
!   ({\bf q},\omega)$ for one ${\bf q}$-point according to
!   $$  \chi^0_{\bf{GG'}}({\bf q},\omega) = \sum_{ou{\bf k}} \left[
!      M^{\bf G}_{ou{\bf k}}({\bf q}) M^{\bf G'}_{ou{\bf k}}({\bf q})^*
!      w_{ou{\bf k}}({\bf q},\omega) +
!      M^{\bf G}_{uo{\bf k}}({\bf q}) M^{\bf G'}_{uo{\bf k}}({\bf q})^*
!      w_{uo{\bf k}}({\bf q},\omega) \right]
!   $$
!   It is related to the Kohn-Sham response function
!    $\chi^0_{\bf{GG'}}({\bf q},\omega)$ by
!   $$ \chi^0_{\bf{GG'}}({\bf q},\omega) = v^{-\frac{1}{2}}_{\bf G}({\bf q})
!      \bar{\chi}^0_{\bf{GG'}}({\bf q},\omega)
!      v^{-\frac{1}{2}}_{\bf G'}({\bf q}) $$
!   and is well defined in the limit as ${\bf q}$ goes to zero.
!   The symmetrized matrix elements are defined as
!   $$   M^{\bf G}_{ou{\bf k}}({\bf q}) =
!         v^{-\frac{1}{2}}_{\bf G}({\bf q})
!         \bar{M}^{\bf G}_{ou{\bf k}}({\bf q}), $$
!   where
!   $$   \bar{M}^{\bf G}_{ou{\bf k}}({\bf q}) =
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
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq
  ! local variables
      Character (*), Parameter :: thisnam = 'dfq'
      Character (256) :: fnscreen
      Real (8), Parameter :: epstetra = 1.d-8
      Complex (8), Allocatable :: w (:)
      Complex (8), Allocatable :: chi0 (:, :, :), hdg (:, :, :)
      Complex (8), Allocatable :: chi0w (:, :, :, :), chi0h (:, :, :), &
     & eps0 (:, :, :)
      Complex (8), Allocatable :: chi0hAHC(:, :)
      Complex (8), Allocatable :: wou (:,:,:), wuo (:,:,:), wouw (:), wuow (:), &
     & wouh (:), wuoh (:)
      Complex (8), Allocatable :: zvou (:), zvuo (:), chi0hs (:, :, :), &
     & bsedg (:, :), scis12c (:, :), scis21c (:, :),zm(:,:,:)
      Real (8), Allocatable :: wreal (:), cw (:), cwa (:), cwsurf (:)
      Real (8), Allocatable :: cwt (:, :), cw1k (:, :, :), cwa1k (:, :, &
     & :), cwsurf1k (:, :, :)
      Real (8), Allocatable :: scis12 (:, :), scis21 (:, :)
      Complex(8) :: zt1, winv
      Real (8) :: brd, cpu0, cpu1, cpuread, cpuosc, cpuupd, cputot, wintv(2), wplas, wrel
      Integer :: n, j, i1, i2, j1, j2, ik, ikq, igq, iw, wi, wf, ist1, &
     & ist2, nwdfp
      Integer :: oct1, oct2, un
      Logical :: tq0
      Logical, External :: tqgamma, transik, transijst
      real(8) :: ta,tb,tc,td

      call timesec(ta)     

      If (input%xs%tddft%acont .And. tscreen) Then
         Write (*,*)
         Write (*, '("Error(", a, "): analytic continuation does not wo&
        &rk for screening")')
         Write (*,*)
         Call terminate
      End If
      tfxcbse = ((input%xs%tddft%fxctypenumber .Eq. 7) .Or. &
     & (input%xs%tddft%fxctypenumber .Eq. 8)) .And. ( .Not. tscreen)
  ! sampling of Brillouin zone
      bzsampl = 0
      If (input%xs%tetra%tetradf) bzsampl = 1
  ! initial and final w-point
      wi = wpari
      wf = wparf
      nwdfp = wf - wi + 1
  ! matrix size for response function
      n = ngq (iq)
  ! zero broadening for analytic contiunation
      brd = input%xs%broad
      If (input%xs%tddft%acont) brd = 0.d0
  ! zero broadening for dielectric matrix (w=0) for band-gap systems
      If (task .Eq. 430) brd = 0.d0
  ! file extension for q-point
      If ( .Not. tscreen) Call genfilname (iqmt=iq, setfilext=.True.)
  ! filenames for output
      If (tscreen) Then
         Call genfilname (basename='TETW', iq=iq, appfilext=.True., &
        & filnam=fnwtet)
         Call genfilname (basename='PMAT', appfilext=.True., &
        & filnam=fnpmat)
         Call genfilname (basename='SCREEN', bzsampl=bzsampl, iq=iq, &
        & filnam=fnscreen)
         Call genfilname (nodotpar=.True., basename='EMAT_TIMING', &
        & iq=iq, etype=input%xs%emattype, procs=procs, rank=rank, &
        & appfilext=.True., filnam=fnetim)
         Call genfilname (nodotpar=.True., basename='X0_TIMING', iq=iq, &
        & procs=procs, rank=rank, appfilext=.True., filnam=fnxtim)
      Else
         Call genfilname (basename='TETW', iqmt=iq, filnam=fnwtet)
         Call genfilname (basename='PMAT_XS', filnam=fnpmat)
         Call genfilname (basename='EMAT', iqmt=iq, filnam=fnemat)
         Call genfilname (nodotpar=.True., basename='X0_TIMING', &
        & bzsampl=bzsampl, iqmt=iq, procs=procs, rank=rank, &
        & filnam=fnxtim)
         Call genfilname (basename='X0', bzsampl=bzsampl, &
        & acont=input%xs%tddft%acont, nar= .Not. input%xs%tddft%aresdf, &
        & tord=input%xs%tddft%torddf, markfxcbse=tfxcbse, iqmt=iq, &
        & filnam=fnchi0)
         Call genfilname (basename='X0', bzsampl=bzsampl, &
        & acont=input%xs%tddft%acont, nar= .Not. input%xs%tddft%aresdf, &
        & tord=input%xs%tddft%torddf, markfxcbse=tfxcbse, iqmt=iq, &
        & procs=procs, rank=rank, filnam=fnchi0_t)
      End If
  ! remove timing files from previous runs
      Call filedel (trim(fnxtim))
  ! calculate k+q and G+k+q related variables
      Call init1offs (qvkloff(1, iq))
  ! generate link array for tetrahedra
      If (input%xs%tetra%tetradf) Call gentetlinkp (vql(1, iq), &
     & input%xs%tetra%qweights)
  ! find highest (partially) occupied and lowest (partially) unoccupied states
      Call findocclims (iq, istocc0, istocc, istunocc0, istunocc, &
     & isto0, isto, istu0, istu)
  ! find limits for band combinations
      Call ematbdcmbs (input%xs%emattype)
  ! check if q-point is Gamma point
      tq0 = tqgamma (iq)
      If (tq0) Then
         Write (unitout, '(a)') 'Info(' // trim (thisnam) // '): Gamma &
        &q - point: using momentum matrix elements for dielectric funct&
        &ion'
      End If
  ! write out matrix size of response function
      Write (unitout, '(a, i6)') 'Info(' // thisnam // '): number of G &
     &+ q vectors (local field effects):', ngq (iq)
      Write (unitout, '(a, 4i6)') 'Info(' // thisnam // '): lowest (par&
     &tially)  unoccupied state: ', istunocc0
      Write (unitout, '(a, 4i6)') 'Info(' // thisnam // '): highest (pa&
     &rtially) occupied state  : ', istocc0
      Write (unitout, '(a, 4i5)') 'Info(' // thisnam // '): &
     &band-combination limits  nst1,  nst2,  nst3,  nst4:',  nst1,  nst2,  nst3,  nst4 
      Write (unitout, '(a, 4i5)') 'Info(' // thisnam // '): &
     &band-combination limits istl1, istu1, istl2, istu2:', istl1, istu1, istl2, istu2
      Write (unitout, '(a, 4i5)') 'Info(' // thisnam // '): &
     &band-combination limits istl3, istu3, istl4, istu4:', istl3, istu3, istl4, istu4
  ! allocate arrays for eigenvalue and occupation number differences
      If (allocated(deou)) deallocate (deou)
      Allocate (deou(nst1, nst2))
      If (allocated(deuo)) deallocate (deuo)
      Allocate (deuo(nst3, nst4))
      If (allocated(docc12)) deallocate (docc12)
      Allocate (docc12(nst1, nst2))
      If (allocated(docc21)) deallocate (docc21)
      Allocate (docc21(nst3, nst4))
  ! allocate matrix elements arrays
      If (allocated(xiou)) deallocate (xiou)
      Allocate (xiou(nst1, nst2, n))
      If (allocated(xiuo)) deallocate (xiuo)
      Allocate (xiuo(nst3, nst4, n))
      If (allocated(pmou)) deallocate (pmou)
      Allocate (pmou(3, nst1, nst2))
      If (allocated(pmuo)) deallocate (pmuo)
      Allocate (pmuo(3, nst3, nst4))
  ! allocate arrays
      Allocate (hdg(nst1, nst2, nkpt))
      Allocate (scis12(nst1, nst2),scis12c(nst1, nst2))
      Allocate (scis21(nst2, nst1),scis21c(nst2, nst1))
      Allocate (w(nwdf))
      Allocate (wreal(nwdfp))
      Allocate (chi0h(3, 3, nwdfp))
      Allocate (chi0hAHC(3, 3))
      Allocate (chi0w(n, 2, 3, nwdfp))
      Allocate (chi0(n, n, nwdfp))
!      Allocate (wou(nwdf))
!      Allocate (wuo(nwdf))
      Allocate (wouw(nwdf), wuow(nwdf), wouh(nwdf), wuoh(nwdf))
      Allocate (zvou(n), zvuo(n))
      Allocate (bsedg(nst1, nst2))
      scis12 (:, :) = 0.d0
      scis21 (:, :) = 0.d0
      bsedg(:,:)=zzero
      If (input%xs%tetra%tetradf) Then
         Allocate (cw(nwdf), cwa(nwdf), cwsurf(nwdf))
         If (input%xs%tetra%cw1k) allocate (cwt(nstsv, nstsv), &
        & cw1k(nst1, nst2, nwdfp), cwa1k(nst1, nst2, nwdfp), &
        & cwsurf1k(nst1, nst2, nwdfp))
      End If
  ! generate complex energy grid
      wintv(1)=input%xs%energywindow%intv(1)
      wintv(2)=input%xs%energywindow%intv(2)
  ! for calculation of static screening the first frequency point should be zero
      if (task.eq.430) wintv(1)=0.d0
      Call genwgrid (nwdf, wintv, &
     & input%xs%tddft%acont, 0.d0, w_cmplx=w)
      wreal (:) = dble (w(wi:wf))
      If (wreal(1) .Lt. epstetra) wreal (1) = epstetra
  ! initializations
      chi0 (:, :, :) = zzero
      chi0w (:, :, :, :) = zzero
      chi0h (:, :, :) = zzero
      chi0hAHC (:, :) = zzero
      If (tscreen) Then
     ! generate radial integrals wrt. sph. Bessel functions
        call timesec(tc)
         Call ematrad (iq)
        call timesec(td)
!        write(*,*) td-tc

     ! delete timing information of previous runs
         Call filedel (trim(fnetim))
     ! write information
         Write (unitout, '(a, i6)') 'Info(' // thisnam // '): number of&
        & G + q vectors:', ngq (iq)
         Call ematqalloc
      End If
      If (tfxcbse) Then
         Call getbsediag
         Write (unitout, '("Info(", a, "): read diagonal of BSE kernel"&
        &)') trim (thisnam)
         Write (unitout, '(" mean value : ", 2g18.10)') bsed
         bsedg (:, :) = bsed
      End If
      call timesec(tb)
!      write(*,*) tb-ta
!      stop
!write(*,*) wi,wf
!write(*,*) 'starting loop'
  ! loop over k-points
      Do ik = 1, nkpt
     ! k-point analysis

         If ( .Not. transik(ik)) Cycle
         Call chkpt (3, (/ task, iq, ik /), 'dfq: task, q-point index, &
        &k-point index')
         cpuosc = 0.d0
         cpuupd = 0.d0
         Call timesec (cpu0)
         ikq = ikmapikq (ik, iq)

         Call getdevaldoccsv (iq, ik, ikq, istl1, istu1, istl2, istu2, &
        & deou, docc12, scis12)
         Call getdevaldoccsv (iq, ik, ikq, istl2, istu2, istl1, istu1, &
        & deuo, docc21, scis21)
        scis12c(:,:)=scis12(:,:)
        scis21c(:,:)=scis21(:,:)

         If (tscreen) Then
        ! do not use scissors correction for screening
            If (task .Eq. 430) Then
               scis12c (:, :) = zzero
               scis21c (:, :) = zzero
            End If
        ! for screening calculate matrix elements of plane wave on the fly
            Call ematqk1 (iq, ik)
            If ( .Not. allocated(xiuo)) allocate (xiuo(nst3, nst4, n))
            If ( .Not. allocated(pmuo)) allocate (pmuo(3, nst3, nst4))
         End If
     ! add BSE diagonal shift use with BSE-kernel
         scis12c (:, :) = scis12c (:, :) + bsedg (:, :)
         scis21c (:, :) = scis21c (:, :) + transpose (bsedg(:, :))
     ! get matrix elements (exp. expr. or momentum op.)
!write(*,*) 'starting getpemat'
         Call getpemat (iq, ik, trim(fnpmat), trim(fnemat), m12=xiou, &
        & m34=xiuo, p12=pmou, p34=pmuo)
!write(*,*) 'done'
      
     ! set matrix elements to one for Lindhard function
         If (input%xs%tddft%lindhard) Then
       ! set G=0 components to one
            xiou (:, :, 1) = zone
            xiuo (:, :, 1) = zone
       ! set G/=0 components to zero
            If (n .Gt. 1) Then
               xiou (:, :, 2:) = zzero
               xiuo (:, :, 2:) = zzero
            End If
       ! set momentum matrix elements to one
            pmou (:, :, :) = zone
            pmuo (:, :, :) = zone
         End If
         If (input%xs%tetra%cw1k) Then
            Do iw = 1, nwdfp
               Call tetcwifc_1k (ik, nkpt, nstsv, evalsv, efermi, &
              & wreal(iw), 2, cwt)
               cw1k (:, :, iw) = cwt (istl1:istu1, istl2:istu2)
               Call tetcwifc_1k (ik, nkpt, nstsv, evalsv, &
              & efermi,-wreal(iw), 2, cwt)
               cwa1k (:, :, iw) = cwt (istl1:istu1, istl2:istu2)
               Call tetcwifc_1k (ik, nkpt, nstsv, evalsv, efermi, &
              & wreal(iw), 4, cwt)
               cwsurf1k (:, :, iw) = cwt (istl1:istu1, istl2:istu2)
            End Do
         End If
         If (tscreen) Then
        ! we don't need anti-resonant parts here, assign them the same
        ! value as for resonant parts, resulting in a factor of two.
            Do igq = 1, n
               xiuo (:, :, igq) = transpose (xiou(:, :, igq))
            End Do
            Do j = 1, 3
               pmuo (j, :, :) = transpose (pmou(j, :, :))
            End Do
            deuo (:, :) = transpose (deou(:, :))
            docc21 (:, :) = transpose (docc12(:, :))
            scis21c (:, :) = transpose (scis12c(:, :))
         End If
     ! turn off antiresonant terms (type 2-1 band combiantions) for Kohn-Sham
     ! response function
         If (( .Not. input%xs%tddft%aresdf) .And. ( .Not. tscreen)) &
        & Then
            xiuo (:, :, :) = zzero
            pmuo (:, :, :) = zzero
         End If
         Do ist1 = 1, istocc0 - istunocc0 + 1
            Do ist2 = 1, istocc0 - istunocc0 + 1
               j = ist1 + istunocc0 - 1
           ! set lower triangle of first block to zero
               If (ist1 .Gt. ist2) Then
                  xiou (j, ist2, :) = zzero
                  pmou (:, j, ist2) = zzero
               End If
           ! set diagonal to zero (project out intraband contributions)
               If (( .Not. input%xs%tddft%intraband) .And. (ist1 .Eq. &
              & ist2)) Then
                  xiou (j, ist2, :) = zzero
                  pmou (:, j, ist2) = zzero
               End If
           ! set upper triangle of second block to zero
           ! also set diagonal to zero to avoid double counting
               If (ist1 .Ge. ist2) Then
                  xiuo (ist2, j, :) = zzero
                  pmuo (:, ist2, j) = zzero
               End If
            End Do
         End Do
         Call timesec (cpu1)
         cpuread = cpu1 - cpu0


!**********************************************************************************

if (.false.) then
      Allocate (wou(nwdf,1,1))
      Allocate (wuo(nwdf,1,1))
!write(*,*) nst1,nst2
         Do ist1 = 1, nst1
!write(*,*) 'ist2 loop'
!call timesec(ta)
            Do ist2 = 1, nst2

           !---------------------!
           !     denominator     !
           !---------------------!
           ! absolute band indices
               i1 = ist1
               i2 = istunocc0 + ist2 - 1
           ! band analysis
               If ( .Not. transijst(ik, i1, i2)) Cycle
               Call timesec (cpu0)
           ! user request termination
               Call terminateqry ('dfq')
               If (input%xs%tetra%tetradf) Then
              ! mirror index pair on diagonal if necessary
                  If (i1 .Gt. i2) Then
                     j1 = ist2
                     j2 = ist1 - istunocc0 + 1
                  Else
                     j1 = ist1
                     j2 = ist2
                  End If
              ! read weights for tetrahedron method
                  If (input%xs%tetra%cw1k) Then
                     cw (wi:wf) = cw1k (ist1, ist2, :)
                     cwa (wi:wf) = cwa1k (ist1, ist2, :)
                     cwsurf (wi:wf) = cwsurf1k (ist1, ist2, :)
                  Else
                     Call gettetcw (iq, ik, j1, j2, nst1, nst2, nwdf, &
                    & trim(fnwtet), cw, cwa, cwsurf)
                  End If
              ! include occupation number differences
                  wou (wi:wf,1,1) = docc12 (ist1, ist2) * cmplx (cw(wi:wf), &
                 & cwsurf(wi:wf), 8) / omega
                  wuo (wi:wf,1,1) = - docc21 (ist2, ist1) * cmplx &
                 & (cwa(wi:wf), 0.d0, 8) / omega
                  If (tq0) Then
                 ! rescale: use delta-function delta(e_nmk + scis_nmk - w)
                 ! take real part of BSE diagonal (being contained in scis12c)
                 ! since tetrahedron method formalism implemented does not allow
                 ! otherwise
                     wouw (wi:wf) = cmplx (dble(wou(wi:wf,1,1)), &
                    & aimag(wou(wi:wf,1,1))*deou(ist1, &
                    & ist2)/(-wreal(:)-dble(scis12c(ist1, ist2))))
                     wuow (wi:wf) = cmplx (dble(wuo(wi:wf,1,1)), &
                    & aimag(wuo(wi:wf,1,1))*deuo(ist2, &
                    & ist1)/(-wreal(:)-dble(scis21c(ist2, ist1))))
                     wouh (wi:wf) = cmplx (dble(wou(wi:wf,1,1)), &
                    & aimag(wou(wi:wf,1,1))*deou(ist1, &
                    & ist2)**2/(-wreal(:)-dble(scis12c(ist1, ist2)))**2)
                     wuoh (wi:wf) = cmplx (dble(wuo(wi:wf,1,1)), &
                    & aimag(wuo(wi:wf,1,1))*deuo(ist2, &
                    & ist1)**2/(-wreal(:)-dble(scis21c(ist2, ist1)))**2)
                  End If
               Else
              ! include occupation number differences
                  do iw=wi,wf
                     ! check for vanishing denominators in case of screening
                     ! (no broadening)
                     zt1=w(iw)+deou(ist1, ist2)+scis12c(ist1,ist2)+zi*brd
                     if (abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
                     wou (iw,1,1) = docc12 (ist1, ist2) * wkpt (ik) / omega / zt1
                     zt1=w(iw)+deuo(ist2, ist1)+scis21c(ist2,ist1)+tordf*zi*brd
                     if (abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
                     wuo (iw,1,1) = docc21 (ist2, ist1) * wkpt (ik) / omega / zt1
                  end do
                  wouw (wi:wf) = wou (wi:wf,1,1)
                  wuow (wi:wf) = wuo (wi:wf,1,1)
                  wouh (wi:wf) = wou (wi:wf,1,1)
                  wuoh (wi:wf) = wuo (wi:wf,1,1)
               End If
               call timesec(cpu1)
               cpuosc=cpuosc+cpu1-cpu0
           !----------------------------------!
           !     update response function     !
           !----------------------------------!
!write(*,*) 'zgerc'
!call timesec(ta)

               zvou (:) = xiou (ist1, ist2, :)
               zvuo (:) = xiuo (ist2, ist1, :)

!write(*,*) 'zgerc'
!call timesec(ta)

               Do iw = wi, wf
!write(*,*) 'zgerc'
!call timesec(ta)
                  ! body
                 Call zgerc (n, n, wou(iw,1,1), zvou, 1, zvou, 1, chi0(:, &
                       & :, iw-wi+1), n)
                 Call zgerc (n, n, wuo(iw,1,1), zvuo, 1, zvuo, 1, chi0(:, &
                      & :, iw-wi+1), n)
!call timesec(tb)
!write(*,*) tb-ta
!write(*,*) 'zgerc done'
!read(*,*)
                  If (tq0) Then
                     Do oct1 = 1, 3
                        ! wings
                        chi0w (2:, 1, oct1, iw-wi+1) = chi0w (2:, 1, &
                             & oct1, iw-wi+1) + wouw (iw) * pmou (oct1, ist1, &
                             & ist2) * conjg (zvou(2:)) + wuow (iw) * pmuo &
                             & (oct1, ist2, ist1) * conjg (zvuo(2:))
                        chi0w (2:, 2, oct1, iw-wi+1) = chi0w (2:, 2, &
                             & oct1, iw-wi+1) + wouw (iw) * zvou (2:) * conjg &
                             & (pmou(oct1, ist1, ist2)) + wuow (iw) * zvuo &
                             & (2:) * conjg (pmuo(oct1, ist2, ist1))
                        Do oct2 = 1, 3
                           ! head
                           If(.Not.input%xs%tddft%ahc) Then
                              chi0h (oct1, oct2, iw-wi+1) = chi0h (oct1, &
                                   & oct2, iw-wi+1) + wouh (iw) * pmou (oct1, &
                                   & ist1, ist2) * conjg (pmou(oct2, ist1, &
                                   & ist2)) + wuoh (iw) * pmuo (oct1, ist2, &
                                   & ist1) * conjg (pmuo(oct2, ist2, ist1))
                           Else
                              winv=1.0d0/(w(iw)+zi*brd)
                              If (Abs(w(iw)).Lt.1.d-8) winv=1.d0
                              chi0h (oct1, oct2, iw-wi+1) = chi0h (oct1, oct2, iw-wi+1) + &
                                   & wouh (iw) * pmou (oct1, ist1, ist2) * conjg (pmou(oct2, ist1, ist2))*&
                                   & (deou(ist1, ist2)*winv) + &
                                   & wuoh (iw) * pmuo (oct1, ist2, ist1) * conjg (pmuo(oct2, ist2, ist1))*&
                                   & (deuo(ist2, ist1)*winv) 
                              
                           End If

                        End Do
                     End Do
                  End If
               End Do
!call timesec(tb)
!write(*,*) tb-ta
!write(*,*) 'zgerc done'
!read(*,*)
               Call timesec (cpu0)
               cpuupd = cpuupd + cpu0 - cpu1
           ! end loop over states combinations
            End Do
!call timesec(tb)
!write(*,*) tb-ta
!write(*,*) 'loop done'
!read(*,*)
         End Do
!write(*,*) sum(chi0(:,:,1))
deallocate(wuo,wou)
!*****************************************************************************************************
else
      Allocate (wou(nwdf,nst1,nst2))
      Allocate (wuo(nwdf,nst1,nst2))
         Do ist1 = 1, nst1
            Do ist2 = 1, nst2

           !---------------------!
           !     denominator     !
           !---------------------!
           ! absolute band indices
               i1 = ist1
               i2 = istunocc0 + ist2 - 1
           ! band analysis
               If ( .Not. transijst(ik, i1, i2)) Cycle
               Call timesec (cpu0)
           ! user request termination
               Call terminateqry ('dfq')
               If (input%xs%tetra%tetradf) Then
              ! mirror index pair on diagonal if necessary
                  If (i1 .Gt. i2) Then
                     j1 = ist2
                     j2 = ist1 - istunocc0 + 1
                  Else
                     j1 = ist1
                     j2 = ist2
                  End If
              ! read weights for tetrahedron method
                  If (input%xs%tetra%cw1k) Then
                     cw (wi:wf) = cw1k (ist1, ist2, :)
                     cwa (wi:wf) = cwa1k (ist1, ist2, :)
                     cwsurf (wi:wf) = cwsurf1k (ist1, ist2, :)
                  Else
                     Call gettetcw (iq, ik, j1, j2, nst1, nst2, nwdf, &
                    & trim(fnwtet), cw, cwa, cwsurf)
                  End If
              ! include occupation number differences
                  wou (wi:wf,ist1,ist2) = docc12 (ist1, ist2) * cmplx (cw(wi:wf), &
                 & cwsurf(wi:wf), 8) / omega
                  wuo (wi:wf,ist1,ist2) = - docc21 (ist2, ist1) * cmplx &
                 & (cwa(wi:wf), 0.d0, 8) / omega
                  If (tq0) Then
                 ! rescale: use delta-function delta(e_nmk + scis_nmk - w)
                 ! take real part of BSE diagonal (being contained in scis12c)
                 ! since tetrahedron method formalism implemented does not allow
                 ! otherwise
                     wouw (wi:wf) = cmplx (dble(wou(wi:wf,ist1,ist2)), &
                    & aimag(wou(wi:wf,ist1,ist2))*deou(ist1, &
                    & ist2)/(-wreal(:)-dble(scis12c(ist1, ist2))))
                     wuow (wi:wf) = cmplx (dble(wuo(wi:wf,ist1,ist2)), &
                    & aimag(wuo(wi:wf,ist1,ist2))*deuo(ist2, &
                    & ist1)/(-wreal(:)-dble(scis21c(ist2, ist1))))
                     wouh (wi:wf) = cmplx (dble(wou(wi:wf,ist1,ist2)), &
                    & aimag(wou(wi:wf,ist1,ist2))*deou(ist1, &
                    & ist2)**2/(-wreal(:)-dble(scis12c(ist1, ist2)))**2)
                     wuoh (wi:wf) = cmplx (dble(wuo(wi:wf,ist1,ist2)), &
                    & aimag(wuo(wi:wf,ist1,ist2))*deuo(ist2, &
                    & ist1)**2/(-wreal(:)-dble(scis21c(ist2, ist1)))**2)
                  End If
               Else
              ! include occupation number differences
                  do iw=wi,wf
                     ! check for vanishing denominators in case of screening
                     ! (no broadening)
                     zt1=w(iw)+deou(ist1, ist2)+scis12c(ist1,ist2)+zi*brd
                     if (abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
                     wou (iw,ist1,ist2) = docc12 (ist1, ist2) * wkpt (ik) / omega / zt1
                     zt1=w(iw)+deuo(ist2, ist1)+scis21c(ist2,ist1)+tordf*zi*brd
                     if (abs(zt1).lt. input%xs%epsdfde) zt1=1.d0
                     wuo (iw,ist1,ist2) = docc21 (ist2, ist1) * wkpt (ik) / omega / zt1
                  end do
                  wouw (wi:wf) = wou (wi:wf,ist1,ist2)
                  wuow (wi:wf) = wuo (wi:wf,ist1,ist2)
                  wouh (wi:wf) = wou (wi:wf,ist1,ist2)
                  wuoh (wi:wf) = wuo (wi:wf,ist1,ist2)
               End If
               call timesec(cpu1)
               cpuosc=cpuosc+cpu1-cpu0
           !----------------------------------!
           !     update response function     !
           !----------------------------------!

               zvou (:) = xiou (ist1, ist2, :)
               zvuo (:) = xiuo (ist2, ist1, :)

               Do iw = wi, wf
                  ! body
!                  Call zgerc (n, n, wou(iw), zvou, 1, zvou, 1, chi0(:, &
!                       & :, iw-wi+1), n)
!                  Call zgerc (n, n, wuo(iw), zvuo, 1, zvuo, 1, chi0(:, &
!                       & :, iw-wi+1), n)
                  If (tq0) Then
                     Do oct1 = 1, 3
                        ! wings
                        chi0w (2:, 1, oct1, iw-wi+1) = chi0w (2:, 1, &
                             & oct1, iw-wi+1) + wouw (iw) * pmou (oct1, ist1, &
                             & ist2) * conjg (zvou(2:)) + wuow (iw) * pmuo &
                             & (oct1, ist2, ist1) * conjg (zvuo(2:))
                        chi0w (2:, 2, oct1, iw-wi+1) = chi0w (2:, 2, &
                             & oct1, iw-wi+1) + wouw (iw) * zvou (2:) * conjg &
                             & (pmou(oct1, ist1, ist2)) + wuow (iw) * zvuo &
                             & (2:) * conjg (pmuo(oct1, ist2, ist1))
                        Do oct2 = 1, 3
                           ! head
                           If(.Not.input%xs%tddft%ahc) Then
                              chi0h (oct1, oct2, iw-wi+1) = chi0h (oct1, &
                                   & oct2, iw-wi+1) + wouh (iw) * pmou (oct1, &
                                   & ist1, ist2) * conjg (pmou(oct2, ist1, &
                                   & ist2)) + wuoh (iw) * pmuo (oct1, ist2, &
                                   & ist1) * conjg (pmuo(oct2, ist2, ist1))
                           Else
                              winv=1.0d0/(w(iw)+zi*brd)
                              If (Abs(w(iw)).Lt.1.d-8) winv=1.d0
                              chi0h (oct1, oct2, iw-wi+1) = chi0h (oct1, oct2, iw-wi+1) + &
                                   & wouh (iw) * pmou (oct1, ist1, ist2) * conjg (pmou(oct2, ist1, ist2))*&
                                   & (deou(ist1, ist2)*winv) + &
                                   & wuoh (iw) * pmuo (oct1, ist2, ist1) * conjg (pmuo(oct2, ist2, ist1))*&
                                   & (deuo(ist2, ist1)*winv) 
                              
                           End If

                        End Do
                     End Do
                  End If
               End Do
               Call timesec (cpu0)
               cpuupd = cpuupd + cpu0 - cpu1
           ! end loop over states combinations
            End Do
         End Do


         allocate(zm(n,nst1,nst2))
         do iw=wi,wf
           Do ist2 = 1, nst2
             Do ist1 = 1, nst1
              zm(:,ist1,ist2)=conjg(wou(iw,ist1,ist2)*xiou(ist1,ist2,:))
!              zm(:,ist1,ist2)=wou(iw,ist1,ist2)*xiou(ist1,ist2,:)
             enddo
           enddo
           
           call zgemm('n', &           ! TRANSA = 'C'  op( A ) = A**H.
                      'n', &           ! TRANSB = 'N'  op( B ) = B.
                       n, &          ! M ... rows of op( A ) = rows of C
                       n, &           ! N ... cols of op( B ) = cols of C
                       nst1*nst2, &          ! K ... cols of op( A ) = rows of op( B )
                       zone, &          ! alpha
!                       xiou(1,1,1), &
!                       nst1*nst2,&
!                       zm(1,1,1), &
!                       n, &
                       zm(1,1,1), &           ! B
                       n, &          ! LDB ... leading dimension of B
                       xiou(1,1,1), &           ! A
                       nst1*nst2,&           ! LDA ... leading dimension of A
                       zone, &          ! beta
                       chi0(1,1,iw-wi+1), &  ! C
                       n & ! LDC ... leading dimension of C
                      )
         enddo
         deallocate(zm)
         allocate(zm(n,nst2,nst1))

         do iw=wi,wf
           Do ist2 = 1, nst2
             Do ist1 = 1, nst1
              zm(:,ist2,ist1)=conjg(wuo(iw,ist1,ist2)*xiuo(ist2,ist1,:))
!             zm(:,ist2,ist1)=wuo(iw,ist1,ist2)*xiuo(ist2,ist1,:)
             enddo
           enddo

           
           call zgemm('n', &           ! TRANSA = 'C'  op( A ) = A**H.
                      'n', &           ! TRANSB = 'N'  op( B ) = B.
                       n, &          ! M ... rows of op( A ) = rows of C
                       n, &           ! N ... cols of op( B ) = cols of C
                       nst1*nst2, &          ! K ... cols of op( A ) = rows of op( B )
                       zone, &          ! alpha
!                       xiuo(1,1,1), &
!                       nst1*nst2,&
!                       zm(1,1,1), &
!                       n, & 
                       zm(1,1,1), &           ! B
                       n, &          ! LDB ... leading dimension of B
                       xiuo(1,1,1), &           ! A
                       nst1*nst2,&           ! LDA ... leading dimension of A
                       zone, &          ! beta
                       chi0(1,1,iw-wi+1), &  ! C
                       n & ! LDC ... leading dimension of C
                      )
         enddo

!         chi0(:,:,:)=conjg(chi0(:,:,:))

         Call timesec (cpu1)
         cpuupd = cpuupd + cpu1 - cpu0

!write(*,*) sum(chi0(:,:,1))

deallocate(wou,wuo,zm)


endif
!do iw=1,n
!  write(*,*) chi0(iw,2,1)
!enddo
!write(*,*) dble(sum(chi0(:,:,1)))
!stop

!*****************************************************************************************************
!write(*,*) sum(chi0h)
!write(*,*) sum(chi0w)
!write(*,*) sum(chi0)
!read(*,*)
         cputot = cpuread + cpuosc + cpuupd
     ! timing information
         Call dftim (iq, ik, trim(fnxtim), cpuread, cpuosc, cpuupd, &
        & cputot)
     ! synchronize
         If ( .Not. tscreen) Call barrier
     ! end loop over k-points
      End Do
      chi0(:,:,:)=conjg(chi0(:,:,:))


      wplas = input%xs%tddft%drude(1)
      wrel = input%xs%tddft%drude(2)
      If ((wplas>1.d-8).And.(wrel>1.d-8)) then
         Do iw = wi, wf
            winv=1.0d0/(w(iw)+zi*brd)
            If (Abs(w(iw)).Lt.1.d-8) winv=1.d0
            Do oct1 = 1, 3
               chi0h (oct1, oct1, iw-wi+1) = chi0h (oct1, oct1, iw-wi+1) + wplas**2/(w(iw)+zi*wrel)*winv
            End Do
         End Do
      End If

      If (tscreen) Call ematqdealloc
  ! symmetrize head
      If (tq0) Then
         Allocate (chi0hs(3, 3, nwdfp), eps0(3, 3, nwdf))
     ! write dielectric tensor to file (unsymmetrized)
         Forall (iw=1:nwdf)
            eps0 (:, :, iw) = dble (krondelta) - chi0h (:, :, iw)
         End Forall
         If (rank .Eq. 0) Call writedielt ('DIELTENS0_NOSYM', 1, 0.d0, &
        & eps0(:, :, 1), 0)
     ! symmetrize the macroscopic dielectric function tensor
         Do oct1 = 1, 3
            Do oct2 = 1, 3
               Call symt2app (oct1, oct2, nwdfp, symt2, chi0h, &
              & chi0hs(oct1, oct2, :))
            End Do
         End Do
     ! re-assign the symmetrized head
         chi0h (:, :, :) = chi0hs (:, :, :)
     ! write dielectric tensor to file
         Forall (iw=1:nwdf)
            eps0 (:, :, iw) = dble (krondelta) - chi0hs (:, :, iw)
         End Forall
         If (rank .Eq. 0) Call writedielt ('DIELTENS0', 1, 0.d0, &
        & eps0(:, :, 1), 0)
         Deallocate (chi0hs, eps0)
      End If
  ! write response function to file
      If (tscreen) Then
     ! write out screening
         Call getunit (un)
         Open (un, File=trim(fnscreen), Form='formatted', Action='write&
        &', Status='replace')
         Call putscreen (un, tq0, n, chi0(:, :, 1), chi0h(:, :, 1), &
        & chi0w(:, :, :, 1))
         Call writevars (un, iq, 0)
         Close (un)
      Else
         Do j = 0, procs - 1
            If (rank .Eq. j) Then
               Do iw = wi, wf
                  Call putx0 (tq0, iq, iw-wi+1, trim(fnchi0_t), '', &
                 & chi0(:, :, iw-wi+1), chi0w(:, :, :, iw-wi+1), &
                 & chi0h(:, :, iw-wi+1))
               End Do
            End If
            Call barrier
         End Do
      End If
!write(*,*) sum(chi0(1:10,1:20,1))


      Deallocate (chi0, chi0h, chi0w)
      Deallocate (docc12, docc21, scis12, scis21, scis12c, scis21c)
      Deallocate (deou, deuo, wouw, wuow, wouh, wuoh, zvou, &
     & zvuo)
      Deallocate (xiou, xiuo, pmou, pmuo)
      Deallocate (bsedg)
      Deallocate (w, wreal)
      If (input%xs%tetra%tetradf) Then
         Deallocate (cw, cwa, cwsurf)
         If (input%xs%tetra%cw1k) deallocate (cwt, cw1k, cwa1k, &
        & cwsurf1k)
      End If
!      write(*,*)
End Subroutine dfq
!EOC
