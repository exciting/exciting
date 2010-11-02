! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xsinit
      Use modinput
      Use modmain
      Use modmpi
      Use modxs
      Use modfxcifc
      Use m_getunit
      Use m_genfilname
      Implicit None
  ! local variables
      Character (10) :: dat, tim
      Integer :: i
      Real (8) :: tv (3)
!
  !---------------------------!
  !     initialize timing     !
  !---------------------------!
  ! remember how often this routine is called
      calledxs = calledxs + 1
  ! only recalculate symmetries in init0
      If (calledxs .Gt. 1) init0symonly = .True.
  ! initialize global counters
      Call cpu_time (cputim0i)
      Call system_clock (COUNT_RATE=cntrate)
      Call system_clock (COUNT=systim0i)
      Call date_and_time (date=dat, time=tim)
      If (calledxs .Eq. 1) Call system_clock (COUNT=systimcum)
!
  !---------------------!
  !     output file     !
  !---------------------!
  ! name of output file
      Call genfilname (nodotpar=.True., basename='INFOXS', procs=procs, &
     & rank=rank, filnam=xsfileout)
  ! reset or append to output file
      Call getunit (unitout)
      If (input%xs%tappinfo .Or. (calledxs .Gt. 1)) Then
         Open (unitout, File=trim(xsfileout), Action='write', &
        & Position='append')
      Else
         Open (unitout, File=trim(xsfileout), Action='write', Status='r&
        &eplace')
      End If
  ! write to info file
      If (calledxs .Eq. 1) Then
         Write (unitout,*)
         Write (unitout, '("+-----------------------------------------------------------+")')
         Write (unitout, '("| EXCITING helium    (",I2.2,".",I2.2,".",I2.2,") started                     |")') version
         Write (unitout, '("| version hash id: ",a," |")') githash
#ifdef MPI
         Write (unitout, '("| MPI version using ",i6," processor(s)                     |")') procs
         if (rank .ne. 0) &
         Write (unitout, '("|  rank of current processor: ",i6,"                        |")') rank
#ifndef MPI1
         Write (unitout, '("|  using MPI-2 features                                     |")')
#endif
#endif
         Write (unitout, '("+ ----------------------------------------------------------+")')
         If (notelns .Gt. 0) Then
            Write (unitout,*)
            Write (unitout, '("Notes :")')
            Do i = 1, notelns
               Write (unitout, '(A)') notes (i)
            End Do
         End If
      End If
      Write (unitout,*)
      Write (unitout, '("Date (YYYY-MM-DD) : ", A4, "-", A2, "-", A2)') &
     & dat (1:4), dat (5:6), dat (7:8)
      Write (unitout, '("Time (hh:mm:ss)   : ", A2, ":", A2, ":", A2)') &
     & tim (1:2), tim (3:4), tim (5:6)
      Write (unitout, '("Info(xsinit): task Nr.", i6, " started")') &
     & task
      Call flushifc (unitout)
!
  !--------------------------------------------!
  !     map xs parameters associated to GS     !
  !--------------------------------------------!
      Call mapxsparameters
!
!
  !-----------------------------------!
  !     parallelization variables     !
  !-----------------------------------!
      If ((procs .Lt. 1) .Or. (procs .Gt. maxproc)) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): Error in parallel initializa&
        &tion: number of processes out of range: ", i6)') procs
         Write (unitout,*)
         Call terminate
      End If
      If ((rank .Gt. procs) .Or. (rank .Lt. 0)) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): Error in parallel initializa&
        &tion: rank out of range: ", i6)') rank
         Write (unitout,*)
         Call terminate
      End If
!
  !------------------------!
  !     spin variables     !
  !------------------------!
  ! warn for spin polarized calculations
      If (associated(input%groundstate%spin)) Then
         Write (unitout,*)
         Write (unitout, '("Warning(xsinit): calculation is spin-polari&
        &zed - formalism may be incomplete")')
         Write (unitout,*)
      End If
  ! no spin-spirals
      If (isspinspiral()) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): xs-part not working for spin&
        &-spirals")')
         Write (unitout,*)
         Call terminate
      End If
!
  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
      If (input%xs%lmaxapwwf .Eq.-1) input%xs%lmaxapwwf = &
     & input%groundstate%lmaxmat
      lmmaxapwwf = (input%xs%lmaxapwwf+1) ** 2
      lmmaxemat = (input%xs%lmaxemat+1) ** 2
      lmmaxdielt = (input%xs%BSE%lmaxdielt+1) ** 2
      If (input%xs%lmaxapwwf .Gt. input%groundstate%lmaxapw) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): lmaxapwwf > lmaxapw: ", i6)') input%xs%lmaxapwwf
         Write (unitout,*)
         Call terminate
      End If
      If (input%xs%lmaxemat .Gt. input%groundstate%lmaxapw) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): lmaxemat > lmaxapw: ", i6)') &
        & input%xs%lmaxemat
         Write (unitout,*)
         Call terminate
      End If
      If (input%xs%tddft%lmaxalda .Gt. input%groundstate%lmaxapw) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): lmaxalda > lmaxapw: ", i6)') &
        & input%xs%tddft%lmaxalda
         Write (unitout,*)
         Call terminate
      End If
      If (input%xs%lmaxemat .Gt. input%xs%lmaxapwwf) Then
         Write (unitout,*)
         Write (unitout, '("Warning(xsinit): lmaxemat > lmaxapwwf: ", i&
        &6)') input%xs%lmaxemat
         Write (unitout,*)
      End If
!
  !---------------------!
  !     k-point set     !
  !---------------------!
      If (any(input%xs%screening%ngridk .Eq. 0)) &
     & input%xs%screening%ngridk(:) = input%groundstate%ngridk(:)
      If (any(input%xs%screening%vkloff .Eq.-1.d0)) &
     & input%xs%screening%vkloff(:) = input%groundstate%vkloff(:)
      If (any(input%xs%BSE%vkloff .Eq.-1.d0)) input%xs%BSE%vkloff(:) = &
     & input%groundstate%vkloff(:)
!
  !---------------------!
  !     G+k vectors     !
  !---------------------!
      If (input%xs%screening%rgkmax .Eq. 0.d0) &
     & input%xs%screening%rgkmax = input%groundstate%rgkmax
      If (input%xs%BSE%rgkmax .Eq. 0.d0) input%xs%BSE%rgkmax = &
     & input%groundstate%rgkmax
!
  !------------------------------------!
  !     secular equation variables     !
  !------------------------------------!
      If (input%xs%screening%nempty .Eq. 0) input%xs%screening%nempty = &
     & input%groundstate%nempty
  ! set splittfile parameter for splitting of eigenvector files in
  ! parallelization of SCF cycle
      If ((task .Ne. 301) .And. (task .Ne. 401)) splittfile = .False.
!
  !----------------------------!
  !     response functions     !
  !----------------------------!
  ! set time-ordering
      tordf = 1.d0
      If (input%xs%tddft%torddf) tordf = - 1.d0
      tscreen = .False.
      If ((task .Ge. 400) .And. (task .Le. 499)) tscreen = .True.
  ! tetrahedron method not implemented for analytic continuation
      If (input%xs%tetra%tetradf .And. input%xs%tddft%acont) Then
         Write (unitout,*)
         Write (unitout, '("Error(xsinit): tetrahedron method does not &
        &work in  combination with analytic continuation")')
         Write (unitout,*)
         Call terminate
      End If
  ! if imaginary frequencies intervals are not specified
 ! if (input%xs%tddft%nwacont.eq.0) input%xs%tddft%nwacont=input%xs%energywindow%points
 ! nwdf=input%xs%energywindow%points
!
      If (input%xs%tddft%acont) Then
         nwdf = input%xs%tddft%nwacont
      Else
         nwdf = input%xs%energywindow%points
      End If
  ! get exchange-correlation kernel functional data
      Call getfxcdata (input%xs%tddft%fxctypenumber, fxcdescr, fxcspin)
!
  !-----------------------------!
  !     xc-kernel variables     !
  !-----------------------------!
  ! set time-ordering
      torfxc = 1.d0
      If (input%xs%tddft%tordfxc) torfxc = - 1.d0
!
  !-----------------------!
  !     miscellaneous     !
  !-----------------------!
  ! scaling factor for output of energies
      escale = 1.d0
      If (input%xs%tevout) escale = h2ev
      tleblaik = .True.
      If (input%xs%BSE%nleblaik .Eq. 0) tleblaik = .False.
!
  !----------------------------------!
  !     task dependent variables     !
  !----------------------------------!
      tgqmaxg=.false.
      if ((input%xs%xstype.eq."TDDFT").and.(input%xs%gqmaxtype.eq."|G|")) tgqmaxg=.true.
      tfxcbse = .False.
      If (input%xs%tddft%fxctypenumber .Eq. 5) Then
         If (input%groundstate%gmaxvr .Lt. 2.d0*input%xs%gqmax) Then
            Write (unitout,*)
            Write (unitout, '("Error(xsinit): 2*gqmax > gmaxvr", 2g18.1&
           &0)') 2.d0 * input%xs%gqmax, input%groundstate%gmaxvr
            Write (unitout,*)
            Call terminate
         End If
      End If
      If ((input%xs%tddft%fxctypenumber .Eq. 7) .Or. &
     & (input%xs%tddft%fxctypenumber .Eq. 8)) tfxcbse = .True.
      If ((task .Ge. 401) .And. (task .Le. 439)) Then
     ! screening
         input%groundstate%nosym = input%xs%screening%nosym
         input%groundstate%reducek = input%xs%screening%reducek
         input%groundstate%rgkmax = input%xs%screening%rgkmax
         input%groundstate%nempty = input%xs%screening%nempty
         input%groundstate%ngridk (:) = input%xs%screening%ngridk(:)
         input%groundstate%vkloff (:) = input%xs%screening%vkloff(:)
         Write (unitout,*)
         Write (unitout, '("Info(xsinit): mapping screening-specific pa&
        &rameters")')
         Write (unitout,*)
         tv (:) = dble (input%xs%screening%ngridk(:)) / dble &
        & (ngridq(:))
         tv (:) = tv (:) - Int (tv(:))
         If (sum(tv) .Gt. input%structure%epslat) Then
            Write (unitout,*)
            Write (unitout, '("Error(xsinit): ngridkscr must be an inte&
           &ger multiple of ngridq")')
            Write (unitout, '(" ngridkscr : ", 3i6)') &
           & input%xs%screening%ngridk
            Write (unitout, '(" ngridq    : ", 3i6)') ngridq
            Write (unitout,*)
            Call terminate
         End If
      Else If ((task .Ge. 440) .And. (task .Le. 459)) Then
     ! BSE
         input%groundstate%nosym = input%xs%BSE%nosym
         input%groundstate%reducek = input%xs%BSE%reducek
         input%groundstate%rgkmax = input%xs%BSE%rgkmax
         input%groundstate%vkloff (:) = input%xs%BSE%vkloff(:)
         ngridq (:) = input%xs%ngridq(:)
         Write (unitout,*)
         Write (unitout, '("Info(xsinit): mapping BSE-specific paramete&
        &rs")')
         Write (unitout,*)
         If (any(input%groundstate%ngridk .Ne. ngridq)) Then
            Write (unitout,*)
            Write (unitout, '("Error(xsinit): ngridk must be equal ngri&
           &dq for the BSE-Hamiltonian")')
            Write (unitout, '(" ngridk : ", 3i6)') &
           & input%groundstate%ngridk
            Write (unitout, '(" ngridq : ", 3i6)') ngridq
            Write (unitout,*)
            Call terminate
         End If
      End If
!
  !--------------------!
  !     file names     !
  !--------------------!
  ! revert file names to default
      Call revert_names
!
  !---------------------!
  !     checkpoints     !
  !---------------------!
      If (procs .Gt. 1) Then
         Call genfilname (basename='resume', rank=rank, procs=procs, &
        & dotext='', filnam=fnresume)
      Else
         Call genfilname (basename='resume', dotext='', &
        & filnam=fnresume)
      End If
  ! check for stale checkpoint file
      Call chkptchk
!
  ! define checkpoint
      Call chkpt (1, (/ task /), 'passed xsinit')
End Subroutine xsinit
