! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine xsinit(j, plan)
  use modinput, only: input, isspinspiral, plan_type
  use mod_names,only: revert_names
  use mod_qpoint,only: ngridq
  use mod_constants,only: h2ev
  use mod_spin,only: ncmag
  use mod_misc,only: versionname, githash, notelns, notes
  use modmpi
  use modxs,only: calledxs, init0symonly, cputim0i, cntrate, &
                & systim0i, systimcum, xsfileout, fnetim, &
                & fnchi0_t, fnxtim, unitout, maxproc, &
                & lmmaxemat, lmmaxapwwf, lmmaxdielt, tordf, &
                & tscreen, nwdf, fxcdescr, fxcspin, &
                & torfxc, escale, tleblaik, tgqmaxg, &
                & tfxcbse, temat, fnresume
  use modfxcifc,only: getfxcdata
  use m_getunit,only: getunit
  use m_genfilname,only: genfilname

  implicit none
  integer, intent(in) :: j
  type(plan_type), intent(in) :: plan
  ! local variables
  character(10) :: dat, tim
  character(500) :: taskname
  integer :: i, task
  real(8) :: tv(3)
  real(8), parameter :: eps=1.d-7
  character(77) :: string

  !--------------------------------!
  !     set taskname and number    !
  !--------------------------------!
  task=plan%doonlyarray(j)%doonly%tasknumber
  taskname=plan%doonlyarray(j)%doonly%task
  
  ! Backups of groundstate variables
  call backup0
  call backup1
  call backup2

  !---------------------------!
  !     initialize timing     !
  !---------------------------!
  ! remember how often this routine is called
  calledxs = calledxs + 1
  ! only recalculate symmetries in init0
  if(calledxs .gt. 1) init0symonly = .true.
  ! initialize global counters
  call cpu_time(cputim0i)
  call system_clock(count_rate=cntrate)
  call system_clock(count=systim0i)
  call date_and_time(date=dat, time=tim)
  if(calledxs .eq. 1) call system_clock(count=systimcum)

  !---------------------!
  !     output file     !
  !---------------------!
  ! name of output file
  call genfilname(nodotpar=.true., basename='INFOXS', procs=procs,&
    & rank=rank, filnam=xsfileout)
  call genfilname(basename='EMAT_TIMING', procs=procs, rank=rank, filnam=fnetim)
  call genfilname(basename='X0', procs=procs, rank=rank, filnam=fnchi0_t)
  call genfilname(basename='X0_TIMING', procs=procs, rank=rank, filnam=fnxtim)

  ! reset or append to output file
  call getunit(unitout)
  if(input%xs%tappinfo .or. (calledxs .gt. 1)) then
    open(unitout, file=trim(xsfileout), action='write', position='append')
  else
    open(unitout, file=trim(xsfileout), action='write', status='replace')
  end if
  
  ! write to info file
  call printline(unitout,"=")
  write(string,'("EXCITING ", a, " started for task ",a," (",i3,")")') trim(versionname), trim(taskname), task
  call printtext(unitout,"=",string)
  if(calledxs .eq. 1) then
    write(string,'("version hash id: ",a)') githash
    call printtext(unitout,"=",string)
  end if
#ifdef MPI
  write(string,'("MPI version using ",i6," processor(s)")') procs
  call printtext(unitout,"=",string)
#ifndef MPI1   
  write(string,'("|  using MPI-2 features")') 
  call printtext(unitout,"=",string)
#endif
#endif
  if(calledxs .eq. 1) then
    if(notelns .gt. 0) then
      write(unitout,*)
      write(unitout, '("Notes :")')
      do i = 1, notelns
        write(unitout, '(a)') notes(i)
      end do
    end if
  end if
  if(calledxs .eq. 1) call printtext(unitout,"=","")
  write(string,'("Date (DD-MM-YYYY) : ", a2, "-", a2, "-", a4)')&
    & dat(7:8), dat(5:6), dat(1:4)
  call printtext(unitout,"=",string)
  write(string,'("Time (hh:mm:ss)   : ", a2, ":", a2, ":", a2)')&
    & tim(1:2), tim(3:4), tim(5:6)
  call printline(unitout,"=")
  call flushifc(unitout)

  !--------------------------------------------!
  !     map xs parameters associated to gs     !
  !--------------------------------------------!

  if(input%xs%rgkmax .eq. 0.d0) input%xs%rgkmax = input%groundstate%rgkmax

  ! This sets input%groundstate%*=input%xs%*, where * are the following
  ! nosym, ngridk, reducek, vkloff, maxscl
  ! If phonons also: rgkmax, swidth, lmaxapw, lmaxmat, nempty
  ! If spin: bfieldc
  call mapxsparameters

  !-----------------------------------!
  !     parallelization variables     !
  !-----------------------------------!
  if((procs .lt. 1) .or. (procs .gt. maxproc)) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): Error in parallel &
      &initialization: number of processes out of range: ", i6)') procs
    write(unitout,*)
    call terminate
  end if
  if((rank .gt. procs) .or. (rank .lt. 0)) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): Error in parallel &
      &initialization: rank out of range: ", i6)') rank
    write(unitout,*)
    call terminate
  end if

  !------------------------!
  !     spin variables     !
  !------------------------!
  ! no spin-spirals
  if(isspinspiral()) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): xs-part not working&
      & for spin-spirals")')
    write(unitout,*)
    call terminate
  end if

  !-----------------------------!
  !     Core Non-TDA calc.s     !
  !-----------------------------!
  if((input%xs%BSE%xas) .and. (input%xs%BSE%coupling)) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): Calculations of Core BSE&
      & spectra beyond the Tamm-Dancoff approximation not implemented yet. &
      &Please contact the developers at exciting-code.org")')
    write(unitout,*)
    call terminate
  end if

  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
  if(input%xs%lmaxapwwf .eq.-1) input%xs%lmaxapwwf = input%groundstate%lmaxmat
  lmmaxapwwf = (input%xs%lmaxapwwf+1) ** 2
  lmmaxemat = (input%xs%lmaxemat+1) ** 2
  lmmaxdielt = (input%xs%bse%lmaxdielt+1) ** 2
  if(input%xs%lmaxapwwf .gt. input%groundstate%lmaxapw) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): lmaxapwwf > lmaxapw: ", i6)') input%xs%lmaxapwwf
    write(unitout,*)
    call terminate
  end if
  if(input%xs%lmaxemat .gt. input%groundstate%lmaxapw) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): lmaxemat > lmaxapw: ", i6)')&
      & input%xs%lmaxemat
    write(unitout,*)
    call terminate
  end if
  if(input%xs%tddft%lmaxalda .gt. input%groundstate%lmaxapw) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): lmaxalda > lmaxapw: ", i6)')&
      & input%xs%tddft%lmaxalda
    write(unitout,*)
    call terminate
  end if
  if(input%xs%lmaxemat .gt. input%xs%lmaxapwwf) then
     write(unitout,*)
     write(unitout, '("Warning(xsinit): lmaxemat > lmaxapwwf: ", i6)') input%xs%lmaxemat
     write(unitout,*)
  end if

  !---------------------!
  !     k-point set     !
  !---------------------!
  if(any(input%xs%screening%ngridk .eq. 0))&
    & input%xs%screening%ngridk(:) = input%groundstate%ngridk(:)
  if(any(input%xs%screening%vkloff .eq. -1.d0))&
    & input%xs%screening%vkloff(:) = input%groundstate%vkloff(:)
  if(any(input%xs%bse%vkloff .eq. -1.d0))&
    & input%xs%bse%vkloff(:) = input%groundstate%vkloff(:)

  !---------------------!
  !     G+k vectors     !
  !---------------------!
  if(input%xs%screening%rgkmax .eq. 0.d0)&
    & input%xs%screening%rgkmax = input%groundstate%rgkmax
  if(input%xs%bse%rgkmax .eq. 0.d0)&
    & input%xs%bse%rgkmax = input%groundstate%rgkmax

  !------------------------------------!
  !     secular equation variables     !
  !------------------------------------!
  if(input%xs%screening%nempty .eq. 0)&
    & input%xs%screening%nempty = input%groundstate%nempty
  ! set splittfile parameter for splitting of eigenvector files in
  ! parallelization of scf cycle
  ! Task 301 corresponds to "xsgeneigvec" plan
  ! Task 401 corresponds to "scrgeneigvec" plan
  if((task .ne. 301) .and. (task .ne. 401)) splittfile = .false.

  !----------------------------!
  !     response functions     !
  !----------------------------!
  ! set time-ordering
  tordf = 1.d0
  if(input%xs%tddft%torddf) tordf = - 1.d0
  tscreen = .false.
  ! Any task related to screening and bse is in the 400 range
  if((task .ge. 400) .and. (task .le. 499)) tscreen = .true.
  ! tetrahedron method not implemented for analytic continuation
  if(input%xs%tetra%tetradf .and. input%xs%tddft%acont) then
    write(unitout,*)
    write(unitout, '("Error(xsinit): tetrahedron method does not work&
     & in combination with analytic continuation")')
    write(unitout,*)
    call terminate
  end if

  if(input%xs%tddft%acont) then
    nwdf = input%xs%tddft%nwacont
  else
    nwdf = input%xs%energywindow%points
  end if

  ! get exchange-correlation kernel functional data
  call getfxcdata(input%xs%tddft%fxctypenumber, fxcdescr, fxcspin)

  !-----------------------------!
  !     xc-kernel variables     !
  !-----------------------------!
  ! set time-ordering
  torfxc = 1.d0
  if(input%xs%tddft%tordfxc) torfxc = -1.d0

  !-----------------------!
  !     miscellaneous     !
  !-----------------------!
  ! scaling factor for output of energies
  escale = 1.d0
  if(input%xs%tevout) escale = h2ev
  tleblaik = .true.
  if(input%xs%bse%nleblaik .eq. 0) tleblaik = .false.

  !----------------------------------!
  !     task dependent variables     !
  !----------------------------------!
  tgqmaxg=.false.
  if((input%xs%xstype.eq."TDDFT") .and. (input%xs%gqmaxtype.eq."|G|")) tgqmaxg=.true.
  tfxcbse = .false.
  if(input%xs%tddft%fxctypenumber .eq. 5) then
    if(input%groundstate%gmaxvr .lt. 2.d0*input%xs%gqmax) then
      write(unitout,*)
      write(unitout, '("Error(xsinit): 2*gqmax > gmaxvr", 2g18.10)')&
       & 2.d0 * input%xs%gqmax, input%groundstate%gmaxvr
      write(unitout,*)
      call terminate
    end if
  end if

  if((input%xs%tddft%fxctypenumber .eq. 7) .or.&
    & (input%xs%tddft%fxctypenumber .eq. 8)) tfxcbse = .true.

  if((task .ge. 401) .and. (task .le. 439)) then
    ! screening
    input%groundstate%nosym = input%xs%screening%nosym
    input%groundstate%reducek = input%xs%screening%reducek
    input%groundstate%rgkmax = input%xs%screening%rgkmax
    input%groundstate%nempty = input%xs%screening%nempty
    input%groundstate%ngridk(:) = input%xs%screening%ngridk(:)
    input%groundstate%vkloff(:) = input%xs%screening%vkloff(:)
    write(unitout,*)
    write(unitout, '("Info(xsinit): mapping screening-specific parameters")')
    write(unitout,*)
    tv(:) = dble(input%xs%screening%ngridk(:)) / dble(ngridq(:))
    tv(:) = tv(:) - int(tv(:))
    if(sum(tv) .gt. input%structure%epslat) then
      write(unitout,*)
      write(unitout, '("Error(xsinit): ngridkscr must be an&
        & integer multiple of ngridq")')
      write(unitout, '(" ngridkscr : ", 3i6)') input%xs%screening%ngridk
      write(unitout, '(" ngridq    : ", 3i6)') ngridq
      write(unitout,*)
      call terminate
    end if
  else if((task .ge. 440) .and. (task .le. 459)) then
    ! bse
    input%groundstate%nosym = input%xs%bse%nosym
    input%groundstate%reducek = input%xs%bse%reducek
    input%groundstate%rgkmax = input%xs%bse%rgkmax
    input%groundstate%vkloff(:) = input%xs%bse%vkloff(:)
    ngridq(:) = input%xs%ngridq(:)
    write(unitout,*)
    write(unitout, '("Info(xsinit): mapping BSE-specific parameters")')
    write(unitout,*)
    if(any(input%groundstate%ngridk .ne. ngridq)) then
      write(unitout,*)
      write(unitout, '("Error(xsinit): ngridk must be equal ngridq&
        & for the BSE-hamiltonian")')
      write(unitout, '(" ngridk : ", 3i6)') input%groundstate%ngridk
      write(unitout, '(" ngridq : ", 3i6)') ngridq
      write(unitout,*)
      call terminate
    end if
  end if

  temat=.true.
  if(input%xs%xstype.eq."TDDFT") then
    if((size(input%xs%qpointset%qpoint, 2).eq.1).and.(input%xs%gqmax.lt.eps)) then
      if(sum(abs(input%xs%qpointset%qpoint(:, 1))) .lt. eps) then 
        temat = .false.
      end if
    end if
  end if

  !--------------------!
  !     file names     !
  !--------------------!
  ! revert file names to default
  call revert_names

  !---------------------!
  !     checkpoints     !
  !---------------------!
  if(procs .gt. 1) then
    call genfilname(basename='resume', rank=rank, procs=procs,&
      & dotext='', filnam=fnresume)
  else
    call genfilname(basename='resume', dotext='', filnam=fnresume)
  end if
  ! check for stale checkpoint file
  call chkptchk
  ! define checkpoint
  call chkpt(1, (/ task /), 'passed xsinit')

  ! Some xas specific init
  if(input%xs%bse%xas) call xasinit

end subroutine xsinit
