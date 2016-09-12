! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine xsinit
  use modinput,only: input, isspinspiral
  use mod_names,only: revert_names
  use mod_qpoint,only: ngridq
  use mod_constants,only: h2ev
  use mod_spin,only: ncmag
  use mod_misc,only: versionname, task, githash, notelns, notes
  use modmpi,only: procs,rank,splittfile
  use modxs,only: calledxs, init0symonly, cputim0i, cntrate, &
                & systim0i, systimcum, xsfileout, fnetim, &
                & fnchi0_t, fnxtim, unitout, maxproc, &
                & lmmaxemat, lmmaxapwwf, lmmaxdielt, tordf, &
                & tscreen, nwdf, fxcdescr, fxcspin, &
                & torfxc, escale, tleblaik, tgqmaxg, &
                & tfxcbse, temat, fnresume,&
                & kset,qset,kqset,gset,gkset,gqset
  use modfxcifc,only: getfxcdata
  use m_getunit,only: getunit
  use m_genfilname,only: genfilname

  use mod_kpointset
  use mod_Gkvector, only: gkmax ! init1
  use mod_Gvector, only: intgv  ! init1
  use mod_lattice, only: bvec   ! init0
  
  implicit none

  ! local variables
  character(10) :: dat, tim
  integer :: i, un
  real(8) :: tv(3)
  real(8), parameter :: eps=1.d-7
  character(77) :: string, istring, kfname

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
  call genfilname(nodotpar=.true., basename='INFOXS', procs=procs, rank=rank, filnam=xsfileout)
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
  write(string,'("exciting ", a, " started for task ",i6)') trim(versionname), task
  call printtext(unitout,"=",string)
  if(calledxs .eq. 1) then
    write(string,'("version hash id: ",a)') githash
    call printtext(unitout,"=",string)
  end if
#ifdef MPI
  write(string,'("mpi version using ",i6," processor(s)")') procs
  call printtext(unitout,"=",string)
#ifndef MPI1   
  write(string,'("|  using mpi-2 features")') 
  call printtext(unitout,"=",string)
#endif
#endif
  if(calledxs .eq. 1) then
    if(notelns .gt. 0) then
      write(unitout,*)
      write(unitout, '("notes :")')
      do i = 1, notelns
        write(unitout, '(a)') notes(i)
      end do
    end if
  end if
  if(calledxs .eq. 1) call printtext(unitout,"=","")
  write(string,'("date (dd-mm-yyyy) : ", a2, "-", a2, "-", a4)') &
    & dat(7:8), dat(5:6), dat(1:4)
  call printtext(unitout,"=",string)
  write(string,'("time (hh:mm:ss)   : ", a2, ":", a2, ":", a2)') &
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

  !---------------------!
  !     k-point set     !
  !---------------------!
  if(any(input%xs%screening%ngridk .eq. 0)) &
    & input%xs%screening%ngridk(:) = input%groundstate%ngridk(:)
  if(any(input%xs%screening%vkloff .eq.-1.d0)) &
    & input%xs%screening%vkloff(:) = input%groundstate%vkloff(:)
  if(any(input%xs%bse%vkloff .eq.-1.d0)) &
    & input%xs%bse%vkloff(:) = input%groundstate%vkloff(:)

  !---------------------!
  !     g+k vectors     !
  !---------------------!
  if(input%xs%screening%rgkmax .eq. 0.d0) &
    & input%xs%screening%rgkmax = input%groundstate%rgkmax
  if(input%xs%bse%rgkmax .eq. 0.d0) &
    & input%xs%bse%rgkmax = input%groundstate%rgkmax

  ! Make grids

  call init0
  call init1

 ! write(*,*) "kset"
 ! write(*,*) kset
 ! write(*,*) "bvec"
 ! write(*,*) bvec
 ! write(*,*) "ngridk"
 ! write(*,*) input%groundstate%ngridk
 ! write(*,*) "vkloff"
 ! write(*,*) input%groundstate%vkloff
 ! write(*,*) "reducek"
 ! write(*,*) input%groundstate%reducek

  call generate_k_vectors(kset, bvec,&
    & input%groundstate%ngridk,&
    & input%groundstate%vkloff,&
    & input%groundstate%reducek) ! reducek is false in xs
  call getunit(un)
  open(un, file='k_set.out', action='write', status='replace')
  call print_k_vectors(kset, un)
  close(un)

  call generate_k_vectors(qset, bvec,&
    & input%xs%ngridq,&
    & [0.0d0, 0.0d0, 0.0d0],&    
    & input%xs%reduceq) ! reduceq is true in xs
  call getunit(un)
  open(un, file='q_set.out', action='write', status='replace')
  call print_k_vectors(qset, un)
  close(un)

  call generate_kq_vectors(kqset, bvec,&
    & input%groundstate%ngridk,&
    & input%groundstate%vkloff,&
    & input%groundstate%reducek)
  call getunit(un)
  open(un, file='kq_set.out', action='write', status='replace')
  call print_kq_vectors(kqset, un)
  close(un)

  call generate_G_vectors(gset, bvec,&
    & intgv,&
    & input%groundstate%gmaxvr)
  call getunit(un)
  open(un, file='g_set.out', action='write', status='replace')
  call print_G_vectors(gset, un)
  close(un)

  call generate_Gk_vectors(gkset, kset, gset, gkmax) 
  do i=1, kset%nkpt
    istring=''
    kfname=''
    write(istring,'(I12)') i
    write(kfname,'("gk_set_ik",a,".out")') trim(adjustl(istring))
    call getunit(un)
    open(un, file=kfname, action='write', status='replace')
    call print_Gk_vectors(gkset, i, un)
    close(un)
  end do
  close(un)

  call generate_Gk_vectors(gqset, qset, gset, input%xs%gqmax) 
  do i=1, qset%nkpt
    istring=''
    kfname=''
    write(istring,'(I12)') i
    write(kfname,'("gq_set_iq",a,".out")') trim(adjustl(istring))
    call getunit(un)
    open(un, file=kfname, action='write', status='replace')
    call print_Gk_vectors(gqset, i, un)
    close(un)
  end do

  !-----------------------------------!
  !     parallelization variables     !
  !-----------------------------------!
  if((procs .lt. 1) .or. (procs .gt. maxproc)) then
    write(unitout,*)
    write(unitout, '("error(xsinit): error in parallel &
      &initialization: number of processes out of range: ", i6)') procs
    write(unitout,*)
    call terminate
  end if
  if((rank .gt. procs) .or. (rank .lt. 0)) then
    write(unitout,*)
    write(unitout, '("error(xsinit): error in parallel &
      &initialization: rank out of range: ", i6)') rank
    write(unitout,*)
    call terminate
  end if

  !------------------------!
  !     spin variables     !
  !------------------------!
  ! warn for non-collinear spin polarized calculations
  if(ncmag) then
    write(unitout,*)
    write(unitout, '("warning(xsinit): calculation is spin polarized &
      &non-collinear. formalism may be incomplete.")')
    write(unitout,*)
  end if
  if(associated(input%groundstate%spin) .and. (input%xs%gqmax .gt. eps)) then
    write(unitout,*)
    write(unitout, '("warning(xsinit): spin-polarized&
    & calculation with local field effects (gqmax > 0). &
    &formalism may be incomplete.")')
    write(unitout,*)
  end if
  ! no spin-spirals
  if(isspinspiral()) then
    write(unitout,*)
    write(unitout, '("error(xsinit): xs-part not working &
      &for spin-spirals")')
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
    write(unitout, '("error(xsinit): lmaxapwwf > lmaxapw: ", i6)') input%xs%lmaxapwwf
    write(unitout,*)
    call terminate
  end if
  if(input%xs%lmaxemat .gt. input%groundstate%lmaxapw) then
    write(unitout,*)
    write(unitout, '("error(xsinit): lmaxemat > lmaxapw: ", i6)') &
      & input%xs%lmaxemat
    write(unitout,*)
    call terminate
  end if
  if(input%xs%tddft%lmaxalda .gt. input%groundstate%lmaxapw) then
    write(unitout,*)
    write(unitout, '("error(xsinit): lmaxalda > lmaxapw: ", i6)') &
      & input%xs%tddft%lmaxalda
    write(unitout,*)
    call terminate
  end if
  if(input%xs%lmaxemat .gt. input%xs%lmaxapwwf) then
     write(unitout,*)
     write(unitout, '("warning(xsinit): lmaxemat > lmaxapwwf: ", i6)') input%xs%lmaxemat
     write(unitout,*)
  end if


  !------------------------------------!
  !     secular equation variables     !
  !------------------------------------!
  if(input%xs%screening%nempty .eq. 0) &
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
    write(unitout, '("error(xsinit): tetrahedron method does not work &
     &in  combination with analytic continuation")')
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
  if(input%xs%tddft%tordfxc) torfxc = - 1.d0

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
  if((input%xs%xstype.eq."tddft").and.(input%xs%gqmaxtype.eq."|g|")) tgqmaxg=.true.
  tfxcbse = .false.
  if(input%xs%tddft%fxctypenumber .eq. 5) then
     if(input%groundstate%gmaxvr .lt. 2.d0*input%xs%gqmax) then
        write(unitout,*)
        write(unitout, '("error(xsinit): 2*gqmax > gmaxvr", 2g18.10)') &
         & 2.d0 * input%xs%gqmax, input%groundstate%gmaxvr
        write(unitout,*)
        call terminate
     end if
  end if

  if((input%xs%tddft%fxctypenumber .eq. 7) .or. &
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
    write(unitout, '("info(xsinit): mapping screening-specific parameters")')
    write(unitout,*)
    tv(:) = dble(input%xs%screening%ngridk(:)) / dble(ngridq(:))
    tv(:) = tv(:) - int(tv(:))
    if(sum(tv) .gt. input%structure%epslat) then
      write(unitout,*)
      write(unitout, '("error(xsinit): ngridkscr must be an integer multiple of ngridq")')
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
    write(unitout, '("info(xsinit): mapping bse-specific parameters")')
    write(unitout,*)
    if(any(input%groundstate%ngridk .ne. ngridq)) then
      write(unitout,*)
      write(unitout, '("error(xsinit): ngridk must be equal ngridq for the bse-hamiltonian")')
      write(unitout, '(" ngridk : ", 3i6)') input%groundstate%ngridk
      write(unitout, '(" ngridq : ", 3i6)') ngridq
      write(unitout,*)
      call terminate
    end if
  end if

  temat=.true.
  if(input%xs%xstype.eq."tddft") then
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
    call genfilname(basename='RESUME', rank=rank, procs=procs, dotext='', filnam=fnresume)
  else
    call genfilname(basename='RESUME', dotext='', filnam=fnresume)
  end if
  ! check for stale checkpoint file
  call chkptchk

  ! define checkpoint
  call chkpt(1, (/ task /), 'passed xsinit')

end subroutine xsinit
