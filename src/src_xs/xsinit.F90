



! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine xsinit
#include "../version.inc"
use modinput
  use modmain
  use modmpi
  use modxs
  use modfxcifc
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(10)::dat, tim
  integer :: i
  real(8) :: tv(3)

  !---------------------------!
  !     initialize timing     !
  !---------------------------!
  ! remember how often this routine is called
  calledxs=calledxs+1
  ! only recalculate symmetries in init0
  if (calledxs.gt.1) init0symonly=.true.
  ! initialize global counters
  call cpu_time(cputim0i)
  call system_clock(COUNT_RATE=cntrate)
  call system_clock(COUNT=systim0i)
  call date_and_time(date=dat, time=tim)
  if (calledxs.eq.1) call system_clock(COUNT=systimcum)

  !---------------------!
  !     output file     !
  !---------------------!
  ! set version of XS part
  call xssetversion
  ! name of output file
  call genfilname(nodotpar = .true., basename = 'INFOXS', procs = procs, rank = rank, &
       filnam = xsfileout)
  ! reset or append to output file
  call getunit(unitout)
  if (input%xs%tappinfo.or.(calledxs.gt.1)) then
     open(unitout, file=trim(xsfileout), action='write', position='append')
  else
     open(unitout, file=trim(xsfileout), action='write', status='replace')
  end if
  ! write to info file
  if (calledxs.eq.1) then
     write(unitout, *)
     write(unitout, '("+-------------------------------------------------------&
	  &---+")')
     write(unitout, '("| EXCITING version ", I1.1, ".", I1.1, ".", I3.3, " (eXcited &
	  &States ", I1.1, ".", I3.3, ") started  |")') version, versionxs
     write(unitout, '("| git hash id : ", a20, "		       |")') GITHASH
#ifdef LOCALCHG
     write(unitout, '("| Warning     : source codes deviates from the git hash &
	  &id |")')
#endif
#ifdef MPI
     write(unitout, '("| compiled for MPI execution			       &
	  &   |")')
#endif
#ifndef MPI
     write(unitout, '("| compiled for serial execution			       &
	  &   |")')
#endif
     write(unitout, '("+ ------------------------------------------------------&
	  &---+")')
     if ((procs.gt.1).and.(rank.eq.0)) write(unitout, '("(parallel) master, &
	&rank/number of processes:", 2i8)') rank, procs
     if ((procs.gt.1).and.(rank.ne.0)) write(unitout, '("(parallel) slave,  &
	&rank/number of processes:", 2i8)') rank, procs
     if (notelns.gt.0) then
	write(unitout, *)
	write(unitout, '("Notes :")')
	do i=1, notelns
	   write(unitout, '(A)') notes(i)
	end do
     end if
  end if
  write(unitout, *)
  write(unitout, '("Date (YYYY-MM-DD) : ", A4, "-", A2, "-", A2)') dat(1:4), &
	dat(5:6), dat(7:8)
  write(unitout, '("Time (hh:mm:ss)   : ", A2, ":", A2, ":", A2)') tim(1:2), &
	tim(3:4), tim(5:6)
  write(unitout, '("Info(xsinit): task Nr.", i6, " started")') task
  call flushifc(unitout)

  !--------------------------------------------!
  !     map xs parameters associated to GS     !
  !--------------------------------------------!
  call mapxsparameters


  !-----------------------------------!
  !     parallelization variables     !
  !-----------------------------------!
  if ((procs.lt.1).or.(procs.gt.maxproc)) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): Error in parallel initialization: number &
     &of processes out of range: ", i6)') procs
     write(unitout, *)
     call terminate
  end if
  if ((rank.gt.procs).or.(rank.lt.0)) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): Error in parallel initialization: rank &
     &out of range: ", i6)') rank
     write(unitout, *)
     call terminate
  end if

  !------------------------!
  !     spin variables     !
  !------------------------!
  ! warn for spin polarized calculations
  if (associated(input%groundstate%spin)) then
     write(unitout, *)
     write(unitout, '("Warning(xsinit): calculation is spin-polarized - &
     &formalism may be incomplete")')
     write(unitout, *)
  end if
  ! no spin-spirals
  if (isspinspiral()) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): xs-part not working for spin-spirals")')
     write(unitout, *)
     call terminate
  end if

  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
  if (input%xs%lmaxapwwf.eq.-1) input%xs%lmaxapwwf=input%groundstate%lmaxmat
  lmmaxapwwf=(input%xs%lmaxapwwf+1)**2
  lmmaxemat=(input%xs%lmaxemat+1)**2
  lmmaxdielt=(input%xs%BSE%lmaxdielt+1)**2
  if (input%xs%lmaxapwwf.gt.input%groundstate%lmaxapw) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): lmaxapwwf > lmaxapw: ", i6)') input%xs%lmaxapwwf
     write(unitout, *)
     call terminate
  end if
  if (input%xs%lmaxemat.gt.input%groundstate%lmaxapw) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): lmaxemat > lmaxapw: ", i6)') input%xs%lmaxemat
     write(unitout, *)
     call terminate
  end if
  if (input%xs%tddft%lmaxalda.gt.input%groundstate%lmaxapw) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): lmaxalda > lmaxapw: ", i6)') input%xs%tddft%lmaxalda
     write(unitout, *)
     call terminate
  end if
  if (input%xs%lmaxemat.gt.input%xs%lmaxapwwf) then
     write(unitout, *)
     write(unitout, '("Warning(xsinit): lmaxemat > lmaxapwwf: ", i6)') input%xs%lmaxemat
     write(unitout, *)
  end if

  !---------------------!
  !     k-point set     !
  !---------------------!
  if (any(input%xs%screening%ngridk.eq.0)) input%xs%screening%ngridk(:)=input%groundstate%ngkgrid(:)
  if (any(input%xs%screening%vkloff.eq. - 1.d0)) input%xs%screening%vkloff(:) = input%groundstate%vkloff(:)
  if (any(input%xs%BSE%vkloff.eq.-1.d0)) input%xs%BSE%vkloff(:)=input%groundstate%vkloff(:)

  !---------------------!
  !     G+k vectors     !
  !---------------------!
  if (input%xs%screening%rgkmax.eq.0.d0) input%xs%screening%rgkmax=input%groundstate%rgkmax
  if (input%xs%BSE%rgkmax.eq.0.d0) input%xs%BSE%rgkmax=input%groundstate%rgkmax

  !------------------------------------!
  !     secular equation variables     !
  !------------------------------------!
  if (input%xs%screening%nempty.eq.0) input%xs%screening%nempty=input%groundstate%nempty
  ! set splittfile parameter for splitting of eigenvector files in
  ! parallelization of SCF cycle
  if ((task.ne.301).and.(task.ne.401)) splittfile=.false.

  !----------------------------!
  !     response functions     !
  !----------------------------!
  ! set time-ordering
  tordf=1.d0
  if (input%xs%tddft%torddf) tordf=-1.d0
  tscreen=.false.
  if ((task.ge.400).and.(task.le.499)) tscreen=.true.
  ! tetrahedron method not implemented for analytic continuation
  if (input%xs%tetra%tetradf.and.input%xs%tddft%acont) then
     write(unitout, *)
     write(unitout, '("Error(xsinit): tetrahedron method does not work in &
	& combination with analytic continuation")')
     write(unitout, *)
     call terminate
  end if
  ! if imaginary frequencies intervals are not specified
  if (input%xs%tddft%nwacont.eq.0) input%xs%tddft%nwacont=input%properties%dos%nwdos
  nwdf=input%properties%dos%nwdos
  if (input%xs%tddft%acont) nwdf=input%xs%tddft%nwacont
  ! get exchange-correlation kernel functional data
  call getfxcdata(input%xs%tddft%fxctypenumber, fxcdescr, fxcspin)

  !-----------------------------!
  !     xc-kernel variables     !
  !-----------------------------!
  ! set time-ordering
  torfxc=1.d0
  if (input%xs%tddft%tordfxc) torfxc=-1.d0

  !-----------------------!
  !     miscellaneous     !
  !-----------------------!
  ! scaling factor for output of energies
  escale=1.d0
  if (input%xs%tevout) escale=27.2114d0
  tleblaik=.true.
  if (input%xs%BSE%nleblaik.eq.0) tleblaik=.false.

  !----------------------------------!
  !     task dependent variables     !
  !----------------------------------!
  tfxcbse=.false.
  if (input%xs%tddft%fxctypenumber.eq.5) then
     if (input%groundstate%gmaxvr.lt.2.d0*input%xs%gqmax) then
	write(unitout, *)
	write(unitout, '("Error(xsinit): 2*gqmax > gmaxvr", 2g18.10)') 2.d0* &
	  input%xs%gqmax, input%groundstate%gmaxvr
	write(unitout, *)
	call terminate
     end if
  end if
  if ((input%xs%tddft%fxctypenumber.eq.7).or.(input%xs%tddft%fxctypenumber.eq.8)) tfxcbse=.true.
  if ((task.ge.401).and.(task.le.439)) then
     ! screening
     input%groundstate%nosym=input%xs%screening%nosym
     input%groundstate%reducek=input%xs%screening%reducek
     input%groundstate%rgkmax=input%xs%screening%rgkmax
     input%groundstate%nempty=input%xs%screening%nempty
     input%groundstate%ngkgrid(:)=input%xs%screening%ngridk(:)
     input%groundstate%vkloff(:)=input%xs%screening%vkloff(:)
     write(unitout, *)
     write(unitout, '("Info(xsinit): mapping screening-specific parameters")')
     write(unitout, *)
     tv(:)=dble(input%xs%screening%ngridk(:))/dble(ngridq(:))
     tv(:)=tv(:)-int(tv(:))
     if (sum(tv).gt.input%structure%epslat) then
	write(unitout, *)
	write(unitout, '("Error(xsinit): ngridkscr must be an integer multiple &
	     &of ngridq")')
	write(unitout, '(" ngridkscr : ", 3i6)') input%xs%screening%ngridk
	write(unitout, '(" ngridq    : ", 3i6)') ngridq
	write(unitout, *)
	call terminate
     end if
  else if ((task.ge.440).and.(task.le.459)) then
     ! BSE
     input%groundstate%nosym=input%xs%BSE%nosym
     input%groundstate%reducek=input%xs%BSE%reducek
     input%groundstate%rgkmax=input%xs%BSE%rgkmax
     input%groundstate%vkloff(:)=input%xs%BSE%vkloff(:)
     write(unitout, *)
     write(unitout, '("Info(xsinit): mapping BSE-specific parameters")')
     write(unitout, *)
     if (any(input%groundstate%ngkgrid.ne.ngridq)) then
	write(unitout, *)
	write(unitout, '("Error(xsinit): ngridk must be equal ngridq for the &
	     &BSE-Hamiltonian")')
	write(unitout, '(" ngridk : ", 3i6)') input%groundstate%ngkgrid
	write(unitout, '(" ngridq : ", 3i6)') ngridq
	write(unitout, *)
	call terminate
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
  if (procs.gt.1) then
     call genfilname(basename = 'input%structureoptimization%resume', rank = rank, procs = procs, dotext = '', &
	  filnam = fnresume)
  else
     call genfilname(basename='.input%structureoptimization%resume', dotext='', filnam=fnresume)
  end if
  ! check for stale checkpoint file
  call chkptchk

  ! define checkpoint
  call chkpt(1, (/task/), 'passed xsinit')
end subroutine xsinit
