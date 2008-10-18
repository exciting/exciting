
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xsinit
#include "../version.inc"
  use modmain
  use modmpi
  use modxs
  use modfxcifc
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(10) dat, tim
  integer :: i

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
  call date_and_time(date=dat,time=tim)
  if (calledxs.eq.1) call system_clock(COUNT=systimcum)  

  !-----------------------------------!
  !     parallelization variables     !
  !-----------------------------------!
  if ((procs.lt.1).or.(procs.gt.maxproc)) then
     write(unitout,*)
     write(unitout,'("Error(xsinit): Error in parallel initialization: number &
     &of processes out of range: ",i6)') procs
     write(unitout,*)
     call terminate
  end if
  if ((rank.gt.procs).or.(rank.lt.0)) then
     write(unitout,*)
     write(unitout,'("Error(xsinit): Error in parallel initialization: rank &
     &out of range: ",i6)') rank
     write(unitout,*)
     call terminate
  end if

  !------------------------!
  !     spin variables     !
  !------------------------!
  ! warn for spin polarized calculations
  if (spinpol) then
     write(unitout,*)
     write(unitout,'("Warning(xsinit): calculation is spin-polarized - &
     &formalism may be incomplete")')
     write(unitout,*)
  end if
  ! no spin-spirals
  if (spinsprl) then
     write(unitout,*)
     write(unitout,'("Error(xsinit): xs-part not working for spin-spirals")')
     write(unitout,*)
     call terminate
  end if

  !------------------------------------!
  !     angular momentum variables     !
  !------------------------------------!
  if (lmaxapwwf.eq.-1) lmaxapwwf=lmaxmat
  lmmaxapwwf=(lmaxapwwf+1)**2
  lmmaxemat=(lmaxemat+1)**2
  if (lmaxapwwf.gt.lmaxapw) then
     write(unitout,*)
     write(unitout,'("Error(xsinit): lmaxapwwf > lmaxapw: ",i6)') lmaxapwwf
     write(unitout,*)
     call terminate
  end if
  if (lmaxemat.gt.lmaxapw) then
     write(unitout,*)
     write(unitout,'("Error(xsinit): lmaxemat > lmaxapw: ",i6)') lmaxemat
     write(unitout,*)
     call terminate
  end if
  if (lmaxemat.gt.lmaxapwwf) then
     write(unitout,*)
     write(unitout,'("Warning(xsinit): lmaxemat > lmaxapwwf: ",i6)') lmaxemat
     write(unitout,*)
  end if

  !---------------------!
  !     k-point set     !
  !---------------------!
  if (any(ngridkscr.eq.0)) ngridkscr(:)=ngridk(:)
  if (any(vkloffscr.eq.-1.d0)) vkloffscr(:)=vkloff(:)
  if (any(vkloffbse.eq.-1.d0)) vkloffbse(:)=vkloff(:)

  !---------------------!
  !     G+k vectors     !
  !---------------------!
  if (rgkmaxscr.eq.0.d0) rgkmaxscr=rgkmax
  if (rgkmaxbse.eq.0.d0) rgkmaxbse=rgkmax

  !------------------------------------!
  !     secular equation variables     !
  !------------------------------------!
  if (nemptyscr.eq.0) nemptyscr=nempty
  ! set splittfile parameter for splitting of eigenvector files in
  ! parallelization of SCF cycle
  if ((task.ne.301).and.(task.ne.401)) splittfile=.false.
  
  !----------------------------!
  !     response functions     !
  !----------------------------!
  tscreen=.false.
  if ((task.ge.400).and.(task.le.499)) tscreen=.true.
  ! type of response functions
  if (rsptype.ne.'reta') then
     write(unitout,*)
     write(unitout,'("Error(xsinit): only retarded response functions &
     &implemented")')
     write(unitout,*)
     call terminate
  end if
  ! tetrahedron method not implemented for analytic continuation
  if (tetradf.and.acont) then
     write(unitout,*)
     write(unitout,'("Error(xsinit): tetrahedron method does not work in &
     	& combination with analytic continuation")')
     write(unitout,*)
     call terminate
  end if
  ! if imaginary frequencies intervals are not specified
  if (nwacont.eq.0) nwacont=nwdos
  nwdf=nwdos
  if (acont) nwdf=nwacont
  ! get exchange-correlation kernel functional data
  call getfxcdata(fxctype,fxcdescr,fxcspin)
  
  !-----------------------!
  !     miscellaneous     !
  !-----------------------!
  ! scaling factor for output of energies
  escale=1.d0
  if (tevout) escale=27.2114d0

  !----------------------------------!
  !     task dependent variables     !
  !----------------------------------!
  if ((task.ge.401).and.(task.le.430)) then
     nosym=nosymscr
     reducek=reducekscr
     rgkmax=rgkmaxscr
     nempty=nemptyscr
     ngridk(:)=ngridkscr(:)
     vkloff(:)=vkloffscr(:)
     write(unitout,*)
     write(unitout,'("Info(xsinit): mapping screening-specific parameters")')
     write(unitout,*)
  else if ((task.ge.440).and.(task.le.445)) then
     nosym=nosymbse
     reducek=reducekbse
     rgkmax=rgkmaxbse
     vkloff(:)=vkloffbse(:)
     write(unitout,*)
     write(unitout,'("Info(xsinit): mapping BSE-specific parameters")')
     write(unitout,*)
  end if

  !---------------------!
  !     checkpoints     !
  !---------------------!
  if (procs.gt.1) then
     call genfilname(basename='resume',rank=rank,procs=procs,dotext='',&
          filnam=fnresume)
  else
     call genfilname(basename='.resume',dotext='',filnam=fnresume)
  end if
  ! check for stale checkpoint file
  call chkptchk
  
  !---------------------!
  !     output file     !
  !---------------------!
  ! set version of XS part
  call xssetversion
  ! name of output file
  call genfilname(nodotpar=.true.,basename='XSINFO',procs=procs,rank=rank, &
       filnam=xsfileout)
  ! reset or append to output file
  call getunit(unitout)
  if (tappinfo.or.(calledxs.gt.1)) then
     open(unitout,file=trim(xsfileout),action='write',position='append')
  else
     open(unitout,file=trim(xsfileout),action='write',status='replace')
  end if
  ! write to info file
  if (calledxs.eq.1) then
     write(unitout,*)
     write(unitout,'("+-------------------------------------------------------&
          &---+")')
     write(unitout,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3," (eXcited &
          &States ",I1.1,".",I3.3,") started  |")') version,versionxs
     write(unitout,'("| git hash id : ",2a20,"   |")') GITHASH,GITHASH2
#ifdef LOCALCHG
     write(unitout,'("| Warning     : source codes deviates from the git hash &
          &id |")')
#endif
     write(unitout,'("+ ------------------------------------------------------&
          &---+")')
#ifdef MPI
     write(unitout,'("compiled for MPI execution")') 
#endif
#ifndef MPI
     write(unitout,'("compiled for serial execution")') 
#endif
     if ((procs.gt.1).and.(rank.eq.0)) write(unitout,'("(parallel) master, &
     	&rank/number of processes:")') rank,procs
     if ((procs.gt.1).and.(rank.eq.0)) write(unitout,'("(parallel) slave,  &
     	&rank/number of processes:")') rank,procs
     if (notelns.gt.0) then
        write(unitout,*)
        write(unitout,'("Notes :")')
        do i=1,notelns
           write(unitout,'(A)') notes(i)
        end do
     end if
     write(unitout,*)
  end if
  write(unitout,*)
  write(unitout,'("Date (YYYY-MM-DD) : ",A4,"-",A2,"-",A2)') dat(1:4), &
  	dat(5:6),dat(7:8)
  write(unitout,'("Time (hh:mm:ss)   : ",A2,":",A2,":",A2)') tim(1:2), &
  	tim(3:4),tim(5:6)
  write(unitout,'("Info(xsinit): task Nr.",i6," started")') task
  call flushifc(unitout)
  
  ! define checkpoint
  call chkpt(1,(/task/),'passed xsinit')
end subroutine xsinit
