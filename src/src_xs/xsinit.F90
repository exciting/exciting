
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xsinit(cnt)
#include "../version.inc"
#include "../version2.inc"
#include "../version3.inc"
  use modmain
  use modmpi
  use modxs
  use modfxcifc
  use m_getunit
  use m_genfilname
  implicit none
  ! arguments
  integer, intent(inout) :: cnt
  ! local variables
  character(*), parameter :: thisnam = 'xsinit'
  character(10) dat, tim
  integer :: i
  ! set version of XS part
  call xssetversion
  ! remember how often this routine is called
  cnt=cnt+1
  ! initialize global counters
  call cpu_time(cputim0i)
  call system_clock(COUNT_RATE=cntrate)
  call system_clock(COUNT=systim0i)
  call date_and_time(date=dat,time=tim)
  if (calledxs.eq.1) call system_clock(COUNT=systimcum)
  ! check consistency of rank and procs
  if ((procs.lt.1).or.(procs.gt.maxproc)) then
     write(*,*) 'Error('//trim(thisnam)//'): Error in parallel &
          &initialization: number of processes out of range:',procs
     call terminate
  end if
  if ((rank.gt.procs).or.(rank.lt.0)) then
     write(*,*) 'Error('//trim(thisnam)//'): Error in parallel &
          &initialization: rank out of range:',rank
     call terminate
  end if
  ! set splittfile parameter for splitting of eigenvector files in
  ! parallelization of SCF cycle
  if ((task.ne.301).and.(task.ne.401)) splittfile=.false.
  ! generate resume file
  if (procs.gt.1) then
     call genfilname(basename='resume',rank=rank,procs=procs,dotext='',&
          filnam=fnresume)
  else
     call genfilname(basename='.resume',dotext='',filnam=fnresume)
  end if
  ! name of output file
  call genfilname(nodotpar=.true.,basename='XSINFO',procs=procs,rank=rank, &
       filnam=tdfileout)
  ! reset or append to output file
  call getunit(unitout)
  if (tappinfo.or.(calledxs.gt.1)) then
     open(unitout,file=trim(tdfileout),action='write',position='append')
  else
     open(unitout,file=trim(tdfileout),action='write',status='replace')
  end if
  ! scaling factor for output of energies
  escale=1.d0
  if (tevout) escale=27.2114d0
  ! get exchange-correlation functional data
  call getfxcdata(fxctype,fxcdescr,fxcspin)
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
     if ((procs.gt.1).and.(rank.eq.0)) write(unitout,'(a,2i6)') 'Info('// &
          trim(thisnam)//'):(parallel) master, rank/number of processes:',&
          rank,procs
     if ((procs.gt.1).and.(rank.ne.0)) write(unitout,'(a,2i6)') 'Info('// &
          trim(thisnam)//'):(parallel) slave, rank/number of processes:',&
          rank,procs
     if (notelns.gt.0) then
        write(unitout,*)
        write(unitout,'("Notes :")')
        do i=1,notelns
           write(unitout,'(A)') notes(i)
        end do
     end if
     write(unitout,*)
  end if
  write(unitout,'(a)') 'Date (YYYY-MM-DD) : '//dat(1:4)//'-'//dat(5:6)//'-'// &
       dat(7:8)
  write(unitout,'(a)') 'Time (hh:mm:ss)   : '//tim(1:2)//':'//tim(3:4)//':'// &
       tim(5:6)
  write(unitout,*)
  write(unitout,'(a,i6,a)') 'Info('//trim(thisnam)//'): task Nr.', &
       task,' started'
  call flushifc(unitout)
end subroutine xsinit
