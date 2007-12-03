
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tdinit
  use modmain
  use modxs
  use modfxcifc
  use modmpi
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tdinit'
  character(10) dat, tim
  integer :: un,i
  logical :: ex

  ! assign xs code version
  versionxs=(/0,85/)

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

  ! generate resume file
  if (procs.gt.1) then
     call genfilname(basename='resume',rank=rank,procs=procs,dotext='',&
          filnam=fnresume)
  else
     call genfilname(basename='.resume',dotext='',filnam=fnresume)
  end if

!!$  call getunit(un)
!!$  ! initialize for first call to main routine
!!$  if (calledtd.eq.1) resumechkpts=0
!!$  ! read in checkpoint if present
!!$  call resread(un,resumetask,resumechkpts,tresume)
!!$  ! checkpointing starts
!!$  if (.not.tresume) then
!!$     resumechkpts(:,1)=0
!!$     call resupd(un,task,resumechkpts,' : prolog')
!!$  end if

  !initialize global counters
  call cpu_time(cputim0i)
  call system_clock(COUNT_RATE=cntrate)
  call system_clock(COUNT=systim0i)
  call date_and_time(date=dat,time=tim)
  if (calledtd.eq.1) call system_clock(COUNT=systimcum)

  ! name of output file
  call genfilname(nodotpar=.true.,basename='XSINFO',&
       procs=procs,rank=rank,filnam=tdfileout)

  ! reset or append to output file
  call getunit(unitout)
  if ( tappinfo.or.(calledtd.gt.1)) then
     open(unitout,file=trim(tdfileout),action='write',position='append')
  else
     open(unitout,file=trim(tdfileout),action='write',status='replace')
  end if

  ! scaling factor for output of energies
  escale=1.d0
  if (tevout) escale=27.2117

  ! get exchange-correlation functional data
  call getfxcdata(fxctype,fxcdescr,fxcspin)

  ! write to info file
  if (calledtd == 1) then
     write(unitout,*)
     write(unitout,'("+-------------------------------------------------------&
          &---+")')
     write(unitout,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3,&
          " (eXcited States "I1.1,".",I3.3," ) started |")') version,versionxs
     write(unitout,'("+-------------------------------------------------------&
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

end subroutine tdinit
