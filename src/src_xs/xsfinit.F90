
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine xsfinit
  use modmain
  use modxs
  use modmpi
  use m_filedel
  implicit none
  ! local variables
  character(*), parameter :: thisnam='xsfinit'
  character(10) dat, tim
  real(8) :: cput,wallt,cputcum,walltcum
  real(8) :: hrs
  integer :: days,hours,minutes,seconds
  character(256) :: str,str2
  character(256), external :: stringtim,r2str
  ! finalize global counters
  call date_and_time(date=dat,time=tim)
  call cpu_time(cputim0f)
  call system_clock(COUNT=systim0f)
  cput=cputim0f-cputim0i
  wallt=dble(systim0f-systim0i)/dble(cntrate)
  cputcum=cputim0f-cputimcum
  walltcum=dble(systim0f-systimcum)/dble(cntrate)
  ! write out information
  write(unitout,'(a,i8,a)') 'Info('//thisnam//'): task Nr.', task, &
       ' stopped gracefully'
  !call showunits(unitout)
  write(unitout,'(a)') 'Timings: '
  write(unitout,'(a)') '  Date (YYYY-MM-DD) : '//dat(1:4)//'-'//dat(5:6)//'-' &
       //dat(7:8)
  write(unitout,'(a)') '  Time (hh:mm:ss)   : '//tim(1:2)//':'//tim(3:4)//':' &
       //tim(5:6)
  call gentim(cput,hrs,days,hours,minutes,seconds)
  str=stringtim(cput,hrs,days,hours,minutes,seconds)
  write(unitout,'(a,4g18.6)') '  CPU time               : '//trim(str)
  if (procs.eq.1) then
     call gentim(dble(wallt),hrs,days,hours,minutes,seconds)
     str=stringtim(dble(wallt),hrs,days,hours,minutes,seconds)
     str2=r2str(cput/wallt*100,'(f12.2)')
     write(unitout,'(a,4g18.6)') '  wall time              : '//trim(str)
     write(unitout,'(a,g18.6 )') '  CPU load               : '//trim(str2)//' %'
  end if
  call gentim(cputcum,hrs,days,hours,minutes,seconds)
  str=stringtim(cputcum,hrs,days,hours,minutes,seconds)
  write(unitout,'(a,4g18.6)') '  CPU time  (cumulative) : '//trim(str)
  if (procs.eq.1) then
     call gentim(dble(walltcum),hrs,days,hours,minutes,seconds)
     str=stringtim(dble(walltcum),hrs,days,hours,minutes,seconds)
     str2=r2str(cput/wallt*100,'(f12.2)')
     write(unitout,'(a,4g18.6)') '  wall time (cumulative) : '//trim(str)
     write(unitout,'(a,g18.6)')  '  CPU load  (cumulative) : '//trim(str2)//' %'
  end if
  write(unitout,*)
  write(unitout,'("+----------------------------------------------------------&
       &+")')
  write(unitout,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3," (eXcited &
       &States ",I1.1,".",I3.3," ) stopped |")') version,versionxs
  write(unitout,'("+----------------------------------------------------------&
       &+")')
  write(unitout,*)
  close(unitout)
!!$  ! remove tag
!!$  if (.not.tresume) call filedel(trim(fnresume))
!!$  call filedel(trim(fnresume))
end subroutine xsfinit
