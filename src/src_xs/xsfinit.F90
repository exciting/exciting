! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine xsfinit
  use modinput, only: input
  use mod_misc, only: task, versionname
  use modxs, only: cputim0f, systim0f, cputim0i, systim0i,&
                  & cntrate, systimcum, unitout,&
                  & fnresume, xsfileout, fnetim, fnchi0_t,&
                  & fnchi0, fnxtim
  use modmpi
  use m_filedel, only: filedel

  implicit none

  ! Local variables
  character(*), parameter :: thisname = 'xsfinit'
  character(10) :: dat, tim
  real(8) :: cput, wallt, cputcum, walltcum
  real(8) :: hrs
  integer :: days, hours, minutes, seconds
  character(256) :: str1, str2
  character(77) :: string

  ! External functions
  character(256), external :: stringtim, r2str

  ! Some xas specific finalizations
  if(input%xs%bse%xas) call xasfinit

  ! Finalize global counters
  call date_and_time(date=dat, time=tim)
  call cpu_time(cputim0f)
  call system_clock(count=systim0f)
  cput = cputim0f - cputim0i
  wallt = dble(systim0f-systim0i) / dble(cntrate)
  cputcum = cputim0f
  walltcum = dble(systim0f-systimcum) / dble(cntrate)

  ! Write out information
  write(unitout, '(a,i8,a)') 'Info(' // thisname // '): task Nr.', task,&
    & ' stopped gracefully'
  write(unitout,*)
  write(unitout, '(a)') ' Timings: '
  write(unitout, '(a)') '     Date (DD-MM-YYYY)      : '&
    & // dat(7:8) // '-' // dat(5:6) // '-' // dat(1:4)
  write(unitout, '(a)') '     Time (hh:mm:ss)        : '&
    & // tim(1:2) // ':' // tim(3:4) // ':' // tim(5:6)

  call gentim(cput, hrs, days, hours, minutes, seconds)

  str1 = stringtim(cput, hrs, days, hours, minutes, seconds)

  write(unitout, '(a, 4g18.6)')    '     CPU time               : ' // trim(str1)

  call gentim(dble(wallt), hrs, days, hours, minutes, seconds)

  str1 = stringtim(dble(wallt), hrs, days, hours, minutes, seconds)
  str2 = r2str(cput/wallt*100, '(f12.2)')

  write(unitout, '(a, 4g18.6)') '     wall time              : ' // trim(str1)
  write(unitout, '(a,  g18.6)') '     CPU load               : ' // trim(str2) // ' %'

  call gentim(cputcum, hrs, days, hours, minutes, seconds)

  str1 = stringtim(cputcum, hrs, days, hours, minutes, seconds)

  write(unitout, '(a, 4g18.6)') '     CPU time  (cumulative) : ' // trim(str1)

  if(procs .eq. 1) then
    call gentim(dble(walltcum), hrs, days, hours, minutes, seconds)
    str1 = stringtim(dble(walltcum), hrs, days, hours, minutes, seconds)
    str2 = r2str(cput/wallt*100, '(f12.2)')
    write(unitout, '(a, 4g18.6)') '     wall time (cumulative) : ' // trim(str1)
    write(unitout, '(a,  g18.6)') '     CPU load  (cumulative) : ' // trim(str2) // ' %'
  end if

  write(string,'("EXCITING ", a, " stopped for task ",i6)') trim(versionname), task

  call printbox(unitout,"=",string)
  
  write(unitout,*)
  write(unitout,*)
  close(unitout)

  ! Restore global variables
  call restore0
  call restore1
  call restore2

  ! Remove checkpoint file
  call filedel(trim(fnresume))

  if(rank .ne. 0) call filedel(trim(xsfileout))
  if(rank .ne. 0) call filedel(trim(fnetim))
  if(trim(fnchi0_t) .ne. trim(fnchi0)) call filedel(trim(fnchi0_t))
  if(rank .ne. 0) call filedel(trim(fnxtim))

  call barrier(callername=trim(thisname))
      
end subroutine xsfinit
