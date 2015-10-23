!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine xsfinit
      Use modmain
      Use modxs
      Use modmpi
      Use m_filedel
      Implicit None

  ! local variables
      Character (*), Parameter :: thisnam = 'xsfinit'
      Character (10) :: dat, tim
      Real (8) :: cput, wallt, cputcum, walltcum
      Real (8) :: hrs
      Integer :: days, hours, minutes, seconds
      Character (256) :: str1, str2
      Character (256), External :: stringtim, r2str
      character*(77) :: string

  ! finalize global counters
      Call date_and_time (date=dat, time=tim)
      Call cpu_time (cputim0f)
      Call system_clock (COUNT=systim0f)
      cput = cputim0f - cputim0i
      wallt = dble (systim0f-systim0i) / dble (cntrate)
      cputcum = cputim0f - cputimcum
      walltcum = dble (systim0f-systimcum) / dble (cntrate)

  ! write out information
      Write (unitout, '(a,i8,a)') 'Info(' // thisnam // '): task Nr.', task, ' stopped gracefully'
      write (unitout,*)
      Write (unitout, '(a)') ' Timings: '
      Write (unitout, '(a)') '     Date (DD-MM-YYYY)      : ' // dat (7:8) // '-' // dat (5:6) // '-' // dat (1:4)
      Write (unitout, '(a)') '     Time (hh:mm:ss)        : ' // tim (1:2) // ':' // tim (3:4) // ':' // tim (5:6)
      Call gentim (cput, hrs, days, hours, minutes, seconds)
      str1 = stringtim (cput, hrs, days, hours, minutes, seconds)
      Write (unitout, '(a, 4g18.6)')    '     CPU time               : ' // trim (str1)
!      If (procs .Eq. 1) Then
         Call gentim (dble(wallt), hrs, days, hours, minutes, seconds)
         str1 = stringtim (dble(wallt), hrs, days, hours, minutes, seconds)
         str2 = r2str (cput/wallt*100, '(f12.2)')
         Write (unitout, '(a, 4g18.6)') '     wall time              : ' // trim (str1)
         Write (unitout, '(a,  g18.6)') '     CPU load               : ' // trim (str2) // ' %'
!      End If
      Call gentim (cputcum, hrs, days, hours, minutes, seconds)
      str1 = stringtim (cputcum, hrs, days, hours, minutes, seconds)
      Write (unitout, '(a, 4g18.6)') '     CPU time  (cumulative) : ' // trim (str1)
      If (procs .Eq. 1) Then
         Call gentim (dble(walltcum), hrs, days, hours, minutes, seconds)
         str1 = stringtim (dble(walltcum), hrs, days, hours, minutes, seconds)
         str2 = r2str (cput/wallt*100, '(f12.2)')
         Write (unitout, '(a, 4g18.6)') '     wall time (cumulative) : ' // trim (str1)
         Write (unitout, '(a,  g18.6)') '     CPU load  (cumulative) : ' // trim (str2) // ' %'
      End If

      write(string,'("EXCITING ", a, " stopped for task ",i6)') trim(versionname), task
      call printbox(unitout,"=",string)
      write(unitout,*)
      write(unitout,*)
      close(unitout)

  ! restore global variables
      Call restore0
      Call restore1
      Call restore2

  ! remove checkpoint file
      Call filedel (trim(fnresume))
      If (rank .Ne. 0) Call filedel (trim(xsfileout))
      If (rank .Ne. 0) Call filedel (trim(fnetim))
      If (trim(fnchi0_t) .Ne. trim(fnchi0)) Call filedel (trim(fnchi0_t))
      If (rank .Ne. 0) Call filedel (trim(fnxtim))
      
End Subroutine xsfinit
