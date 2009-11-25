!
!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module modtimer
      Implicit None
      Type timer
         Private
     ! name of timer
         Character (256) :: name
     ! type of timer (system clock or CPU clock)
         Logical :: systemclock
     ! true if the timer is running
         Logical :: running
     ! mark the current time
         Real (8) :: mark
     ! the elapsed time
         Real (8) :: elapsed
     ! number of measured intervals
         Integer :: intervals
     ! last interval
         Real (8) :: int_last
     ! the smallest measured interval
         Real (8) :: int_small
     ! the largest measured interval
         Real (8) :: int_large
     ! the mean value of the measured intervals
         Real (8) :: int_mean
      End Type timer
      Private :: get_system_time
      Private :: get_cpu_time
      Private :: get_time
      Private :: new_interval
Contains
!
!
      Subroutine new_timer (t, typ, name)
         Implicit None
         Type (timer) :: t
         Character (*), Intent (In) :: typ
         Character (*), Intent (In), Optional :: name
         t%name = trim (adjustl(typ))
         If (present(name)) t%name = trim (name)
         t%running = .False.
         t%elapsed = 0.d0
         t%intervals = 0
         t%int_last = 0.d0
         t%int_small = 1.d10
         t%int_large = - 1.d10
         t%int_mean = 0.d0
         Select Case (trim(adjustl(typ)))
         Case ('system', 'wall')
            t%systemclock = .True.
         Case ('cpu')
            t%systemclock = .False.
         Case Default
            Write (*,*)
            Write (*, '("Error(timing:new_timer): unknown clocktype: ",&
           & a)') trim (typ)
            Write (*,*)
            Stop
         End Select
      End Subroutine new_timer
!
!
      Subroutine tic_timer (t)
         Implicit None
         Type (timer) :: t
         If (t%running) Then
            Call new_interval (t, get_time(t))
         Else
            t%running = .True.
            t%mark = get_time (t)
         End If
      End Subroutine tic_timer
!
!
      Subroutine toc_timer (t)
         Implicit None
         Type (timer) :: t
         If ( .Not. t%running) Then
            Write (*, '("Error(timing:toc_timer): timer already stopped&
           &")')
            Stop
         End If
         Call new_interval (t, get_time(t))
         t%running = .False.
      End Subroutine toc_timer
!
!
      Subroutine report_timer (t, un)
         Use modinput
         Implicit None
!
         Type (timer) :: t
         Integer, Intent (In) :: un
         Write (un, '("Timer: ", a)') trim (t%name)
         If (t%running) Then
            Write (un, '(" timer is running")')
         Else
            Write (un, '(" timer is stopped")')
         End If
         If (t%systemclock) Then
            Write (un, '(" timer is using system-clock")')
         Else
            Write (un, '(" timer is using CPU-clock")')
         End If
         Write (un, '(" elapsed time	    : ", g18.10)') t%elapsed
         Write (un, '(" intervals	    : ", i8)') t%intervals
         Write (un, '(" last interval     : ", g18.10)') t%int_last
         Write (un, '(" mean interval     : ", g18.10)') t%int_mean
         Write (un, '(" smallest interval : ", g18.10)') t%int_small
         Write (un, '(" largest  interval : ", g18.10)') t%int_large
         Write (un, '(" current mark	    : ", g18.10)') t%mark
      End Subroutine report_timer
!
      Character (256) Function get_name_timer (t)
         Implicit None
         Type (timer) :: t
         get_name_timer = trim (t%name)
      End Function get_name_timer
!
      Logical Function is_systemclock_timer (t)
         Implicit None
         Type (timer) :: t
         is_systemclock_timer = t%systemclock
      End Function is_systemclock_timer
!
      Logical Function is_cpuclock_timer (t)
         Implicit None
         Type (timer) :: t
         is_cpuclock_timer = .Not. is_systemclock_timer (t)
      End Function is_cpuclock_timer
!
      Logical Function is_running_timer (t)
         Implicit None
         Type (timer) :: t
         is_running_timer = t%running
      End Function is_running_timer
!
      Real (8) Function get_elapsed_timer (t)
         Implicit None
         Type (timer) :: t
         get_elapsed_timer = t%elapsed
      End Function get_elapsed_timer
!
      Real (8) Function get_mark_timer (t)
         Implicit None
         Type (timer) :: t
         get_mark_timer = t%mark
      End Function get_mark_timer
!
      Integer Function get_intervals_timer (t)
         Implicit None
         Type (timer) :: t
         get_intervals_timer = t%intervals
      End Function get_intervals_timer
!
      Real (8) Function get_intlast_timer (t)
         Implicit None
         Type (timer) :: t
         get_intlast_timer = t%int_last
      End Function get_intlast_timer
!
      Real (8) Function get_intsmall_timer (t)
         Implicit None
         Type (timer) :: t
         get_intsmall_timer = t%int_small
      End Function get_intsmall_timer
!
      Real (8) Function get_intlarge_timer (t)
         Implicit None
         Type (timer) :: t
         get_intlarge_timer = t%int_large
      End Function get_intlarge_timer
!
      Real (8) Function get_intmean_timer (t)
         Implicit None
         Type (timer) :: t
         get_intmean_timer = t%int_mean
      End Function get_intmean_timer
!
  !
  !private routines
  !
!
!
      Subroutine new_interval (t, s)
         Implicit None
         Type (timer) :: t
         Real (8), Intent (In) :: s
         Real (8) :: sint
         sint = s - t%mark
         t%elapsed = t%elapsed + sint
         t%mark = s
         t%intervals = t%intervals + 1
         t%int_last = sint
         t%int_small = Min (t%int_small, sint)
         t%int_large = Max (t%int_large, sint)
         t%int_mean = (dble(t%intervals-1)*t%int_mean+sint) / dble &
        & (t%intervals)
      End Subroutine new_interval
!
      Real (8) Function get_time (t)
         Implicit None
         Type (timer) :: t
         If (t%systemclock) Then
            get_time = get_system_time ()
         Else
            get_time = get_cpu_time ()
         End If
      End Function get_time
!
      Real (8) Function get_system_time ()
         Implicit None
         Integer :: cnt, cntr
         Call system_clock (count=cnt, count_rate=cntr)
         get_system_time = dble (cnt) / dble (cntr)
      End Function get_system_time
!
      Real (8) Function get_cpu_time ()
         Implicit None
         Call cpu_time (get_cpu_time)
      End Function get_cpu_time
!
End Module modtimer
!
!
!///////////////////////////////////////////////////////////////////////////
!
!
Module modtimer2
      Use modtimer
      Type timer2
         Private
         Type (timer) :: sys
         Type (timer) :: cpu
         Logical :: tsys
         Logical :: tcpu
         Character (256) :: name
         Logical :: tsysdefault
      End Type timer2
!
Contains
!
!
      Subroutine new_timer2 (t, typ, name)
         Implicit None
         Type (timer2) :: t
         Character (*), Intent (In), Optional :: typ
         Character (*), Intent (In), Optional :: name
         t%tsysdefault = .True.
         t%tsys = .True.
         t%tcpu = .True.
         If (present(typ)) Then
            Select Case (trim(adjustl(typ)))
            Case ('combined')
            Case ('system', 'wall')
               t%tcpu = .False.
            Case ('cpu')
               t%tsys = .False.
            Case Default
               Write (*,*)
               Write (*, '("Error(timing:new_timer2): unknown type: ", &
              &a)') trim (typ)
               Write (*,*)
               Stop
            End Select
         End If
         t%name = 'unnamed'
         If (present(name)) t%name = trim (name)
         If (t%tsys) Call new_timer (t%sys, 'system', 'combined:System')
         If (t%tcpu) Call new_timer (t%cpu, 'cpu', 'combined:CPU')
      End Subroutine new_timer2
!
!
      Subroutine tic_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Call tic_timer (t%sys)
         If (t%tcpu) Call tic_timer (t%cpu)
      End Subroutine tic_timer2
!
!
      Subroutine toc_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Call toc_timer (t%sys)
         If (t%tcpu) Call toc_timer (t%cpu)
      End Subroutine toc_timer2
!
!
      Subroutine report_timer2 (t, un, string)
         Implicit None
         Type (timer2) :: t
         Integer, Intent (In) :: un
         Character (*), Intent (In), Optional :: string
         Real (8) :: r
         Write (un,*)
         Write (un, '("Combined timer; summary below")')
         If (present(string)) write (un, '(" ID		 : ", a)') trim &
        & (string)
         If (t%tsys) Call report_timer (t%sys, un)
         If (t%tcpu) Call report_timer (t%cpu, un)
         If (t%tsys .And. t%tcpu) Then
            r = get_elapsed_timer (t%sys)
            If (r .Ne. 0.d0) r = get_elapsed_timer (t%cpu) / r
            Write (un, '(" CPU/System (%)    : ", g18.10)') r * 100.d0
         End If
         Write (un,*)
      End Subroutine report_timer2
!
      Character (256) Function get_name_timer2 (t)
         Implicit None
         Type (timer2) :: t
         get_name_timer2 = trim (t%name)
      End Function get_name_timer2
!
      Logical Function is_systemclock_timer2 (t)
         Implicit None
         Type (timer2) :: t
         is_systemclock_timer2 = t%tsys
      End Function is_systemclock_timer2
!
      Logical Function is_cpuclock_timer2 (t)
         Implicit None
         Type (timer2) :: t
         is_cpuclock_timer2 = t%tcpu
      End Function is_cpuclock_timer2
!
      Logical Function is_running_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) is_running_timer2 = is_running_timer (t%sys)
         If (t%tcpu) is_running_timer2 = is_running_timer (t%cpu)
      End Function is_running_timer2
!
      Real (8) Function get_elapsed_system_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Then
            get_elapsed_system_timer2 = get_elapsed_timer (t%sys)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_elapsed_system_timer2): no Sy&
           &stem clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_elapsed_system_timer2
!
      Real (8) Function get_elapsed_cpu_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tcpu) Then
            get_elapsed_cpu_timer2 = get_elapsed_timer (t%cpu)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_elapsed_cpu_timer2): no CPU c&
           &lock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_elapsed_cpu_timer2
!
      Real (8) Function get_elapsed_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsysdefault .And. t%tsys) Then
            get_elapsed_timer2 = get_elapsed_system_timer2 (t)
         Else
            get_elapsed_timer2 = get_elapsed_cpu_timer2 (t)
         End If
      End Function get_elapsed_timer2
!
      Real (8) Function get_mark_system_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Then
            get_mark_system_timer2 = get_mark_timer (t%sys)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_mark_system_timer2): no Syste&
           &m clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_mark_system_timer2
!
      Real (8) Function get_mark_cpu_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tcpu) Then
            get_mark_cpu_timer2 = get_mark_timer (t%cpu)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_mark_cpu_timer2): no CPU cloc&
           &k timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_mark_cpu_timer2
!
      Real (8) Function get_mark_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsysdefault .And. t%tsys) Then
            get_mark_timer2 = get_mark_system_timer2 (t)
         Else
            get_mark_timer2 = get_mark_cpu_timer2 (t)
         End If
      End Function get_mark_timer2
!
      Integer Function get_intervals_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) get_intervals_timer2 = get_intervals_timer (t%sys)
         If (t%tcpu) get_intervals_timer2 = get_intervals_timer (t%cpu)
      End Function get_intervals_timer2
!
      Real (8) Function get_intlast_system_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Then
            get_intlast_system_timer2 = get_intlast_timer (t%sys)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intlast_system_timer2): no Sy&
           &stem clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intlast_system_timer2
!
      Real (8) Function get_intlast_cpu_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tcpu) Then
            get_intlast_cpu_timer2 = get_intlast_timer (t%cpu)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intlast_cpu_timer2): no CPU c&
           &lock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intlast_cpu_timer2
!
      Real (8) Function get_intlast_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsysdefault .And. t%tsys) Then
            get_intlast_timer2 = get_intlast_system_timer2 (t)
         Else
            get_intlast_timer2 = get_intlast_cpu_timer2 (t)
         End If
      End Function get_intlast_timer2
!
      Real (8) Function get_intsmall_system_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Then
            get_intsmall_system_timer2 = get_intsmall_timer (t%sys)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intsmall_system_timer2): no S&
           &ystem clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intsmall_system_timer2
!
      Real (8) Function get_intsmall_cpu_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tcpu) Then
            get_intsmall_cpu_timer2 = get_intsmall_timer (t%cpu)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intsmall_cpu_timer2): no CPU &
           &clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intsmall_cpu_timer2
!
      Real (8) Function get_intsmall_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsysdefault .And. t%tsys) Then
            get_intsmall_timer2 = get_intsmall_system_timer2 (t)
         Else
            get_intsmall_timer2 = get_intsmall_cpu_timer2 (t)
         End If
      End Function get_intsmall_timer2
!
      Real (8) Function get_intlarge_system_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Then
            get_intlarge_system_timer2 = get_intlarge_timer (t%sys)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intlarge_system_timer2): no S&
           &ystem clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intlarge_system_timer2
!
      Real (8) Function get_intlarge_cpu_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tcpu) Then
            get_intlarge_cpu_timer2 = get_intlarge_timer (t%cpu)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intlarge_cpu_timer2): no CPU &
           &clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intlarge_cpu_timer2
!
      Real (8) Function get_intlarge_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsysdefault .And. t%tsys) Then
            get_intlarge_timer2 = get_intlarge_system_timer2 (t)
         Else
            get_intlarge_timer2 = get_intlarge_cpu_timer2 (t)
         End If
      End Function get_intlarge_timer2
!
      Real (8) Function get_intmean_system_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsys) Then
            get_intmean_system_timer2 = get_intmean_timer (t%sys)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intmean_system_timer2): no Sy&
           &stem clock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intmean_system_timer2
!
      Real (8) Function get_intmean_cpu_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tcpu) Then
            get_intmean_cpu_timer2 = get_intmean_timer (t%cpu)
         Else
            Write (*,*)
            Write (*, '("Error(timing:get_intmean_cpu_timer2): no CPU c&
           &lock timer running")')
            Write (*,*)
            Stop
         End If
      End Function get_intmean_cpu_timer2
!
      Real (8) Function get_intmean_timer2 (t)
         Implicit None
         Type (timer2) :: t
         If (t%tsysdefault .And. t%tsys) Then
            get_intmean_timer2 = get_intmean_system_timer2 (t)
         Else
            get_intmean_timer2 = get_intmean_cpu_timer2 (t)
         End If
      End Function get_intmean_timer2
!
End Module modtimer2
