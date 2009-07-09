



! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtimer
  implicit none
  type timer
     private
     ! name of timer
     character(256) :: name
     ! type of timer (system clock or CPU clock)
     logical :: systemclock
     ! true if the timer is running
     logical :: running
     ! mark the current time
     real(8) :: mark
     ! the elapsed time
     real(8) :: elapsed
     ! number of measured intervals
     integer :: intervals
     ! last interval
     real(8) :: int_last
     ! the smallest measured interval
     real(8) :: int_small
     ! the largest measured interval
     real(8) :: int_large
     ! the mean value of the measured intervals
     real(8) :: int_mean
  end type timer
  private :: get_system_time
  private :: get_cpu_time
  private :: get_time
  private :: new_interval
contains


subroutine new_timer(t, typ, name)
    implicit none
    type(timer) :: t
    character(*), intent(in) :: typ
    character(*), intent(in), optional :: name
    t%name=trim(adjustl(typ))
    if (present(name)) t%name=trim(name)
    t%running=.false.
    t%elapsed=0.d0
    t%intervals=0
    t%int_last=0.d0
    t%int_small=1.d10
    t%int_large=-1.d10
    t%int_mean=0.d0
    select case (trim(adjustl(typ)))
    case('system', 'wall')
       t%systemclock=.true.
    case('cpu')
       t%systemclock=.false.
    case default
       write(*, *)
       write(*, '("Error(timing:new_timer): unknown clocktype: ", a)') &
	    trim(typ)
       write(*, *)
       stop
    end select
  end subroutine new_timer


subroutine tic_timer(t)
    implicit none
    type(timer) :: t
    if (t%running) then
       call new_interval(t, get_time(t))
    else
       t%running=.true.
       t%mark=get_time(t)
    end if
  end subroutine tic_timer


subroutine toc_timer(t)
    implicit none
    type(timer) :: t
    if (.not.t%running) then
       write(*, '("Error(timing:toc_timer): timer already stopped")')
       stop
    end if
    call new_interval(t, get_time(t))
    t%running=.false.
  end subroutine toc_timer


subroutine report_timer(t, un)
use modinput
    implicit none

    type(timer) :: t
    integer, intent(in) :: un
    write(un, '("Timer: ", a)') trim(t%name)
    if (t%running) then
       write(un, '(" timer is running")')
    else
       write(un, '(" timer is stopped")')
    end if
    if (t%systemclock) then
       write(un, '(" timer is using system-clock")')
    else
       write(un, '(" timer is using CPU-clock")')
    end if
    write(un, '(" elapsed time	    : ", g18.10)') t%elapsed
    write(un, '(" intervals	    : ", i8)') t%intervals
    write(un, '(" last interval     : ", g18.10)') t%int_last
    write(un, '(" mean interval     : ", g18.10)') t%int_mean
    write(un, '(" smallest interval : ", g18.10)') t%int_small
    write(un, '(" largest  interval : ", g18.10)') t%int_large
    write(un, '(" current mark	    : ", g18.10)') t%mark
  end subroutine report_timer

  character(256) function get_name_timer(t)
    implicit none
    type(timer) :: t
    get_name_timer=trim(t%name)
  end function get_name_timer

  logical function is_systemclock_timer(t)
    implicit none
    type(timer) :: t
    is_systemclock_timer=t%systemclock
  end function is_systemclock_timer

  logical function is_cpuclock_timer(t)
    implicit none
    type(timer) :: t
    is_cpuclock_timer=.not.is_systemclock_timer(t)
  end function is_cpuclock_timer

  logical function is_running_timer(t)
    implicit none
    type(timer) :: t
    is_running_timer=t%running
  end function is_running_timer

  real(8) function get_elapsed_timer(t)
    implicit none
    type(timer) :: t
    get_elapsed_timer=t%elapsed
  end function get_elapsed_timer

  real(8) function get_mark_timer(t)
    implicit none
    type(timer) :: t
    get_mark_timer=t%mark
  end function get_mark_timer

  integer function get_intervals_timer(t)
    implicit none
    type(timer) :: t
    get_intervals_timer=t%intervals
  end function get_intervals_timer

  real(8) function get_intlast_timer(t)
    implicit none
    type(timer) :: t
    get_intlast_timer=t%int_last
  end function get_intlast_timer

  real(8) function get_intsmall_timer(t)
    implicit none
    type(timer) :: t
    get_intsmall_timer=t%int_small
  end function get_intsmall_timer

  real(8) function get_intlarge_timer(t)
    implicit none
    type(timer) :: t
    get_intlarge_timer=t%int_large
  end function get_intlarge_timer

  real(8) function get_intmean_timer(t)
    implicit none
    type(timer) :: t
    get_intmean_timer=t%int_mean
  end function get_intmean_timer

  !
  !private routines
  !


subroutine new_interval(t, s)
    implicit none
    type(timer) :: t
    real(8), intent(in) :: s
    real(8) :: sint
    sint=s-t%mark
    t%elapsed=t%elapsed+sint
    t%mark=s
    t%intervals=t%intervals+1
    t%int_last=sint
    t%int_small=min(t%int_small, sint)
    t%int_large=max(t%int_large, sint)
    t%int_mean=(dble(t%intervals-1)*t%int_mean+sint)/dble(t%intervals)
  end subroutine new_interval

  real(8) function get_time(t)
    implicit none
    type(timer) :: t
    if (t%systemclock) then
       get_time=get_system_time()
    else
       get_time=get_cpu_time()
    end if
  end function get_time

  real(8) function get_system_time()
    implicit none
    integer :: cnt, cntr
    call system_clock(count=cnt, count_rate=cntr)
    get_system_time=dble(cnt)/dble(cntr)
  end function get_system_time

  real(8) function get_cpu_time()
    implicit none
    call cpu_time(get_cpu_time)
  end function get_cpu_time

end module modtimer


!///////////////////////////////////////////////////////////////////////////


module modtimer2
  use modtimer
  type timer2
     private
     type(timer) :: sys
     type(timer) :: cpu
     logical :: tsys
     logical :: tcpu
     character(256) :: name
     logical :: tsysdefault
  end type timer2

contains


subroutine new_timer2(t, typ, name)
    implicit none
    type(timer2) :: t
    character(*), intent(in), optional :: typ
    character(*), intent(in), optional :: name
    t%tsysdefault=.true.
    t%tsys=.true.
    t%tcpu=.true.
    if (present(typ)) then
    select case (trim(adjustl(typ)))
       case('combined')
       case('system', 'wall')
	  t%tcpu=.false.
       case('cpu')
	  t%tsys=.false.
       case default
	  write(*, *)
	  write(*, '("Error(timing:new_timer2): unknown type: ", a)') trim(typ)
	  write(*, *)
	  stop
    end select
    end if
    t%name='unnamed'
    if (present(name)) t%name=trim(name)
    if (t%tsys) call new_timer(t%sys, 'system', 'combined:System')
    if (t%tcpu) call new_timer(t%cpu, 'cpu', 'combined:CPU')
  end subroutine new_timer2


subroutine tic_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) call tic_timer(t%sys)
    if (t%tcpu) call tic_timer(t%cpu)
  end subroutine tic_timer2


subroutine toc_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) call toc_timer(t%sys)
    if (t%tcpu) call toc_timer(t%cpu)
  end subroutine toc_timer2


subroutine report_timer2(t, un, string)
    implicit none
    type(timer2) :: t
    integer, intent(in) :: un
    character(*), intent(in), optional :: string
    real(8) :: r
    write(un, *)
    write(un, '("Combined timer; summary below")')
    if (present(string)) write(un, '(" ID		 : ", a)') trim(string)
    if (t%tsys) call report_timer(t%sys, un)
    if (t%tcpu) call report_timer(t%cpu, un)
    if (t%tsys.and.t%tcpu) then
      r=get_elapsed_timer(t%sys)
      if (r.ne.0.d0) r=get_elapsed_timer(t%cpu)/r
      write(un, '(" CPU/System (%)    : ", g18.10)') r*100.d0
    end if
    write(un, *)
  end subroutine report_timer2

  character(256) function get_name_timer2(t)
    implicit none
    type(timer2) :: t
    get_name_timer2=trim(t%name)
  end function get_name_timer2

  logical function is_systemclock_timer2(t)
    implicit none
    type(timer2) :: t
    is_systemclock_timer2=t%tsys
  end function is_systemclock_timer2

  logical function is_cpuclock_timer2(t)
    implicit none
    type(timer2) :: t
    is_cpuclock_timer2=t%tcpu
  end function is_cpuclock_timer2

  logical function is_running_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) is_running_timer2=is_running_timer(t%sys)
    if (t%tcpu) is_running_timer2=is_running_timer(t%cpu)
  end function is_running_timer2

  real(8) function get_elapsed_system_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) then
       get_elapsed_system_timer2=get_elapsed_timer(t%sys)
    else
       write(*, *)
       write(*, '("Error(timing:get_elapsed_system_timer2): no System clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_elapsed_system_timer2

  real(8) function get_elapsed_cpu_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tcpu) then
       get_elapsed_cpu_timer2=get_elapsed_timer(t%cpu)
    else
       write(*, *)
       write(*, '("Error(timing:get_elapsed_cpu_timer2): no CPU clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_elapsed_cpu_timer2

  real(8) function get_elapsed_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsysdefault.and.t%tsys) then
       get_elapsed_timer2=get_elapsed_system_timer2(t)
    else
       get_elapsed_timer2=get_elapsed_cpu_timer2(t)
    end if
  end function get_elapsed_timer2

  real(8) function get_mark_system_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) then
       get_mark_system_timer2=get_mark_timer(t%sys)
    else
       write(*, *)
       write(*, '("Error(timing:get_mark_system_timer2): no System clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_mark_system_timer2

  real(8) function get_mark_cpu_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tcpu) then
       get_mark_cpu_timer2=get_mark_timer(t%cpu)
    else
       write(*, *)
       write(*, '("Error(timing:get_mark_cpu_timer2): no CPU clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_mark_cpu_timer2

  real(8) function get_mark_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsysdefault.and.t%tsys) then
       get_mark_timer2=get_mark_system_timer2(t)
    else
       get_mark_timer2=get_mark_cpu_timer2(t)
    end if
  end function get_mark_timer2

  integer function get_intervals_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) get_intervals_timer2=get_intervals_timer(t%sys)
    if (t%tcpu) get_intervals_timer2=get_intervals_timer(t%cpu)
  end function get_intervals_timer2

  real(8) function get_intlast_system_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) then
       get_intlast_system_timer2=get_intlast_timer(t%sys)
    else
       write(*, *)
       write(*, '("Error(timing:get_intlast_system_timer2): no System clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intlast_system_timer2

  real(8) function get_intlast_cpu_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tcpu) then
       get_intlast_cpu_timer2=get_intlast_timer(t%cpu)
    else
       write(*, *)
       write(*, '("Error(timing:get_intlast_cpu_timer2): no CPU clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intlast_cpu_timer2

  real(8) function get_intlast_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsysdefault.and.t%tsys) then
       get_intlast_timer2=get_intlast_system_timer2(t)
    else
       get_intlast_timer2=get_intlast_cpu_timer2(t)
    end if
  end function get_intlast_timer2

  real(8) function get_intsmall_system_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) then
       get_intsmall_system_timer2=get_intsmall_timer(t%sys)
    else
       write(*, *)
       write(*, '("Error(timing:get_intsmall_system_timer2): no System clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intsmall_system_timer2

  real(8) function get_intsmall_cpu_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tcpu) then
       get_intsmall_cpu_timer2=get_intsmall_timer(t%cpu)
    else
       write(*, *)
       write(*, '("Error(timing:get_intsmall_cpu_timer2): no CPU clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intsmall_cpu_timer2

  real(8) function get_intsmall_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsysdefault.and.t%tsys) then
       get_intsmall_timer2=get_intsmall_system_timer2(t)
    else
       get_intsmall_timer2=get_intsmall_cpu_timer2(t)
    end if
  end function get_intsmall_timer2

  real(8) function get_intlarge_system_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) then
       get_intlarge_system_timer2=get_intlarge_timer(t%sys)
    else
       write(*, *)
       write(*, '("Error(timing:get_intlarge_system_timer2): no System clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intlarge_system_timer2

  real(8) function get_intlarge_cpu_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tcpu) then
       get_intlarge_cpu_timer2=get_intlarge_timer(t%cpu)
    else
       write(*, *)
       write(*, '("Error(timing:get_intlarge_cpu_timer2): no CPU clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intlarge_cpu_timer2

  real(8) function get_intlarge_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsysdefault.and.t%tsys) then
       get_intlarge_timer2=get_intlarge_system_timer2(t)
    else
       get_intlarge_timer2=get_intlarge_cpu_timer2(t)
    end if
  end function get_intlarge_timer2

  real(8) function get_intmean_system_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsys) then
       get_intmean_system_timer2=get_intmean_timer(t%sys)
    else
       write(*, *)
       write(*, '("Error(timing:get_intmean_system_timer2): no System clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intmean_system_timer2

  real(8) function get_intmean_cpu_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tcpu) then
       get_intmean_cpu_timer2=get_intmean_timer(t%cpu)
    else
       write(*, *)
       write(*, '("Error(timing:get_intmean_cpu_timer2): no CPU clock &
	&timer running")')
       write(*, *)
       stop
    end if
  end function get_intmean_cpu_timer2

  real(8) function get_intmean_timer2(t)
    implicit none
    type(timer2) :: t
    if (t%tsysdefault.and.t%tsys) then
       get_intmean_timer2=get_intmean_system_timer2(t)
    else
       get_intmean_timer2=get_intmean_cpu_timer2(t)
    end if
  end function get_intmean_timer2

end module modtimer2
