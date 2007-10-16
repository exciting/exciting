
subroutine tdepilog
  use modmain
  use modtddft
  use modmpi
  use m_filedel
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tdepilog'
  character(10) dat, tim
  real(8) :: cput,wallt,cputcum,walltcum

  ! finalize global counters
  call date_and_time(date=dat,time=tim)
  call cpu_time(cputim0f)
  call system_clock(COUNT=systim0f)
  cput = cputim0f - cputim0i
  wallt = dble(systim0f - systim0i)/dble(cntrate)
  cputcum = cputim0f - cputimcum
  walltcum = dble(systim0f - systimcum)/dble(cntrate)

  write(unitout,'(a,i8,a)') '*** Info('//thisnam//'): task Nr.', task, &
       ' stopped gracefully'
  write(unitout,*)
  write(unitout,'(a)') 'Timings [seconds/minutes/hours/days]'
  write(unitout,'(a)') '  Date (YYYY-MM-DD) : '//dat(1:4)//'-'//dat(5:6)//'-' &
       //dat(7:8)
  write(unitout,'(a)') '  Time (hh:mm:ss)   : '//tim(1:2)//':'//tim(3:4)//':' &
       //tim(5:6)
  write(unitout,'(a,g18.6)') '  CPU load              :', cput/wallt*100
  write(unitout,'(a,4g18.6)') '  CPU time              :', &
       cput,cput/60,cput/3600,cput/(24*3600)
  write(unitout,'(a,4g18.6)') '  wall time             :', &
       wallt,wallt/60,wallt/3600,wallt/(24*3600)
  write(unitout,'(a,g18.6)') '  CPU load (cumulative) :', cputcum/walltcum*100
  write(unitout,'(a,4g18.6)') '  CPU time (cumulative) :', &
       cputcum,cputcum/60,cputcum/3600,cputcum/(24*3600)
  write(unitout,'(a,4g18.6)') '  wall time (cumulative):', &
       walltcum,walltcum/60,walltcum/3600,walltcum/(24*3600)
  write(unitout,*)
  write(unitout,'(a)') '============ TDDFT@EXCITING stopped ==================&
       &======================'
  write(unitout,*)
  close(unitout)

  ! remove tag
!!$  if (.not.tresume) call filedel(trim(fnresume))
  call filedel(trim(fnresume))

end subroutine tdepilog
