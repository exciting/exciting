
subroutine tdinit
  use modmain
  use modtddft
  use modfxcifc
  use modpar
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'tdinit'
  character(10) dat, tim, s2
  integer :: un,i
  logical :: ex

  ! rank of process
  call getunit(un)

  ! check consistency of rank and nproc
  if ((nproc.lt.1).or.(nproc.gt.maxproc)) then
     write(*,*) 'Error('//trim(thisnam)//'): Error in parallel &
          &initialization: number of processes out of range:',nproc
     call terminate()
  end if

  if ((rank.gt.nproc).or.(rank.lt.1)) then
     write(*,*) 'Error('//trim(thisnam)//'): Error in parallel &
          &initialization: rank out of range:',rank
     call terminate()
  end if

  ! generate resuem file
  if (nproc.gt.1) then
     write(fnresume,'(".resume",i3.3)') rank
     write(spar,'("proc",i3.3)') rank
  else
     fnresume='.resume'
  end if

  ! initialize for first call to main routine
  if (calledtd.eq.1) resumechkpts=0
  ! read in checkpoint if present
  call resread(un,resumetask,resumechkpts,tresume)
  ! checkpointing starts
  if (.not.tresume) then
     resumechkpts(:,1)=0
     call resupd(un,task,resumechkpts,' : prolog')
  end if

  ! separate file extentsion
  tdfilext = '.OUT'
  tdfileout = 'TDINFO'
  if ((nproc.gt.1).and.(rank.gt.1)) tdfilext='_'//trim(spar)//'.OUT'
  if ((nproc.gt.1).and.(rank.gt.1)) tdfileout = '.'//trim(tdfileout)

  !initialize global counters
  call cpu_time(cputim0i)
  call system_clock(COUNT_RATE=cntrate)
  call system_clock(COUNT=systim0i)
  call date_and_time(date=dat,time=tim)
  if (calledtd.eq.1) call system_clock(COUNT=systimcum)
  call getunit(unitout)
  if ( tappinfo.or.(calledtd.gt.1)) then
     open(unitout, file=trim(tdfileout)//trim(tdfilext), action='write', &
          position='append')
  else
     open(unitout, file=trim(tdfileout)//trim(tdfilext), action='write', &
          status='replace')
  end if

  ! scaling factor for output of energies
  escale=1.d0
  if (tevout) escale=27.2117

  ! get exchange-correlation functional data
  call getfxcdata(fxctype,fxcdescr,fxcspin)

  ! write to info file
  write(unitout,'(a)')
  write(unitout,'(a)') '============ TDDFT@EXCITING started ==================&
       &======================'
  write(unitout,*)
  write(unitout,'("EXCITING version  : ",I1.1,".",I1.1,".",I3.3)') version
  write(unitout,'("TDDFT build       : ",I4.4                  )') &
       build_tddft
#ifdef MPI
  write(unitout,'(" compiled for MPI parallelization")') 
#endif
#ifndef MPI
  write(unitout,'(" compiled for pseudo parallelization ")') 
#endif
  write(unitout,'(a)') 'Date (YYYY-MM-DD) : '//dat(1:4)//'-'//dat(5:6)//'-'// &
       dat(7:8)
  write(unitout,'(a)') 'Time (hh:mm:ss)   : '//tim(1:2)//':'//tim(3:4)//':'// &
       tim(5:6)
  write(unitout,*)
  write(unitout,'(a,i6,a)') '*** Info('//trim(thisnam)//'): task Nr.', &
       task,' started'
  if ((nproc.gt.1).and.(rank.eq.1)) write(unitout,'(a,2i6)') '*** Info('// &
       trim(thisnam)//'): (parallel) master, rank/number of processes:',rank, &
       nproc
  if ((nproc.gt.1).and.(rank.ne.1)) write(unitout,'(a,2i6)') '*** Info('// &
       trim(thisnam)//'): (parallel) slave, rank/number of processes:',rank, &
       nproc
  if ((nproc.gt.1).and.(rank.ne.1)) write(unitout,'(a,i6)') '*** Info('// &
       trim(thisnam)//'): (parallel) baridl:',baridl
  if ((nproc.gt.1).and.(rank.ne.1)) write(unitout,'(a,i6)') '*** Info('// &
       trim(thisnam)//'): (parallel) baridl2:',baridl2
  if (tresume) write(unitout,'(a,i6,a,10i9)') '*** Info('//trim(thisnam)// &
       '): resuming task',resumetask,' from checkpoints',resumechkpts(:,1)
  if (notelns.gt.0) then
     write(unitout,*)
     write(unitout,'("Notes :")')
     do i=1,notelns
        write(unitout,'(A)') notes(i)
     end do
  write(unitout,*)
  end if
  call flushifc(unitout)

end subroutine tdinit
