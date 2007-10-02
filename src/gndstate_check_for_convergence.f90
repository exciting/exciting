subroutine gndstate_check_for_convergence(tlast_,tstop_,dv)
  use modmain
  use modmpi
  real(8),intent(in):: dv
  logical,intent(inout)::tlast_,tstop_
  logical ::exist
  if(rank.eq.0)then
     ! check for convergence
     if (iscl.ge.2) then
        write(60,*)
        write(60,'("RMS change in effective potential (target) : ",G18.10,&
             &" (",G18.10,")")') dv,epspot
        if (dv.lt.epspot) then
           write(60,*)
           write(60,'("Potential convergence target achieved")')
           tlast_=.true.
        end if
        if (xctype.lt.0) then
           write(60,'("Magnitude of OEP residue : ",G18.10)') resoep
        end if
     end if
     ! check for STOP file
     inquire(file='STOP',exist=exist)
     if (exist) then
        write(60,*)
        write(60,'("STOP file exists - stopping self-consistent loop")')
        tstop_=.true.
        tlast_=.true.
        open(50,file='STOP')
        close(50,status='DELETE')
     end if
     ! output the current total CPU time
     timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
     write(60,*)
     write(60,'("Time (CPU seconds) : ",F12.2)') timetot
     ! end the self-consistent loop
  endif
#ifdef MPI
  call MPI_bcast(tstop_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call MPI_bcast(tlast_,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine gndstate_check_for_convergence
