subroutine gndstate_output_timing_information
use modmain
use modmpi
 if(rank.eq.0) then
	   write(60,*)
     write(60,'("Timings (CPU seconds) :")')
     write(60,'(" initialisation                        : ",F12.2)') timeinit
     write(60,'(" Hamiltonian and overlap matrix set up : ",F12.2)') timemat
     write(60,'(" first-variational secular equation    : ",F12.2)') timefv
     if (spinpol) then
        write(60,'(" second-variational calculation        : ",F12.2)') timesv
     end if
     write(60,'(" charge density calculation            : ",F12.2)') timerho
     write(60,'(" potential calculation                 : ",F12.2)') timepot
     if (tforce) then
        write(60,'(" force calculation                     : ",F12.2)') timefor
     end if
     timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
     write(60,'(" total                                 : ",F12.2)') timetot
     write(60,*)
     write(60,'("+----------------------------------+")')
     write(60,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3," stopped |")') version
     write(60,'("+----------------------------------+")')
 endif
end subroutine