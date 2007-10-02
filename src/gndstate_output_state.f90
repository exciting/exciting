subroutine gndstate_output_state
use modmain
use modmpi
logical exist
  if(rank.eq.0) then
        call writeengy(60)
        write(60,*)
        write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
        write(60,'(" (states/Hartree/spin/unit cell)")')
        ! write total energy to TOTENERGY.OUT and flush
        write(61,'(G18.10)') engytot
        call flushifc(61)
        ! write DOS at Fermi energy to FERMIDOS.OUT and flush
        write(62,'(G18.10)') fermidos
        call flushifc(62)
        ! output charges and moments
        call writechg(60)
        ! write total moment to MOMENT.OUT and flush
        if (spinpol) then
           write(63,'(3G18.10)') momtot(1:ndmag)
           call flushifc(63)
        end if
        ! output effective field for fixed spin moment calculations
        if (fixspin) then
           write(60,*)
           write(60,'("FSM effective field      : ",3G18.10)') bfsmc(1:ndmag)
        end if
        ! check for WRITE file
        inquire(file='WRITE',exist=exist)
        if (exist) then
           write(60,*)
           write(60,'("WRITE file exists - writing STATE.OUT")')
           call writestate
           open(50,file='WRITE')
           close(50,status='DELETE')
        end if
        ! write STATE.OUT file if required
        if (nwrite.ge.1) then
           if (mod(iscl,nwrite).eq.0) then
              call writestate
              write(60,*)
              write(60,'("Wrote STATE.OUT")')
           end if
        end if
        ! exit self-consistent loop if last iteration is complete
     endif

end subroutine