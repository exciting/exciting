!BOP
!
! !ROUTINE: task_band

! !INTERFACE:
subroutine task_band

!!DESCRIPTION:
!
! Calculate GW band structure by interpolating the QP energies obtained in GW cycle.
! 
!!USES:
    use modinput
    use modmain
    use modgw

!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ik, ib

!!REVISION HISTORY:
!       
! Created July 2011 by DIN
!
!EOP
!BOC

    if (.not.associated(input%properties%bandstructure)) then
      write(*,*)
      write(*,*) 'ERROR(task_band): Band structure path is not specified!'
      write(*,*) '                  Check your input file (properties%bandstructure).'
      stop
    end if

    ! read QP energies from file and perform Fourier interpolation (if required)
    call getevalqp(nkpt,vkl,evalsv)
      
    ! write QP bandstructure to disk
    open(50,file='BAND-QP.OUT',action='WRITE',form='FORMATTED')
    do ib = ibgw, min(nbgw,nstsv)
      do ik = 1, nkpt
        write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)
      end do !ik
      write(50,*)
    end do !ib
    close(50)
    
    return
end subroutine
!EOC
