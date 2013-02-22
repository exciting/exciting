!BOP
!
! !ROUTINE: gw_main
!
! !INTERFACE:
subroutine task_band
!
! !DESCRIPTION:
!
! This subroutine calculate gw band structure by 
! interpolating the qp energy corrections
! 
! !USES:

    use modinput
    use modmain
    use modgw
    
!
! !DEFINED PARAMETERS:

    implicit none
!
! !LOCAL VARIABLES:

    integer(4) :: ik, ib
    real(8)    :: tstart, tend
    complex(8), allocatable :: de1(:,:), de2(:,:)
!
! !INTRINSIC ROUTINES:
    intrinsic cpu_time

!
! !EXTERNAL ROUTINES: 

!
! !REVISION HISTORY:
!       
! Created July 2011 by DIN
!
!EOP
!BOC
      call cpu_time(tstart)

!     read QP energies from file and perform Fourier interpolation (if required)
      call getevalqp(nkpt)

!     Write QP bandstructure to disk
      open(50,file='BAND-QP.OUT',action='WRITE',form='FORMATTED')
      do ib = ibgw, min(nbgw,nstsv)
        do ik = 1, nkpt
           write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)-efermi
        end do !ik
        write(50,*)
      end do !ib
      close(50)

      call cpu_time(tend)
      if (debug) call write_cputime(55,tend-tstart, 'TASK_BAND')
 
      return
end subroutine
!!EOC
