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
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
!
!     Read the quasiparticle energies
!
      call linmsg(fgw,'-',"BAND STRUCTURE")

! ---- Old way (read data from QPENE-eV.OUT)

!     call readeqp
!     
!     allocate(de1(nkp1,ibgw:nbgw),de2(nkpt,ibgw:nbgw))
!     do ik = 1, nkp1
!        de1(ik,:)=cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
!     enddo

!     de2(:,:)=zzero
!     call fourintp(de1,nkp1,kvecs1,de2,nkpt,vkl,nbgw-ibgw+1)

!     do ib = ibgw, min(nbgw,nstsv)
!        do ik = 1, nkpt
!           evalsv(ib,ik)=evalsv(ib,ik)+real(de2(ik,ib))
!        enddo 
!     enddo

! ----

!     read QP energies and perform Fourier interpolation 
      call getevalqp(nkpt)

!     Write the bandstructure to disk
      open(50,file='BAND-QP.OUT',action='WRITE',form='FORMATTED')
      do ib = ibgw, min(nbgw,nstsv)
        do ik = 1, nkpt
           write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)-efermi
        end do !ik
        write(50,*)
      end do !ib
      close(50)

      write(fgw,*)
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'TASK_BAND')
 
      return
end subroutine
!!EOC
