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
    real(8) :: tstart, tend
    complex(8), allocatable :: dek1(:,:), dek2(:,:)
   
!
! !INTRINSIC ROUTINES:
    external  kgen
    real(8), external :: dostet
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

! ---- Old way
!      call readeqp
!      allocate(dek1(ibgw:nbgw,nkp1),dek2(ibgw:nbgw,nkpt))
!      do ik = 1, nkp1
!         dek1(ibgw:nbgw,ik)=cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
!      enddo
!      dek2(:,:)=zzero
!      call fourintp(nbgw-ibgw+1,nkp1,kvecs1,dek1,nkpt,vkl,dek2)
!      do ib = ibgw, min(nbgw,nstsv)
!         do ik = 1, nkpt
!            evalsv(ib,ik)=evalsv(ib,ik)+real(dek2(ib,ik))
!         enddo 
!      enddo
! ----

!     shift KS energies
      evalsv(:,:)=evalsv(:,:)-efermi

!     read QP energies and perform Fourier interpolation 
      call getevalqp(nkpt)

!     Write the bandstructure to disk
      open(50,file='BAND-QP.OUT',action='WRITE',form='FORMATTED')
      do ib = ibgw, min(nbgw,nstsv)
        do ik = 1, nkpt
           write(50,'(2G18.10)') dpp1d(ik), evalsv(ib,ik)
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
