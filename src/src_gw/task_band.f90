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
    integer(4) :: np, ip
    
    real(8) :: emin, emax
    real(8) :: tstart, tend
    real(8) :: e, dos0, dos1
    
    complex(8), allocatable :: dek1(:,:), dek2(:,:)
    
    integer(4) :: i1, i2
    integer(4) :: nsym, isym, lspl
    integer, allocatable :: symc(:,:,:)
    
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
      call readeqp
!
!     Calculate the correction to the energies
!
      allocate(dek1(nkp1,ibgw:nbgw),dek2(nkp2,ibgw:nbgw))

      write(6,*)
      write(6,*)' band: calc bandstruct using Fourier interpolation'
      write(6,*)

      call boxmsg(6,'-',"BAND STRUCTURE SUMMARY")

      do ik = 1, nkp1
         dek1(ik,:)=cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
      enddo

!--------------------------------------------------
!     Fourier interpolation
!--------------------------------------------------
      dek2(:,:)=zzero
      call fourintp(dek1,nkp1,kvecs1,dek2,nkp2,kvecs2,nbandsgw)

      allocate(eqp2(ibgw:nbgw,nkp2))
      do ib = ibgw, nbgw
         do ik = 1, nkp2
            eqp2(ib,ik)=eks2(ib,ik)+real(dek2(ik,ib))
         enddo 
      enddo

!--------------------------------------------------
!     Write the bandstructure to disk
!--------------------------------------------------
      open(50,file='BAND-QP.OUT',action='WRITE',form='FORMATTED')
      do ib = ibgw, nbgw
        do ik = 1, nkp2
           write(50,'(2G18.10)') dpp1d(ik), eqp2(ib,ik)
        end do !ik
        write(50,*)
      end do !ib
      close(50)

!--------------------------------------------------
!     Initialize k-point grid for BZ integration
!--------------------------------------------------
     
      call init0
     
      nkpt = input%groundstate%ngridk(1)* &
     &       input%groundstate%ngridk(2)* &
     &       input%groundstate%ngridk(3)
      ntet = 6*nkpt

!     allocate the reduced k-point set arrays
      if (allocated(ivk)) deallocate (ivk)
      allocate (ivk(3,nkpt))
      if (allocated(vkl)) deallocate (vkl)
      allocate (vkl(3,nkpt))
      if (allocated(vkc)) deallocate (vkc)
      allocate (vkc(3,nkpt))
      if (allocated(wkpt)) deallocate (wkpt)
      allocate (wkpt(nkpt))
      if (allocated(indkp)) deallocate(indkp)
      allocate(indkp(nkpt))
      if (allocated(iwkp)) deallocate(iwkp)
      allocate(iwkp(nkpt))
      if (allocated(wtet)) deallocate(wtet)
      allocate(wtet(ntet))
      if (allocated(tnodes)) deallocate(tnodes)
      allocate(tnodes(4,ntet))

      nsym=1
      If (input%groundstate%reducek) nsym=nsymcrys
!             
!     get rotational part of crystal symmetries
!
      Allocate (symc(3,3,nsym))
      Do isym = 1, nsym
        lspl = lsplsymc(isym)
!       transpose of rotation for use with the library
        Do i1 = 1, 3
           Do i2 = 1, 3
              symc(i1,i2,isym) = symlat(i2,i1,lspl)
           End Do
        End Do
      End Do

!     suppress debug output in tetrahedron integration library (0)
      call tetrasetdbglv(0)
!
      call factorize(3,input%groundstate%vkloff,ikloff,dkloff)
!
      call kgen(bvec,nsym,symc,input%groundstate%ngridk,ikloff,dkloff,&
     &     nkpt,ivk,dvk,indkp,iwkp,ntet,tnodes,wtet,tvol,mnd)
!
!     Calculate DOS
!      
      open(50,file='TDOS-KS-QP.OUT',action='WRITE',form='FORMATTED')

      emin=minval(eqp1(ibgw:nbgw,:))
      emax=maxval(eqp1(ibgw:nbgw,:))

      np=400
      do ip = 1, np 
        e=emin+dble(ip-1)*(emax-emin)/(np-1)
        dos0=dostet(nkpt,nbandsgw,eks1(ibgw:nbgw,:),ntet,tnodes(1:4,1:ntet),wtet(1:ntet),tvol,e)
        dos1=dostet(nkpt,nbandsgw,eqp1(ibgw:nbgw,:),ntet,tnodes(1:4,1:ntet),wtet(1:ntet),tvol,e)
        write(50,'(3G18.10)') e, dos0, dos1
      enddo 

      close(50)
      
      write(fgw,*)
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart, 'TASK_BAND')
 
      deallocate(eks1,eks2,eqp1,eqp2,dek1,dek2)
      deallocate(symc)
      deallocate(ivk,vkl,vkc,wkpt,indkp,iwkp,wtet,tnodes)

end subroutine
!!EOC
