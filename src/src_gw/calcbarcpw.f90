!BOP
!
! !ROUTINE: calcbarcpw
!
! !INTERFACE:
      subroutine calcbarcpw(iq)
      
! !DESCRIPTION:
!
! This subroutine calculates the bare coulomb potential matrix by
! expanding it in plane waves and then transforming to the mixed basis.
! (For test purposes only)      
!
! !USES:
      
      use modmain
      use modgw
!
! !INPUT PARAMETERS: 

      implicit none
      
      integer(4), intent(in) :: iq ! index of the q-point

!
! !LOCAL VARIABLES:
      
      integer(4) :: i, ipw, ipin

      real(8) :: t1,t2,t3
      real(8) :: exev, qg1len
      real(8), dimension(3) :: qvec
      real(8), dimension(3) :: qg1
      complex(8), allocatable :: mat1(:,:), mat2(:,:)
    
!EOP
!BOC      
      call cpu_time(t1)
      if(t1.lt.0.0d0)write(fgw,*)'warning, t1 < 0'

      if(allocated(barc))deallocate(barc)
      allocate(barc(matsiz,matsiz))
      if(allocated(sqbarc))deallocate(sqbarc)
      allocate(sqbarc(matsiz,matsiz))
      barc=zzero
      sqbarc=zzero

      allocate(mat1(1:matsiz,1:ngbarc(iq)))
      allocate(mat2(1:matsiz,1:ngbarc(iq)))
      mat1=zzero
      mat2=zzero
      
      qvec(1:3)=vqc(1:3,iq)

      call calcmpwmix(iq)
      call cpu_time(t2)
      if(t2.lt.0.0d0)write(fgw,*)'warning, t2 < 0'
      
      ipin=1
      if (Gamma) ipin=2
      
      do ipw=ipin,ngbarc(iq)
!
!       Calculate the G vector          
! 
        qg1(1:3)=qvec(1:3)+vgc(1:3,igqigb(ipw,iq))
        qg1len=qg1(1)*qg1(1)+qg1(2)*qg1(2)+qg1(3)*qg1(3)
        exev=4.0d0*pi/qg1len
        do i=1,matsiz
          mat1(i,ipw)=mpwmix(i,ipw)*exev
          mat2(i,ipw)=mpwmix(i,ipw)*sqrt(exev)
        enddo
      enddo    

      call zgemm('n','c',matsiz,matsiz,ngbarc,zone,mat1,matsiz,mpwmix,  &
     &            matsiz,zzero,barc,matsiz)
     
      call zgemm('n','c',matsiz,matsiz,ngbarc,zone,mat2,matsiz,mpwmix,  &
     &            matsiz,zzero,sqbarc,matsiz)
     
      deallocate(mat1,mat2)
      
      call cpu_time(t3)
      if(t3.lt.0.0d0)write(fgw,*)'warning, t3 < 0'
      call write_cputime(fgw,t3-t1, 'CALCBARCPW')
      call write_cputime(fgw,t2-t1, 'CALCMPWMIX')
      
      end subroutine calcbarcpw
!EOC            
     
               
