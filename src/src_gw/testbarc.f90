!BOP
!
! !ROUTINE: testbarc
!
! !INTERFACE:
      subroutine testbarc(iq)

! !DESCRIPTION:
!
! This subroutine calculates the matrix of the sqare root of the bare coulomb
! potential. The resulting matrix is stored back in barc.
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

      integer(4) :: i,j  ! just some counters

      integer(4) :: info   ! output status of the diagonalization
!                            subroutine       
      integer(4) :: lwork  ! size of the workspace array work
      integer(4) :: rwsize ! size of the workspace array rwork
      integer(4) :: gap
      
      real(8) :: qg1len, kk, t1,t2, x
      real(8), dimension(3) :: qg1, qvec, gvec
      real(8), allocatable :: rwork(:)     ! workspace array for the
!                                            diagonalization subroutine
      real(8), allocatable :: ev(:)        ! eigenvalues of barc in pw
      real(8), allocatable :: exev(:)      ! exact eigenvalues of barc

      complex(8), allocatable :: vtemp(:,:) ! temporary storage
!                                             for the eigenvectors
!                                             of barc
      complex(8), allocatable :: vsq(:,:) ! the square root of barc
      complex(8), allocatable :: work(:)  ! workspace array for
!                                           the diagonalization
!                                           subroutine

      complex(8), allocatable :: barcMB(:,:), sqbarcMB(:,:)
     
      logical :: done
 
! !EXTERNAL ROUTINES: 

      external calcvsq
      external errflg
      external outerr
      external zheev
!
! !INTRINSIC ROUTINES: 
!
      intrinsic abs
      intrinsic conjg
      intrinsic sqrt
      intrinsic cpu_time
!
! !REVISION HISTORY:
! 
! Created 4th. Jan. 2005 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC

!------------------------------------
!     Calculate barc in MB
!------------------------------------
      call cpu_time(t1)
      call calcbarcmb(iq)
      call cpu_time(t2)
      write(fgw,*)'calcbarcns -->',t2-t1
!
!     Save the barc matrices obtained from MB
!      
      allocate(barcMB(matsiz,matsiz),sqbarcMB(matsiz,matsiz))
      barcMB(:,:)=barc(:,:)
      sqbarcMB(:,:)=sqbarc(:,:)
      
      allocate(vsq(matsiz,matsiz))      
      call zgemm('n','n',matsiz,matsiz,matsiz,zone,sqbarc,matsiz,  &
     &              sqbarc,matsiz,zzero,vsq,matsiz)   

      write(99,*)
      write(99,*)'# i j     barc(i,j)         barc(j,i)'
      do i=1,matsiz,matsiz/10
        do j=1,matsiz,matsiz/10
          write(99,'(2i4,2f16.6,2f16.6)') i, j, real(barc(i,j)), aimag(barc(i,j)), real(barc(j,i)), aimag(barc(j,i))
        enddo  
      enddo
      write(99,*)
      write(99,*)'# i j     Re(vsq)       Re(barc)        Re(vsq-barc)'
      do i=1,matsiz,matsiz/10
        do j=1,matsiz,matsiz/10
          write(99,13)i,j,real(vsq(i,j)),real(barc(i,j)),real(vsq(i,j)-barc(i,j))
        enddo  
      enddo
      write(99,*)
      write(99,*)'# i j     Im(vsq)       Im(barc)        Im(vsq-barc)'
      do i=1,matsiz,matsiz/10
        do j=1,matsiz,matsiz/10
          write(99,13)i,j,aimag(vsq(i,j)),aimag(barc(i,j)),aimag(vsq(i,j)-barc(i,j))
        enddo  
      enddo
      deallocate(vsq)

!------------------------------------
!       Calculate barc from PW
!------------------------------------

      call cpu_time(t1)
      call calcbarcpw(iq)
      call cpu_time(t2)
      write(fgw,*)'calcbarcpw -->',t2-t1
!       
!     Set up the workspace for the diagonalization subroutine
!
      allocate(vtemp(matsiz,matsiz))
      allocate(ev(matsiz))

      lwork=2*matsiz
      rwsize=3*matsiz
      allocate(work(lwork))
      allocate(rwork(rwsize))

! --- Check barc
      vtemp(:,:)=barc(1:matsiz,1:matsiz)
      call zheev('v','u',matsiz,vtemp,matsiz,ev,work,lwork,rwork,info)
      if(info.ne.0)then
        write(*,*)'zheev, info =',info
        stop 'error in testbarc'
      endif    

      write(94,*) "### barc from pw ###"
      do i=1,matsiz,matsiz/10
        write(94,*)'eigenvector:',i,'eval =',ev(i)
        do j=1,matsiz,matsiz/10
           write(94,11)j,vtemp(j,i)
        enddo
        write(94,*)
      enddo
      
      write(95,*) "### Difference (real) between two different barc"
      write(96,*) "### Difference (imag) between two different barc"
      write(95,*)'iq =',iq     
      write(96,*)'iq =',iq     
      write(95,*)
      write(96,*)
      write(95,*)'# barc MB vs PW: i j    Re(barcMB)   Re(barcPW)     Re(barcMB-barcPW)'
      write(96,*)'# barc MB vs PW: i j    Im(barcMB)   Im(barcPW)     Im(barcMB-barcPW)'
      do i=1,matsiz,matsiz/10
        do j=1,matsiz,matsiz/10
          x=real(barcMB(i,j))*real(barc(i,j))
          if((x.lt.0.0d0).and.(abs(x).gt.1.0d-10))then
            write(95,16)i,j,real(barcMB(i,j)),real(barc(i,j)),         &
     &                   real(barcMB(i,j)-barc(i,j))
          else
            write(95,13)i,j,real(barcMB(i,j)),real(barc(i,j)),         &
     &                   real(barcMB(i,j)-barc(i,j))
          endif
          x=aimag(barcMB(i,j))*aimag(barc(i,j))
          if((x.lt.0.0d0).and.(abs(x).gt.1.0d-10))then
            write(96,16)i,j,aimag(barcMB(i,j)),aimag(barc(i,j)),         &
     &                   aimag(barcMB(i,j)-barc(i,j))
          else
            write(96,13)i,j,aimag(barcMB(i,j)),aimag(barc(i,j)),         &
     &                   aimag(barcMB(i,j)-barc(i,j))
          endif
        enddo  
      enddo

!================================================!
!       Calculate the exact eigenvalues 
!================================================!
      allocate(exev(matsiz))
      exev(:)=0.d0
      if (Gamma) then
        do i=1,matsiz
          gvec(1:3)=vgc(1:3,igqigb(matsiz-i+2,iq))
          qg1len=gvec(1)*gvec(1)+gvec(2)*gvec(2)+gvec(3)*gvec(3)
          exev(i)=4.0d0*pi/qg1len
        enddo    
      else  
        qvec(1:3)=vqc(1:3,iq)
        do i=1,matsiz
!
!         Calculate the G vector          
!
          gvec(1:3)=vgc(1:3,igqigb(matsiz-i+1,iq))
!
!         Calculate q+G:
!
          qg1(1:3)=qvec(1:3)+gvec(1:3)
!
!         Calculate |q+G|^2:
!
          qg1len=qg1(1)*qg1(1)+qg1(2)*qg1(2)+qg1(3)*qg1(3)
          exev(i)=4.0d0*pi/qg1len
        enddo    
!
!       Sort the exev
!        
        gap=matsiz/2
        do while(gap.ge.1)
          done=.false.
          do while(.not.done)
            done=.true.
            do i=1,matsiz-gap
              if(exev(i).gt.exev(i+gap))then
                kk=exev(i)
                exev(i)=exev(i+gap)
                exev(i+gap)=kk
                done=.false.
              endif
            enddo
          enddo
          gap=gap/2
        enddo
      endif 
!
!     Compare 'mb', 'pw' and 'exact' barc
!
      write(97,*) "### Comparison of barc eigenvalues ###"
      do i=1,matsiz
        write(97,12)i,barcev(i),ev(i),exev(i)
      enddo  

      deallocate(work)
      deallocate(rwork)
      deallocate(vtemp)
      deallocate(ev)
      deallocate(exev)
    
   11 format(i4,'(',g18.10,',',g18.10,')')
   12 format(i4,3g18.10)
   13 format(2i4,3f18.10)
   14 format(2i4,2('(',g18.10,',',g18.10,')'))
   15 format(i4,'(',g18.10,',',g18.10,')')
   16 format(2i4,3f18.10,' err')

      return
         
      end subroutine testbarc
!EOC            
            
            
        
      
