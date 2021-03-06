!BOP
!
! !ROUTINE: task_epsev
!
! !INTERFACE:
      subroutine task_epsev
      
! !DESCRIPTION:
!
! 
!
! !USES:
      use modinput
      use modmain
      use modgw
      
! !LOCAL VARIABLES:

      implicit none

      integer(4) :: iq, iqp
      integer(4) :: i, iom
      
      real(8)    :: tstart,tend
      
      complex(8), allocatable :: eps(:,:)
      
      integer(4) :: info
      integer(4) :: lwork
      real(8), allocatable :: rwork(:)
      complex(8), allocatable :: work(:)
      complex(8), allocatable :: vl(:,:),vr(:,:)

      complex(8), allocatable :: eval(:), evalom(:,:)
      complex(8), allocatable :: ievom1(:,:), ievom2(:,:)
!
! !EXTERNAL ROUTINES: 
!
      external zgeev
      
! !INTRINSIC ROUTINES: 

      intrinsic cpu_time      
      
! !REVISION HISTORY:
!
! Created 17th. April 2004 by RGA
! Revisited: Jan 2012 by DIN
!
!EOP
!BOC            
      call cpu_time(tstart)
!
!     Q-point
!
      iq=input%gw%iik
      iqp=indkpq(iq,1)
!
!     Calculate the integration weights using the linearized tetrahedron method
!
      call kintw
!
!     Calculate the momentum matrix elements
!
      if(.not.input%gw%rpmat)then
        call calcpmat        ! <--- original (RGA's) version
        !call calcpmatgw     ! <--- modified xs (SAG's) version
      end if
!
!     Calculate the dielectric function matrix
!      
      if(allocated(epsilon))deallocate(epsilon)
      allocate(epsilon(matsizmax,matsizmax,nomeg))

      call calcepsilon(iqp)

      if(allocated(inveps))deallocate(inveps)
      allocate(inveps(matsizmax,matsizmax,nomeg))
      
      call calcinveps

!------------------------------------------------------------------------
!     Test the calculated dielectric function
!------------------------------------------------------------------------      

      lwork=64*matsiz
      allocate(work(lwork))
      allocate(rwork(lwork))
      allocate(eps(matsiz,matsiz))
      allocate(eval(matsiz))
      allocate(vl(matsiz,matsiz))
      allocate(vr(matsiz,matsiz))
      allocate(evalom(matsiz,1:nomeg))
      allocate(ievom1(matsiz,1:nomeg))
      allocate(ievom2(matsiz,1:nomeg))

      do iom = 1, nomeg
      
        eps(1:matsiz,1:matsiz) = epsilon(1:matsiz,1:matsiz,iom)
        
        call zgeev('n','n',matsiz,eps,matsiz,eval,vl,matsiz,vr, &
       &  matsiz,work,lwork,rwork,info)        
        if(info.ne.0)then
          write(*,*)'zgeev, info =',info
          stop
        endif
        
        ! \epsilon
        evalom(1:matsiz,iom)=eval(1:matsiz)

        ! \frac{1}{\epsilon}
        do i=1,matsiz
          if(abs(eval(i)).gt.1.0d-10) then
            ievom1(i,iom)=zone/eval(i)
          else
            ievom1(i,iom)=zzero
          end if
        enddo  

        eps(1:matsiz,1:matsiz) = inveps(1:matsiz,1:matsiz,iom)

        call zgeev('n','n',matsiz,eps,matsiz,eval,vl,matsiz,vr, &
       &  matsiz,work,lwork,rwork,info)        
        if(info.ne.0)then
          write(*,*)'zgeev, info =',info
          stop
        endif    
        
        ! \epsilon^{-1}
        ievom2(1:matsiz,iom)=eval(1:matsiz)
      
      end do

      deallocate(eps)
      deallocate(epsilon)
      deallocate(inveps)
      
      deallocate(work)
      deallocate(rwork)
      deallocate(eval)
      deallocate(vl)
      deallocate(vr)

!---------------------------------------
!     Output
!---------------------------------------

!    # Eigenvalues of the dielectric functions
!    # omeg [eV]  \epsilon   \frac{1}{\epsilon}   \epsilon^{-1}

      open(unit=72,file='epsevr.out',form='formatted',status='unknown')
      open(unit=73,file='epsevi.out',form='formatted',status='unknown')

      write(72,*)'# Eigenvalues of the dielectric functions (real part)'
      write(72,*)'# omeg [eV]  \epsilon   \frac{1}{\epsilon}   \epsilon^{-1}'
      write(73,*)'# Eigenvalues of the dielectric functions (imaginary part)'
      write(73,*)'# omeg [eV]  \epsilon   \frac{1}{\epsilon}   \epsilon^{-1}'
      
      do i=1,matsiz,10
        write(72,*)'# \epsilon_{i}, i=', i
        write(73,*)'# \epsilon_{i}, i=', i
        do iom=1,nomeg
          write(72,10)freqs(iom)*hev,real(evalom(i,iom)),real(ievom1(i,iom)),real(ievom2(i,iom))
          write(73,10)freqs(iom)*hev,aimag(evalom(i,iom)),aimag(ievom1(i,iom)),aimag(ievom2(i,iom))
        enddo
        write(71,*)
        write(72,*)
      enddo
   10 format(4(f18.8,' '))      
      close(72)
      close(73)

      deallocate(evalom)
      deallocate(ievom1)
      deallocate(ievom2)

      call cpu_time(tend)
      call write_cputime(fgw,tend-tstart,'TASK_EPSEV')

      return
      end subroutine
!EOC
