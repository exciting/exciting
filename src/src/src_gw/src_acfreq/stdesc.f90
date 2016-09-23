!BOP
!
! !ROUTINE: stdesc
!
! !INTERFACE: 
      subroutine stdesc(x,y,ndata,a,ma,varsq)

! !DESCRIPTION:
!
! This program perform one iteration of the steepest descent method.
!
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: ndata ! Number of datapoints to be fitted
      integer(4), intent(in) :: ma    ! Number of coefficients

      real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
      complex(8),    intent(in) :: y(*)  ! data value of the datapoint

! !INPUT/OUTPUT PARAMETERS:
      
      real(8),    intent(inout) :: varsq
      
      real(8),    intent(inout) :: a(1:ma)
      
      
! !LOCAL VARIABLES:

      integer(4) :: j
      integer(4) :: it
      integer(4) :: i
      
      real(8) :: varsqold
      real(8) :: num
      real(8) :: denomi
      real(8) :: denom
      real(8) :: lambda
      real(8) :: deltach
      
      real(8) :: atemp(1:ma)
      real(8), allocatable :: beta(:)
      real(8), allocatable :: alpha(:,:)
      
      logical :: conv


! !DEFINED PARAMETERS: 

      real(8), parameter :: tolstd = 1.0d-3


! !EXTERNAL ROUTINES: 


      external mrqcof
      
! !REVISION HISTORY:
!
! Created: 11th. Jul. 2005 by RGA
!
 
!EOP
!BOC
      allocate(beta(1:ma))
      allocate(alpha(1:ma,1:ma))
      varsqold=0.0d0
      conv=.false.
      it=0
      do while (.not.conv)
        it=it+1
        call mrqcof(x,y,ndata,a,ma,alpha,beta,ma,varsq)
        num=0.0d0
        denom=0.0d0
        do i=1,ma
          num=num+beta(i)*beta(i)
          denomi=0.0d0
          do j=1,ma
            denomi=denomi+alpha(i,j)*beta(j)
          enddo
          denom=denom+denomi*beta(i)
        enddo
        lambda=num/denom
        do j=1,ma 
          a(j)=a(j)+lambda*beta(j)
        enddo    
        deltach=abs(varsq-varsqold)
        varsqold=varsq
!        write(6,3)it,varsq,lambda,deltach
        atemp(1:ma)=a(1:ma)
        conv=((deltach.lt.tolstd).or.(it.gt.200))
      enddo
      if(it.gt.200)write(*,*)'WARNING: stepest descent did not converge &
     & after',it, 'iterations'
      deallocate(beta)
      deallocate(alpha)

!    3 format(' iter =',i4,' chi^2 = ',g18.10,'lambda = ',g18.10,        &
!     &       ' delta(chi^2) = ',g18.10)
      
      return
      
      end subroutine stdesc
!EOC        
      
