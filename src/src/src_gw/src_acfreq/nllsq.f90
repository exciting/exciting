!BOP
!
! !ROUTINE: nllsq
!
! !INTERFACE: 
      subroutine nllsq(x,y,ndata,a,ma,chisq)

! !DESCRIPTION:
!
! This subroutine fits the set for data points 
!\texttt{x(1:ndata)}, \texttt{y(1:ndata)}, and a nonlinear function
!dependent on \texttt{ma} coefficients \texttt{a(1:ma)} using the
!Levenberg-Marquardt method.
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: ndata ! Number of datapoints to be fitted
      integer(4), intent(in) :: ma    ! Number of coefficients
      real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
      complex(8), intent(in) :: y(*)  ! data value of the datapoint

! !INPUT/OUTPUT PARAMETERS:
      real(8),    intent(inout) :: a(1:ma)
      
! !OUTPUT PARAMETERS:
      real(8),    intent(inout) :: chisq
      
! !LOCAL VARIABLES:      
 
      integer(4) :: it

      real(8) :: alambda
      real(8) :: alambdaold
      real(8) :: chisqold
      real(8) :: deltach
      real(8) :: covar(1:ma,1:ma)
      real(8) :: alpha(1:ma,1:ma)

      logical :: converg

! !DEFINED PARAMETERS: 

      real(8), parameter :: tol = 1.0d-10


! !EXTERNAL ROUTINES: 

      
      external mrqmin
      external stdesc

! !REVISION HISTORY:
!
! Created: 8th. Jul. 2005 by RGA
!
!EOP
!BOC

      converg = .false.
      call stdesc(x,y,ndata,a,ma,chisqold)
      alambda = -0.1d0
      alambdaold=1.0d0
      it=0
      do while (.not.converg)
        it=it+1
        call mrqmin(x,y,ndata,a,ma,covar,alpha,ma,chisq,alambda)
        deltach=abs(chisq-chisqold)
        if(alambda.lt.alambdaold)then
          converg = (deltach.lt.tol).or.(it.gt.1000)
          chisqold=chisq
        else  
          converg = (it.gt.1000)
        endif  
        alambdaold=alambda
      enddo
      if(it.gt.1000)write(*,*)'WARNING: chisq did not converge'
      alambda=0.0d0
      return
!    3 format(' iter =',i4,' chi^2 = ',g18.10,' delta(chi^2) = ',g18.10)
      end subroutine nllsq
      
!EOC      
