!BOP
!
! !ROUTINE: mrqcof
!
! !INTERFACE:
      subroutine mrqcof(x,y,ndata,a,ma,alpha,beta,nca,          &
     &                  chisq)

! !DESCRIPTION:
!
!  Used by \texttt{mrqmin} to evaluate the linearized 
!  fitting matrix \texttt{alpha}, and vector \texttt{beta} as:
!
! \begin{equation}\label{mrqcof-01}
! \begin{aligned}
!\beta_k &\equiv - \frac{1}{2}\frac{\partial\chi^2}{\partial a_k} &
!\alpha_{kl}&\equiv \frac{1}{2}\frac{\partial^2\chi^2}{\partial
!a_k\partial a_l} , 
!\end{aligned}
!\end{equation}
!and calculate $\chi^2$. 
!     
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: ndata ! Number of datapoints to be fitted
      integer(4), intent(in) :: ma    ! Number of coefficients
      integer(4), intent(in) :: nca   ! Size of working-space arrays

      real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
      complex(8),    intent(in) :: y(*)  ! data value of the datapoint

! !INPUT/OUTPUT PARAMETERS:
      
      real(8),    intent(inout) :: chisq
      
      real(8),    intent(inout) :: a(1:ma)
      real(8),    intent(inout) :: alpha(1:nca,1:nca)
      real(8),    intent(inout) :: beta(1:ma)
      
! !LOCAL VARIABLES:      

      integer(4) :: i 
      integer(4) :: j 
      integer(4) :: k 
      integer(4) :: l 
      integer(4) :: m 

      complex(8)    :: wt 
      complex(8)    :: dy 
      complex(8)    :: ymod 
      complex(8)    :: cx
      complex(8)    :: dyda(1:ma) 
      complex(8)    :: d2yda(1:ma,1:ma) 


! !EXTERNAL ROUTINES: 

      
      external ratfun
      
! !REVISION HISTORY:
!
! Original subroutine: mrqcof.for (c) copr. 1986-92 numerical recipes
! software &681i..
! Last modified: 7th. Jul. 2005 by RGA
!
!EOP
!BOC 
      do j=1,ma
!      
!       Initialize (symmetric) alpha, beta. 
!
        do k=1,j
          alpha(j,k)=0. 
        enddo ! k
        beta(j)=0. 
      enddo ! j
      chisq=0.0d0 
!      
!     Summation loop over all data. 
!
      do i=1,ndata 
        cx=cmplx(0.0d0,x(i),8)
        call ratfun(cx,a,ymod,dyda,d2yda,ma/4-1)
        dy=y(i)-ymod 
        do l=1,ma 
          wt=dyda(l) 
          do m=1,l 
            alpha(l,m)=alpha(l,m)+real(wt*conjg(dyda(m)))
!     &           real(dy*conjg(d2yda(l,m)))
          enddo ! m 
          beta(l)=beta(l)+real(dy*conjg(wt))
        enddo ! l
        chisq=chisq+dy*conjg(dy) ! And find \chi^2. 
      enddo ! i
      chisq=chisq/ndata
!      
!     Fill in the symmetric side. 
!
      do j=2,ma
        do k=1,j-1 
          alpha(k,j)=alpha(j,k) 
        enddo ! k
      enddo ! j 
      return 
      end subroutine mrqcof
!EOC  
