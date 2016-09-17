!BOP
!
! !ROUTINE: mrqmin
!
! !INTERFACE:
      subroutine mrqmin(x,y,ndata,a,ma,covar,alpha,nca,chisq,alamda)
     
! !DESCRIPTION:
!     
!Levenberg-Marquardt method, attempting to reduce the value
!$\chi^2$ of a fit between a set of data points
!\texttt{x(1:ndata)}, \texttt{y(1:ndata)}, and a nonlinear function
!dependent on \texttt{ma} coefficients \texttt{a(1:ma)}. The program 
!returns current best-fit values for the
!parameters \texttt{a(1:ma)}, and $\chi^2 = \texttt{chisq}$. The
!arrays \texttt{covar(1:nca,1:nca)}, \texttt{alpha(1:nca,1:nca)}
!with physical dimension \texttt{nca} ($\geq$ the number of fitted
!parameters) are used as working space during most iterations.
! On the first call provide an initial guess for the
!parameters \texttt{a}, and set $\texttt{alamda}<0$ for
!initialization (which then sets $\texttt{alamda}=.001$). If
!a step succeeds \texttt{chisq} becomes smaller and \texttt{alamda}
!decreases by a factor of 10. If a step fails \texttt{alamda} grows
!by a factor of 10. You must call this routine repeatedly until
!convergence is achieved. Then, make one final call with
!$\texttt{alamda}=0$, so that \texttt{covar(1:ma,1:ma)} returns the
!covariance matrix, and \texttt{alpha} the curvature matrix.
!(Parameters held fixed will return zero covariances.)
!
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: ndata ! Number of datapoints to be fitted
      integer(4), intent(in) :: ma    ! Number of coefficients
      integer(4), intent(in) :: nca   ! Size of working-space arrays

      real(8),    intent(in) :: x(*)  ! abscisas of the datapoints
      complex(8),    intent(in) :: y(*)  ! data value of the datapoint

! !INPUT/OUTPUT PARAMETERS:
      
      real(8),    intent(inout) :: alamda
      real(8),    intent(inout) :: chisq
      
      real(8),    intent(inout) :: a(1:ma)
      real(8),    intent(inout) :: covar(1:nca,1:nca)
      real(8),    intent(inout) :: alpha(1:nca,1:nca)
      
! !LOCAL VARIABLES:      
            
      integer(4) :: j
      integer(4) :: k
      integer(4) :: l
      integer(4), parameter :: max = 50
      
      real(8)    :: ochisq

      real(8), dimension(1:max) :: atry
      real(8), dimension(1:max) :: beta
      real(8), dimension(1:max) :: da

      save ochisq,atry,beta,da
      

! !EXTERNAL ROUTINES: 

      
      external covsrt
      external gaussj
      external mrqcof
      
! !REVISION HISTORY:
!
! Original subroutine: mrqmin.for (c) copr. 1986-92 numerical recipes
! software &680i..
! Last modified: 7th. Jul. 2005 by RGA
!
 
!EOP
!BOC
      if(alamda.lt.0.0d0)then ! Initialization.
        alamda=1.0d-3
        call mrqcof(x,y,ndata,a,ma,alpha,beta,nca,chisq)
        ochisq=chisq
        do j=1,ma
          atry(j)=a(j)
        enddo ! j
      endif
!
!     Alter linearized fitting matrix, by augmenting diagonal elements.
!
      do j=1,ma
        do k=1,ma
          covar(j,k)=alpha(j,k)
        enddo ! k
        covar(j,j)=alpha(j,j)*(1.0d0+alamda)
        da(j)=beta(j)
      enddo ! j
!      
      call gaussj(covar,ma,nca,da,1,1)
!      
!      Once converged, evaluate covariance matrix. 
!      
!      if(alamda.eq.0.0d0)then
!        call covsrt(covar,nca,ma,ia,mfit)
!        
!       Spread out alpha to its full size too.
!        
!        call covsrt(alpha,nca,ma,ia,mfit) 
!        return
!      endif
      j=0
!      
!     Did the trial succeed?
! 
      do l=1,ma
        atry(l)=a(l)+da(l)
      enddo ! l
      
      call mrqcof(x,y,ndata,atry,ma,covar,da,nca,chisq)
      
      if(chisq.lt.ochisq)then     ! Success, accept the new solution.
        alamda=0.1d0*alamda
        ochisq=chisq
        do j=1,ma
          do k=1,ma
            alpha(j,k)=covar(j,k)
          enddo ! k
          beta(j)=da(j)
        enddo ! j
        do l=1,ma
          a(l)=atry(l)
        enddo ! l
      else                        ! Failure, increase alamda and return.
        alamda=1.0d+1*alamda
        chisq=ochisq
      endif
      return
      end subroutine mrqmin
      
!EOC
