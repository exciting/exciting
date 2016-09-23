!BOP
!
! !ROUTINE: gaulag
!
! !INTERFACE:
      subroutine gaulag(x,w,n,alfa)
      
! !DESCRIPTION:
!
! Given \texttt{alfa}, the parameter $\alpha$ of the \texttt{Laguerre}
! polynomials, this routine returns arrays \texttt{x(1:n)} and
! \texttt{w(1:n)} containing the abscissas and weights of the
! \texttt{n}-point Gauss-Laguerre quadrature formula. The smallest abscissa
! is returned in \texttt{x(1)}, the largest in \texttt{x(n)}.
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: n    ! Order of the quadrature
      
      real(8),    intent(in) :: alfa ! alpha parameter of the Laguerre
!                                      polynomials

! !OUTPUT PARAMETERS:

      real(8),    intent(out) :: x(n) ! abscissas
      real(8),    intent(out) :: w(n) ! weights
      
! !LOCAL VARIABLES:
      
      integer(4) :: i      
      integer(4) :: its
      integer(4) :: j
      
      real(8) :: ai
      real(8) :: dj
      real(8) :: p1      
      real(8) :: p2      
      real(8) :: p3      
      real(8) :: pp      
      real(8) :: z      
      real(8) :: z1     
      
      logical(1) :: conv 
      
! !DEFINED PARAMETERS:

      integer(4), parameter :: maxit = 10
      
      real(8),    parameter :: eps = 3.0d-14
            

!
! !EXTERNAL ROUTINES: 
!


      real(8), external :: gammln


!
! !INTRINSIC ROUTINES: 


      intrinsic dble
      intrinsic exp
        
! !REVISION HISTORY:
!      
! Original subroutine: gauleg.for (C) copr. 1986-92 copr. 1986-92
! numerical recipes software &146i..  
! Last modified: 20.06.05 by RGA.    
!
!EOP
!BOC
!     
!     Loop over the desired roots.
!
      do i = 1, n
        if(i.eq.1)then              ! Initial guess for the smallest root.

          z = (1.0d0 + alfa) * (3.0d0+ 0.92d0 * alfa) / (1.0d0 +        &
     &         2.4d0 * dble(n) + 1.8 * alfa)

        else if(i.eq.2)then         ! Initial guess for the second root.

          z = z+(1.5d+1+6.25d0*alfa)/(1.0d0+0.9d0*alfa+2.5d0*dble(n))

        else                        ! Initial guess for the other roots.  
        
          ai = dble(i-2)
          z= z + ((1.0d0+2.55d0*ai)/(1.9d0*ai)+1.26d0*ai*alfa/          &
     &            (1.0d0+3.5d0*ai))*(z-x(i-2))/(1.0d0+0.3d0*alfa)
        
        endif     
!       
!       Starting with the above approximation to the \texttt{ith} root, we
!       enter the main loop of refinement by Newton's method 
!
        conv = .false.
        its = 0
        do while (.not.conv)
          its = its + 1
          if(its.gt.10)goto 999
          p1 = 1.0d0
          p2 = 0.0d0
!          
!         Loop up the recurrence relation to get the Laguerre polynomial
!         evaluated at z
!
          do j = 1, n
            dj = dble(j)
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*dj-1.0d0+alfa-z)*p2-(dj-1.0d0+alfa)*p3)/dj
          enddo ! j        
!          
!         p1 is now the desired Laguerre polynomial. We next compute pp,
!         its derivative, by a standard relation involving also p2, the
!         polynomial of one lower order.
!
          pp = (dble(n)*p1-(dble(n)+alfa)*p2)/z          
!
!         Newton's method.
!        
          z1 = z
          z = z1 - p1 / pp
          conv = abs(z-z1).le.eps
        enddo ! conv     
!
!       Store the root and the weight
!        
        x(i) = z
        w(i) = -exp(gammln(alfa+dble(n))-gammln(dble(n)))/(pp*dble(n)*p2)
      enddo ! i
!
!     Successful exit
!        
      return
!
!     Error handling
!
  999 stop 'too many iterations in gaulag'

      end subroutine gaulag
!EOC         
      
