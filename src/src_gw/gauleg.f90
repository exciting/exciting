!BOP
!
! !ROUTINE: gauleg
!
! !INTERFACE:
      subroutine gauleg(x1,x2,x,w,n)
      
! !DESCRIPTION:
!
! Given the lower and upper limits of integration \texttt{x1} and
! \texttt{x2}, and given \texttt{n}, this routine returns arrays
! \texttt{x(1:n)} and \texttt{w(1:n)} of length \texttt{n}, containing the
! abscissas and weights of the Gauss-Legendre \texttt{n}-point quadrature
! formula.

! !INPUT PARAMETERS:

      implicit none
      
      
      integer(4), intent(in) :: n  ! Order of the quadrature
      
      real(8),    intent(in) :: x1 ! Lower integration limit
      real(8),    intent(in) :: x2 ! Upper integration limit
      
! !OUTPUT PARAMETERS:

      real(8),    intent(out) :: x(n) ! abscissas      
      real(8),    intent(out) :: w(n) ! weights
      
! !LOCAL VARIABLES:

      integer(4) :: i    
      integer(4) :: j            
      integer(4) :: m
      
      real(8) :: dj
      real(8) :: p1            
      real(8) :: p2            
      real(8) :: p3            
      real(8) :: pp            
      real(8) :: xl            
      real(8) :: xm            
      real(8) :: z            
      real(8) :: z1            
      
      logical(1) :: conv

! !DEFINED PARAMETERS:

      real(8), parameter :: eps = 3.0d-14
      
      real(8), parameter :: pi  = 3.141592653589793d+0
      

!
! !INTRINSIC ROUTINES: 
!

      
      intrinsic abs      
      intrinsic cos      
      intrinsic dble

! !REVISION HISTORY:
!      
! Original subroutine: gauleg.for (C) copr. 1986-92 copr. 1986-92
! numerical recipes software &145i..  
! Last modified: 20.06.05 by RGA.    
!
!EOP
!BOC
!  
!     The roots are symmetric in the interval, so we only have to find
!     half of them.
!
      m = ( n + 1 ) / 2
      xm = 0.5d0 * ( x2 + x1 )
      xl = 0.5d0 * ( x2 - x1 )
      
      do i = 1, m ! Loop over the desired roots
        z = cos(pi*(dble(i)-0.25d0)/(dble(n)+0.5d0))
!       
!       Starting with the above approximation to the \texttt{ith} root, we
!       enter the main loop of refinement by Newton's method 
!
        conv = .false.
        do while (.not.conv)
          p1 = 1.0d0
          p2 = 0.0d0
!        
!         Loop up the recurrence relation to get the Legendre polynomial
!         evaluated at z 
! 
          do j = 1, n
            dj=dble(j)
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*dj-1.0d0)*z*p2 - (dj-1.0d0)*p3)/dj
          enddo ! j  
!      
!         p1 is now the desired Legendre polynomial. We next compute pp, its
!         derivative, by a standard relation involving also p2, the
!         polynomial of one lower order.          
!
          pp = dble(n) * (z * p1 - p2) / (z * z - 1.0d0) 
!
!         Newton's method.
!        
          z1 = z
          z = z1 - p1 / pp
          conv = abs(z-z1).le.eps
        enddo ! conv  
!
!       Scale the root to the desired interval, and put in its symmetric
!       counterpart.
!        
        x(i) = xm - xl * z
        x(n+1-i) = xm + xl * z
!        
!        Compute the weight and its symmetric counterpart.
!        
        w(i) = 2.0d0 * xl / ((1.0d0 - z * z) * pp * pp)
        w(n+1-i) = w(i)
      enddo ! i  
!
!     Successful exit
!        
      return

      end subroutine gauleg
!EOC      
