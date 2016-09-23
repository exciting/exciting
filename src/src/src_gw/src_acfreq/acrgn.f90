!BOP
!
! !ROUTINE: acrgn
!
! !INTERFACE:
      subroutine acrgn(npar,x,ca,y,dy)
 
! !DESCRIPTION:
!
! Given the abscissa \texttt{x}, and the parameters \texttt{a(1:2n+2)}
!this subroutine returns in \texttt{y} the value of the function:
!
!\begin{equation}\label{acrgn-1}     
!f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
!\end{equation}
!
!where 

!\begin{equation}\label{acrgn-2}     
!P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
!\end{equation}
!
!\begin{equation}\label{acrgn-3}     
!Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
!\end{equation}
!
!and $c_k=\texttt{a(k)}+i\texttt{a(k+2n+2)}$. In \texttt{dy} the 
!derivative of the function with respect to \texttt{x}:
!
!\begin{equation}\label{acrgn-4}     
!\frac{\partial f(x,\{c\})}{\partial x}=\frac{1}{Q(x,\{c\})}\frac{\partial %
!P(x,\{c\})}{\partial x}-\frac{P(x,\{c\})}{[Q(x,\{c\})]^2}\frac{\partial %
!Q(x,\{c\})}{\partial x} %
!\end{equation}
! 
! with:
!
!\begin{subequations}\label{acrgn-5}
!\begin{align}
!\frac{\partial P(x,\{c\})}{\partial x}=&\sum\limits_{k=2}^{n+1}{%
!(k-1)c_{k}x^{k-2}}\\
!\frac{\partial Q(x,\{c\})}{\partial x}=&
!\sum\limits_{k=n+2}^{2n+2}{(k-n-1)c_k x^{k-n-2}}\\
!\end{align}
!\end{subequations}
!

! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: npar        ! number of fitting parameters 
      complex(8), intent(in) :: x,ca(npar)    ! fitting coeffients 
      
! !OUTPUT PARAMETERS:
      complex(8),    intent(out) :: y       ! fitted selfec at real frequency 
      complex(8),    intent(out) :: dy      ! first order derivative of selfec with respect to frequency
      
! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: n
      
      complex(8) :: q
      complex(8) :: p
      complex(8) :: dq
      complex(8) :: dp
      
     
! !DEFINED PARAMETERS:
      
      complex(8), parameter :: czero = (0.0d0,0.0d0)      
      complex(8), parameter :: cone = (1.0d0,0.0d0)      
      complex(8), parameter :: cim = (0.0d0,1.0d0)      
      
! !REVISION HISTORY:
!
! Created: 13th. Jul. 2005 by RGA
!
!EOP
!BOC 
!
!     calculate the qominator and perator of y
!
      n=(npar-2)/2
      q=czero
      p=czero
      dp=czero
      dq=czero
      do i=n,1,-1
        q = (q + ca(n+2+i)) * x
        p = (p + ca(i+1)) * x
      enddo
      q = (q + ca(n+2))* x + cone
      p = p + ca(1)
      do i=n,2,-1
        dq = (dq + ca(n+2+i)*dble(i+1))*x
        dp = (dp + ca(i+1)*dble(i)) * x
      enddo  
      dq = dq + ca(n+2) + 2.0d0*ca(n+3)*x
      dp = dp + ca(2)
      y = p / q
      dy = dp / q - y * dq / q
      return

      end  subroutine acrgn
!EOC      
      
        
      
