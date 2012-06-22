!BOP
!
! !ROUTINE: ratfun
!
! !INTERFACE:
      subroutine ratfun(x,a,y,dyda,hess,n)
 
! !DESCRIPTION:
!
! Given the abscissa \texttt{x}, and the parameters \texttt{a(1:4n+4)}
!this subroutine returns in \texttt{y} the value of the function:
!
!\begin{equation}\label{ratfun-1}     
!f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
!\end{equation}
!
!where 

!\begin{equation}\label{ratfun-2}     
!P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
!\end{equation}
!
!\begin{equation}\label{ratfun-3}     
!Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
!\end{equation}
!
!and $c_k=\texttt{a(k)}+i\texttt{a(k+2n+2)}$. In \texttt{dyda} the partial
!derivatives of the function with respect to the parameters:
!
!\begin{equation}\label{ratfun-4}     
!\frac{\partial f(x,\{c\})}{\partial c_j}=\left\{
!\begin{array}{lc}
!{\displaystyle \frac{x^{j-1}}{Q(x,\{c\})}}& 1 \leq j \leq n+1 \vspace{3mm}\\
!{\displaystyle -\frac{P(x,\{c\})}{[Q(x,\{c\})]^2}x^{j-n-1}}& n+2 %
!\leq j \leq 2n+2 \\
!\end{array}
!\right.
!\end{equation}
! 
! The algorithm for calculating the derivatives is:
!
!\begin{subequations}\label{ratfun-5}
!\begin{align}
!\frac{\partial f(x,\{c\})}{\partial c_1}=&
!\frac{1}{Q(x,\{c\})}\\
!\frac{\partial f(x,\{c\})}{\partial c_{n+2}}=&
!-\frac{f(x,\{c\})}{Q(x,\{c\})}\\
!\frac{\partial f(x,\{c\})}{\partial c_j}=&\left\{
!\begin{array}{lc}
!{\displaystyle \frac{\partial f(x,\{c\})}{\partial c_{j-1}}x}& 2 \leq j \leq n+1 \\
!{\displaystyle -f(x,\{c\})\frac{\partial f(x,\{c\})}{\partial%
!c_{j-n-1}}}& n+3 \leq j \leq 2n+2 \\
!\end{array}
!\right.
!\end{align}
!\end{subequations}
!
!
! The hessian is returned in the array \texttt{hess}, and is calculated
!as:
!
!\begin{equation}\label{ratfun-6}     
!\frac{\partial^2 f(x,\{c\})}{\partial c_k\partial c_j}=\left\{
!\begin{array}{lc}
! 0 & 1 \leq j,k \leq n+1 \\
!{\displaystyle -\frac{x^{j-1+k-n-1}}{[Q(x,\{c\})]^2}}& %
!1 \leq j \leq n+1 < k \leq 2n+2 \vspace{3mm}\\
!{\displaystyle 2\frac{P(x,\{c\})}{[Q(x,\{c\})]^3}x^{j+k-2n-2}}& n+2 %
!\leq j,k \leq 2n+2 \\
!\end{array}
!\right.
!\end{equation}
!
! which are calculated as:
!
!\begin{equation}\label{ratfun-7}     
!\frac{\partial^2 f(x,\{c\})}{\partial c_k\partial c_j}=\left\{
!\begin{array}{lc}
! 0 & 1 \leq j,k \leq n+1 \vspace{3mm}\\
!{\displaystyle -\frac{\partial f(x,\{c\})}{\partial c_{j}}\frac{\partial%
!f(x,\{c\})}{\partial c_{k-n-1}}}&1 \leq j \leq n+1 < k \leq 2n+2 
!\vspace{3mm}\\
!{\displaystyle 2f(x,\{c\})\frac{\partial f(x,\{c\})}{\partial%
!c_{j-n-1}}\frac{\partial f(x,\{c\})}{\partial%
!c_{k-n-1}}}& n+2 %
!\leq j,k \leq 2n+2 \\
!\end{array}
!\right.
!\end{equation}
!
! The derivatives with respect to the real parameters are easily
!calculated by using
!
!\begin{equation}\label{ratfun-8}     
!\frac{\partial c_k}{\partial a_j}=\delta_{k,j}+i\delta_{k,j-2n-2}
!\end{equation}
! and the chain rule.
!\begin{equation}\label{ratfun-9}     
!\frac{\partial f(x,\{c\})}{\partial a_j}=\sum_{k=1}^{2*n+2}%
!{\frac{\partial f(x,\{c\})}{\partial c_k}\frac{\partial c_k}{\partial
!a_j}}
!\end{equation}
!

! !INPUT PARAMETERS:

      implicit none
      
      
      integer(4), intent(in) :: n
      
      complex(8),    intent(in) :: x
      real(8),    intent(in) :: a(1:4*n+4)
      
! !OUTPUT PARAMETERS:
      
      complex(8),    intent(out) :: dyda(1:4*n+4)
      complex(8),    intent(out) :: hess(1:4*n+4,1:4*n+4)
      complex(8),    intent(out) :: y
      
! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: j
      integer(4) :: l
      integer(4) :: m
      
      complex(8) :: q
      complex(8) :: p
      
      complex(8) :: ca(1:2*n+2)
      
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
!     Initializations
!
      dyda(:)=czero
!
!     Convert the real coefficient to complex ones.      
!
      do i=1,2*n+2
        ca(i)=cmplx(a(i),a(2*n+2+i),8)
      enddo  
      
!
!     calculate the qominator and perator of y
!
      q=czero
      p=czero
      do i=1,n
        q = (q + ca(2*n+3-i)) * x
        p = (p + ca(n+2-i)) * x
      enddo
      q = (q + ca(n+2))* x + cone
      p = p + ca(1)
      y = p / q
!
!     Calculate the partial derivatives with respect to the complex
!     coefficients ca(:) (equal to the derivative with respect to the real
!     part of the coefficient
!      
      dyda(1) = 1 / q
!      dyda(n+2) = - y * dyda(1)
      do i=1, n
        dyda(i+1) = dyda(i) * x
        dyda(n+1+i) = - y * dyda(i+1)
      enddo  
      dyda(2*n+2) = - y * dyda(n+1) * x
!
!     Convert to the derivatives with respect to the real coeficients a(:)
!     For the real part of the coefficient it is done, just extend to the
!     imaginary part.  
!      
      do i=1,2*n+2
        dyda(2*n+2+i) = dyda(i) * cim
      enddo  
!      
!     calculate the Hessian
!
      do i=1, n+1
        l= i + n + 1
        do j=1, n+1
          m = j + n + 1
          hess(i,j) = czero
          hess(l,j) = -dyda(i) * dyda(j)
          hess(i,m) = hess(l,j)
          hess(l,m) = -y * dyda(i) * dyda(j)
        enddo
      enddo    
!
!     extend to the imaginary part of the coefficients
!        
      do i=1,2*n+2
        do j=1,2*n+2
          hess(2*n+2+i,j) = hess (i,j) * cim      
          hess(i,2*n+2+j) = hess (i,j) * cim      
          hess(2*n+2+i,2*n+2+j)= - hess(i,j)
        enddo
      enddo    
          
      return

      end  subroutine ratfun
!EOC      
      
        
      
