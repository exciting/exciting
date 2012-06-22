      subroutine init_c(x,y,c,m)

! !DESCRIPTION:
!
! This subroutine calculates the values of the \texttt{2n+2} parameters $c_k$ 
! of a function of the form:
!
!\begin{equation}\label{init_c-1}     
!f(x,\{c\})=\frac{P(x,\{c\})}{Q(x,\{c\})}
!\end{equation}
!
!where 
!
!\begin{equation}\label{init_c-2}     
!P(x,\{c\})=\sum\limits_{k=1}^{n+1}{c_{k}x^{k-1}},
!\end{equation}
!
!\begin{equation}\label{init_c-3}     
!Q(x,\{c\})=1+\sum\limits_{k=n+2}^{2n+2}{c_k x^{k-n-1}}
!\end{equation}
!
!by adjusting them to a set of \texttt{2n+2} $(x,y)$-pairs.
!
!
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: m
      
      complex(8),    intent(in) :: x(*)
      complex(8),    intent(in) :: y(*)
      
! !OUTPUT PARAMETERS:
      
      complex(8),    intent(out) :: c(*)
      
! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: j
      integer(4) :: n
      integer(4) :: ip(1:m)
      integer(4) :: inf
      
      complex(8) :: xtn
      complex(8) :: amat(1:m,1:m)


! !DEFINED PARAMETERS:
      
      complex(8), parameter :: czero = (0.0d0,0.0d0)      
      complex(8), parameter :: cone = (1.0d0,0.0d0)      
      complex(8), parameter :: cim = (0.0d0,1.0d0)      

       
!
! !EXTERNAL ROUTINES: 
!


      external zgetrf
      external zgetrs

! !REVISION HISTORY:
!
! Created: 13th. Jul. 2005 by RGA
!
!EOP
!BOC 
!
!     Initializations
!
      n=m/2-1
      c(1:m)=y(1:m)
!
!     Generation of the matrix
!
      do i=1, m
        xtn=cone
        do j=1, n+1
          amat(i,j) = xtn
          xtn = xtn * x(i)
          amat(i,j+n+1)=-y(i)*xtn
        enddo ! j
      enddo ! i    
      call zgetrf(m,m,amat,m,ip,inf)
      if(inf.ne.0)then
        write(6,*)'acont: zgetrf, inf =',inf
        stop
      endif    
      call zgetrs('n',m,1,amat,m,ip,c,m,inf)
      if(inf.ne.0)then
        write(6,*)'acont: zgetrs, inf =',inf
        stop
      endif    
      
      return
      
      end subroutine init_c

!EOC      
