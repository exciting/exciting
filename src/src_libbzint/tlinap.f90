!BOP
!
! !ROUTINE: tlinap

! !INTERFACE:
      real(8) function tlinap(x,f)
      
! !DESCRIPTION:
!      
! This function evalueates a function $f(\vec{x})$  at a given point \verb"x"
! (in internal coordinates) inside the tetrahedron aproximating it linearly 
! from the known values of the function at the corners of the tetrahedron 
! (\verb"f") using the isoparametrization.
!
! \begin{equation}
! \bar{f}(\vec{x})=f_0+\sum\limits_{i=1}^3{(f_i-f_0)x_i}
! \end{equation}

! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in) :: x(1:3) ! The coordinates of the point
      
      real(8), intent(in) :: f(1:4) ! The values of the function at the 
!                                   corners of the tetrahedron

! !LOCAL VARIABLES:

      integer(4) :: i

      real(8) :: fx

      real(8) :: df

! !REVISION HISTORY:
!
! Created 21st. April 2004 by RGA
!
!EOP
!BOC
      
      fx=f(1)
      do i=1,3
        df=f(i+1)-f(1)
        fx=fx+df*x(i)
      enddo
      tlinap=fx
      end function tlinap
!EOC
