!BOP
!
! !ROUTINE: gammln
!
! !INTERFACE:
      real(8) function gammln(xx)
      
! !DESCRIPTION:
!      
! Returns the value $\ln[\Gamma(\texttt{xx})]$ for $\texttt{xx}>0$.
!
! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in) :: xx

! !LOCAL VARIABLES:

      integer(4) :: j
      
      real(8) :: ser            
      real(8) :: tmp            
      real(8) :: x            
      real(8) :: y

! !DEFINED PARAMETERS: 

      real(8), save :: stp            
      real(8), save :: cof(6)            
      
      data cof /76.18009172947146d0,-86.50532032941677d0,               &
     &          24.01409824083091d0,-1.231739572450155d0,               &
     &          .1208650973866179d-2,-.5395239384953d-5/
      data stp /2.5066282746310005d0/
      

!
! !INTRINSIC ROUTINES: 
!

      
      intrinsic log

! !REVISION HISTORY:
!      
! Original subroutine: gauleg.for (C) copr. 1986-92 copr. 1986-92
! numerical recipes software pp 207..  
! Last modified: 20.06.05 by RGA.    
!
!EOP
!BOC
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      enddo  
      gammln=tmp+log(stp*ser/x)
      return
      end function gammln
!EOC
