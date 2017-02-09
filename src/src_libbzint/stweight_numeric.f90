!BOP
!
! !ROUTINE: edifwt
!
! !INTERFACE:

      subroutine stweight_numeric(iop_omg,de_vert,omg,wt_vert)
!
! !DESCRIPTION:
!
! This subroutine calculates by numerical integration the weight 
! on the whole small tetrahedron 
! in which the $k$ states are fully occupied and $k-q$ states are fully 
! unoccupied. This is for the $sigfreq=3$ case when we consider the 
! imaginary frequency. 
!                                             
! !USES:
      use tetra_internal, only: weighttol, weightwarn,fout,eta_refreq,& 
     &                          vol_small_tetra,n_gauq,x_gauq,w_gauq  

! !INPUT PARAMETERS:
      implicit none
     
      integer, intent(in) :: iop_omg       ! 0 - real frequency 
                                           ! 1 - imaginary frequency  
      real(8), intent(in) :: de_vert(4)    ! difference of the energy 
                                               ! in k-mesh tetrahedron vertices 
                                               ! and k-q mesh tetrahedron vertices.
     

      real(8), intent(in) :: omg              ! the frequency omega to be calculated

! !OUTPUT PARAMETERS:            
      real(8), intent(out) :: wt_vert(4)   ! the weight on the whole tetrahedron.

! !LOCAL VARIABLES:

      integer(4) :: j,k
      integer(4) :: ivert

      real(8)    :: d0,d10,d20,d30
      real(8)    :: x,y,z         !! common variables used by different internal subroutines

      character(20):: sname="stweight_numeric"

      integer:: nmax = 10, nmin = 3 
      real(8):: eps = 1.e-4
      logical:: l_conv

! 
! !SYSTEM ROUTINES:

      intrinsic datan
      intrinsic dlog
      intrinsic dsign
 
! !REVISION HISTORY:
!
! Created 05.11.2004 by XZL.
!

!EOP
!BOC
      d0 = de_vert(1) 
      d10 = de_vert(2)-de_vert(1)
      d20 = de_vert(3)-de_vert(1)
      d30 = de_vert(4)-de_vert(1)
      call gauq_z(1.d0, wt_vert) 

      contains

!
!     This defines the integrand 
!
      function func(x,y,z) result(fv) 
      real(8):: x,y,z,fv(4) 
      real(8):: de,comm 
      
      de = d0 + x*d10 + y*d20 + z*d30
 
      if(iop_omg.eq.0) then 
        comm = (omg-de)/((omg-de)**2+eta_refreq**2)
      else
        comm = -2.d0*de/(omg*omg + de*de)
      endif 

      fv(1) = (1.d0-x-y-z)*comm
      fv(2) = x*comm
      fv(3) = y*comm
      fv(4) = z*comm
      
      end function 

      function f(xx) result(fx) 
      real(8):: xx,fx(4) 
      x = xx
      fx = func(x,y,z)  
      end function

      function g(yy) result(gy) 
      real(8):: yy,gy(4) 
      y = yy
      call gauq_x(1.0-y-z,gy)
      end function  

      function h(zz) result(hz)
      real(8):: zz,hz(4) 
      z = zz 
      call gauq_y(1.0-z,hz)
      end function  

      subroutine gauq_x(a,s)
      real(8),intent(in) :: a
      real(8),intent(out):: s(4) 
      integer :: i
      s=0.d0
      do i=1,n_gauq
        s = s + w_gauq(i)*a*f(a*x_gauq(i))
      enddo 
      end subroutine
 
      subroutine gauq_y(a,s)
      real(8),intent(in) :: a
      real(8),intent(out):: s(4)
      integer :: i
      s=0.d0
      do i=1,n_gauq
        s = s + w_gauq(i)*a*g(a*x_gauq(i))
      enddo
      end subroutine

      subroutine gauq_z(a,s)
      real(8),intent(in) :: a
      real(8),intent(out):: s(4)
      integer :: i
      s=0.d0
      do i=1,n_gauq
        s = s + w_gauq(i)*a*h(a*x_gauq(i))
      enddo
      end subroutine

      end subroutine stweight_numeric
!EOC      
