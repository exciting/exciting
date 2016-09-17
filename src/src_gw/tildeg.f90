!BOP
!
! !ROUTINE: tildeg
!
! !INTERFACE:
      real(8) function tildeg(l1,l2,m1,m2)

!
! !DESCRIPTION: 
!
!This subroutine calculates the coefficients 
! $\tilde{g}_{lm,l'm'}$ according to equation \ref{tildea}. The factor
!$(4\pi)^{\frac{3}{2}}$ is included into $\tilde{g}_{lm,l'm'}$.
!
! !INPUT PARAMETERS:
!
      implicit none

      integer(4) :: l1,l2,m1,m2
!
! !LOCAL VARIABLES:
!
      integer(4), dimension(4) :: jm        ! value of l1+m1
      integer(4) :: tjpo    ! value of 2(l1+l2)+1
      real(8) :: combjpm    ! value of (l1+l2+m1+m2)!/(l1+m1)!(l2+m2)!
      real(8) :: combjmm    ! value of (l1+l2-m1-m2)!/(l1-m1)!(l2-m2)!
      real(8) :: denom
      real(8), parameter :: pi = 3.1415926536d0
      
! !EXTERNAL ROUTINES: 

      real(8), external :: factr 
      real(8), external :: combin

!EOP
!BOC
      jm(1)=l1+m1
      jm(2)=l2+m2
      jm(3)=l1-m1
      jm(4)=l2-m2
      tjpo=2*(l1+l2)+1
      combjpm=combin(jm(1)+jm(2),jm(1),jm(2))
      combjmm=combin(jm(3)+jm(4),jm(3),jm(4))
      denom=dble((2*l1+1)*(2*l2+1)*tjpo)
      tildeg = ((-1.0d0)**l1)*8.0d0*pi*sqrt(pi*combjmm*combjpm/denom)
      end function tildeg
!EOC
