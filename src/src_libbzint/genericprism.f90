!BOP
!
! !ROUTINE: genericprism
! 
! !INTERFACE:
      subroutine genericprism(corners,w)
!     
! !DESCRIPTION:
!
! This subroutine integrates the convolution weight functions inside a
! generic prism. It divides the prism into the corresponding three
! tetrahedra and calls \verb"generictetra" for each of them.
 
! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in) :: corners(6,3) ! The coordinates of the 
!                                           six corners of the prism,
!                                           the first three form the 
!                                           triangle at the basis and
!                                           the last three the triangle 
!                                           at the top, so that the edges
!                                           of the prism are (1,4), (2,5)
!                                           and (3,6)
 
! !OUTPUT PARAMETERS:

      real(8), intent(out) :: w(4) ! the contribution of the prism to 
!                                    the weight at each corner of the 
!                                    containing tetrahedron.      

! !LOCAL VARIABLES:

      integer(4) :: i,itet,inod
      
      real(8), dimension(4,3) :: nodtet
      
      real(8), dimension(4) :: internw

! !REVISION HISTORY:
!
! Created 22nd. April 2004 by RGA, last revised Dec.15th, 2004 by XZL

!EOP
!BOC

      w(1:4)=0.0d0
      do itet=0,2
        do inod=1,4
          nodtet(inod,1:3)=corners(inod+itet,1:3)
        enddo
        call generictetra(nodtet,internw)
        do i=1,4
          w(i)=w(i)+internw(i)
        enddo
      enddo

      end subroutine genericprism

!EOC



