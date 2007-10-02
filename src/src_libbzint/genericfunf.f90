!BOP
!
! !ROUTINE: genericfunf
! 
! !INTERFACE: 

     subroutine genericfunf(corners,w)

     implicit none
!
! !DESCRIPTION:
!This subroutine integrates the convolution weight functions inside a generic pentahedron.
!It divides the pentahedron into the corresponding two tetrahedra and calls \verb"generictetra"
!for each of them

! !INPUT PARAMETERS:

     real(8), intent(in) :: corners(5,3)

! !OUTPUT PARAMETERS:

     real(8), intent(out) :: w(4)

! !LOCAL VARIABLES:

     integer(4) :: i, itet, inod

     real(8), dimension(4) :: internw

     real(8), dimension(4,3) :: nodtet

! !REVISION HISTORY:
!   Created 29nd. August 2004 by XZL; last revised Jan.5th 2005
!EOP
!BOC

     w(1:4)=0.0d0
     do itet=0,1
      do inod=1,4
        nodtet(inod,1:3)=corners(inod+itet,1:3)
      enddo
      call generictetra(nodtet,internw)
      do i=1,4
         w(i)=w(i)+internw(i)
      enddo
     enddo

     end subroutine genericfunf
!EOC

