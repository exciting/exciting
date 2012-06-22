
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: genericprism
! 
! !INTERFACE:
      subroutine genericprism(corners,w,ical,inf)
!     
! !DESCRIPTION:
!
! This subroutine integrates the convolution weight functions inside a
! generic prism. It divides the prism into the corresponding three
! tetrahedra and calls \verb"generictetra" for each of them.
 
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: ical
      real(8), intent(in) :: corners(3,6) ! The coordinates of the 
!                                           six corners of the prism,
!                                           the first three form the 
!                                           triangle at the basis and
!                                           the last three the triangle 
!                                           at the top, so that the edges
!                                           of the prism are (1,4), (2,5)
!                                           and (3,6)
 
! !OUTPUT PARAMETERS:
      integer(4), intent(out) :: inf
      real(8), intent(out) :: w(4) ! the contribution of the prism to 
!                                    the weight at each corner of the 
!                                    containing tetrahedron.      

! !LOCAL VARIABLES:

      integer(4) :: i,itet,inod,info,infl
      
      real(8), dimension(3,4) :: nodtet
      
      real(8), dimension(4) :: internw

! !REVISION HISTORY:
!
! Created 22nd. April 2004 by RGA, last revised Dec.15th, 2004 by XZL

!EOP
!BOC
      inf=0
      infl=0
      w(1:4)=0.0d0
      do itet=0,2
        do inod=1,4
          nodtet(1:3,inod)=corners(1:3,inod+itet)
        enddo
        call generictetra(nodtet,internw,6,info)
        infl=infl+info*(itet+1)
        do i=1,4
          w(i)=w(i)+internw(i)
        enddo
      enddo
      if(infl.ne.0)then
        inf=1
        write(*,'(a7,i4,a8,i4)')'infl = ',infl,' icap = ',ical
        do inod=1,6
          write(*,'(3f13.8)')corners(inod,1:3)
        enddo
      endif    
      end subroutine genericprism

!EOC



