
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: genericfunf
! 
! !INTERFACE: 

     subroutine genericfunf(corners,w,ical,info)

     implicit none
!
! !DESCRIPTION:
!This subroutine integrates the convolution weight functions inside a generic pentahedron.
!It divides the pentahedron into the corresponding two tetrahedra and calls \verb"generictetra"
!for each of them

! !INPUT PARAMETERS:

     real(8), intent(in) :: corners(3,5)
     integer(4), intent(in) :: ical

! !OUTPUT PARAMETERS:

     real(8), intent(out) :: w(4)
     integer(4), intent(out) ::info

! !LOCAL VARIABLES:

     integer(4) :: i, itet, inod,inf,insp

     real(8), dimension(4) :: internw

     real(8), dimension(3,4) :: nodtet

! !REVISION HISTORY:
!   Created 29nd. August 2004 by XZL; last revised Jan.5th 2005
!EOP
!BOC
     insp=0
     info=0
     w(1:4)=0.0d0
     do itet=0,1
      do inod=1,4
        nodtet(1:3,inod)=corners(1:3,inod+itet)
      enddo
      call generictetra(nodtet,internw,5,inf)
      if(inf.ne.0)then
        insp=insp+inf*(itet+1)
      endif  
      do i=1,4
         w(i)=w(i)+internw(i)
      enddo
     enddo
     if(insp.ne.0)then
        info=1
        write(*,'(a7,i4,a8,i4)')'insp = ',insp,' icap = ',ical
        do inod=1,5
          write(*,'(3f13.8)')corners(inod,1:3)
        enddo
      endif    
  
     end subroutine genericfunf
!EOC

