
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: gbass
!
! !INTERFACE:
      subroutine gbass(rbs,gbs)

! !DESCRIPTION:
!
! Calculates reciprocal lattice vectors from real space
! lattice vectors
!
! !INPUT PARAMETERS:

      implicit none

      real(8), intent(in), dimension(3,3) :: rbs ! The real space 
!                                                  lattice vectors

! !OUTPUT PARAMETERS:

      real(8), intent(out), dimension(3,3) :: gbs ! The reciprocal space
!                                                   lattice vectors

! !REVISION HISTORY:
! Last modified February, 25th. 2004 (RGA)
!
! !LOCAL VARIABLES:
!
      integer(4) :: i
      integer(4) :: j

      real(8) ::  det
      
! !DEFINED PARAMETERS:

      real(8), parameter :: pi = 3.14159265358979323846      

!EOP
!BOC
!
! Arguments:
!
!
! Local variables:
!
!
! Procedure:
!
      gbs(1,1)=rbs(2,2)*rbs(3,3)-rbs(3,2)*rbs(2,3)
      gbs(2,1)=rbs(3,2)*rbs(1,3)-rbs(1,2)*rbs(3,3)
      gbs(3,1)=rbs(1,2)*rbs(2,3)-rbs(2,2)*rbs(1,3)
      gbs(1,2)=rbs(2,3)*rbs(3,1)-rbs(3,3)*rbs(2,1)
      gbs(2,2)=rbs(3,3)*rbs(1,1)-rbs(1,3)*rbs(3,1)
      gbs(3,2)=rbs(1,3)*rbs(2,1)-rbs(2,3)*rbs(1,1)
      gbs(1,3)=rbs(2,1)*rbs(3,2)-rbs(3,1)*rbs(2,2)
      gbs(2,3)=rbs(3,1)*rbs(1,2)-rbs(1,1)*rbs(3,2)
      gbs(3,3)=rbs(1,1)*rbs(2,2)-rbs(2,1)*rbs(1,2)
      det=0.0d0
      do i=1,3
        det=det+gbs(i,1)*rbs(i,1)
      enddo
      do i=1,3
        do j=1,3
          gbs(i,j)=gbs(i,j)*2.0d0*pi/det
        enddo
      enddo

      end subroutine gbass
!EOC
