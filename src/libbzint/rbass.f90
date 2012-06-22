
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: rbass
!
! !INTERFACE:
      subroutine rbass(gbs,rbs)

! !DESCRIPTION:
!
! Calculates real lattice vectors from reciprocal space
! lattice vectors or vice versa (without $2\pi$ factor)
!
! !INPUT PARAMETERS:

      implicit none
      
      real(8), intent(in), dimension(3,3) :: gbs ! The reciprocal space 
!                                                   lattice vectors

! !OUTPUT PARAMETERS:

      real(8), intent(out), dimension(3,3) :: rbs ! The real space 
!                                                   lattice vectors

! !REVISION HISTORY:
! Last modified February, 25th. 2004 (RGA)
!
! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: j

      real(8) ::  det
      real(8), dimension(3,3) ::  help

!EOP
!
!BOC
!
! Procedure:
!
      rbs(1,1)=gbs(2,2)*gbs(3,3)-gbs(3,2)*gbs(2,3)
      rbs(2,1)=gbs(3,2)*gbs(1,3)-gbs(1,2)*gbs(3,3)
      rbs(3,1)=gbs(1,2)*gbs(2,3)-gbs(2,2)*gbs(1,3)
      rbs(1,2)=gbs(2,3)*gbs(3,1)-gbs(3,3)*gbs(2,1)
      rbs(2,2)=gbs(3,3)*gbs(1,1)-gbs(1,3)*gbs(3,1)
      rbs(3,2)=gbs(1,3)*gbs(2,1)-gbs(2,3)*gbs(1,1)
      rbs(1,3)=gbs(2,1)*gbs(3,2)-gbs(3,1)*gbs(2,2)
      rbs(2,3)=gbs(3,1)*gbs(1,2)-gbs(1,1)*gbs(3,2)
      rbs(3,3)=gbs(1,1)*gbs(2,2)-gbs(2,1)*gbs(1,2)
      det=0.0d0
      do i=1,3
        det=det+rbs(i,1)*gbs(i,1)
      enddo
      do i=1,3
        do j=1,3
          help(i,j)=rbs(i,j)/det
        enddo
      enddo
      do i=1,3
        do j=1,3
          rbs(i,j)=help(j,i)
        enddo
      enddo

      end subroutine rbass

!EOC
