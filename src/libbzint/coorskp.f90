
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: coorskp
!
! !INTERFACE:
      subroutine coorskp(id,k)

! !DESCRIPTION:
!
! Given the identification number of the k-point it returns its
! coordinates in the submesh

! !USES:

      use kgen_internals, only: div      
!<sag>
      use control, only: tetraifc
!</sag>
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: id     ! Id. number of the k-point
      
! !OUTPUT PARAMETERS:

      integer(4), intent(out) :: k(3)  ! Coordinates of the k-point
      
!
! !REVISION HISTORY:
! 
! Created 24th. Feb. 2004
!
!EOP
!
!BOC

!<sag>
      if (trim(tetraifc)=='wien2k') then
         ! original code
         k(3)=mod(id-1,(div(3)))
         k(2)=mod(id-1,(div(3))*(div(2)))/(div(3))
         k(1)=(id-1)/((div(3))*(div(2)))
      else if (trim(tetraifc)=='exciting') then
         ! new code
         k(1)=mod(id-1,(div(1)))
         k(2)=mod(id-1,(div(1))*(div(2)))/(div(1))
         k(3)=(id-1)/((div(1))*(div(2)))
      end if
!</sag>

      return
      
      end subroutine coorskp

!EOC      
