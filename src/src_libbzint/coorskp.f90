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

      ! <contribution person="Stephan Sagmeister" date="2007/05/01", topic="arbitrary k-point offset">
!!$      k(3)=mod(id-1,(div(3)))
!!$      k(2)=mod(id-1,(div(3))*(div(2)))/(div(3))
!!$      k(1)=(id-1)/((div(3))*(div(2)))

      k(1)=mod(id-1,(div(1)))
      k(2)=mod(id-1,(div(1))*(div(2)))/(div(1))
      k(3)=(id-1)/((div(1))*(div(2)))
      ! </contribution>

      return
      
      end subroutine coorskp

!EOC      
