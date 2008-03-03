!BOP
!
! !ROUTINE: idkp
!
! !INTERFACE:
      integer(4) function idkp(k)
       
! !DESCRIPTION:
!
! Given the coordinates in the submesh of the k-point it returns its 
! identification number, which is the opposite process of 'coorskp'.

!      
! !USES:
      use kgen_internals, only: div      
!<sag>
      use control, only: tetraifc
!</sag>
      
! !INPUT PARAMETERS:

      implicit none
      
      integer(4), intent(in) :: k(3)  ! Coordinates of the k-point
!

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
         idkp=k(1)*div(2)*div(3)+k(2)*div(3)+k(3)+1
      else if (trim(tetraifc)=='exciting') then
         ! new code
         idkp=k(3)*div(2)*div(1)+k(2)*div(1)+k(1)+1
      end if
!</sag>      
      return
      
      end function idkp
!EOC
