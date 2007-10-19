!bop
!
! !FUNCTION: jget
!
! !INTERFACE:
      integer(4) function jget(a,ind)

! !DESCRIPTION:
! 
! Reads the value of the bit ind of vector a
!
! !ARGUMENTS:
      implicit none

      integer(4), intent(in) :: a(*)

      integer(4), intent(in) :: ind
! !RETURN VALUE:
!     
!     integer(4) :: jget            

! !LOCAL VARIABLES:

      integer(4) :: i,off

! !SYSTEM ROUTINES:

      intrinsic ibits
      
!EOP
!BOC
      i=(ind-1)/32+1
      off=mod(ind-1,32)
      jget=ibits(a(i),off,1)
      return
      end function jget
!EOC      

