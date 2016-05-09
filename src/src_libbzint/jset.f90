!BOP
!
! !ROUTINE: jset
!
! !INTERFACE:
      subroutine jset(a,ind,val)
      
!
! !DESCRIPTION:
!
! Sets the bit number ind of the vector a to the value val
!
! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: val
      integer(4), intent(in) :: ind
      
! !INPUT/OUTPUT PARAMETERS:
      integer(4), intent(inout) :: a(*)

! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: off
      
! !SYSTEM ROUTINES:

      intrinsic mod
      intrinsic ibclr
      intrinsic ibset      

!EOP
!
!BOC
    
      i=(ind-1)/32+1
      off=mod(ind-1,32)
      if(val.eq.0)then
        a(i)=ibclr(a(i),off)
      else
        a(i)=ibset(a(i),off)
      endif
      return
      end subroutine jset

!EOC
