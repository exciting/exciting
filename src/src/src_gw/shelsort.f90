!BOP
!
! !ROUTINE: shelsort
!
! !INTERFACE:
      subroutine shelsort(n,v,len)


! !DESCRIPTION:
!
! Sorts a set of vectors (v) by increasing length (len) using the shell
! algorithm.

! !INPUT PARAMETERS:

      implicit none     
      integer(4), intent(in) :: n       ! number of vectors to sort      

! !INPUT/OUTPUT PARAMETERS:

      integer(4), intent(inout) :: v(3,n) ! set of vectors
      real(8),    intent(inout) :: len(n) ! set of lengths 

! !LOCAL VARIABLES:

      integer(4) :: gap
      integer(4) :: i
      integer(4) :: vtemp(3)

      real(8) :: ltemp

      logical :: done

! !REVISION HISTORY:
!
! Created Nov 2006 by RGA
!
!EOP
!BOC
      gap=n/2
      do while(gap.ge.1)
        done=.false.
        do while(.not.done)
          done=.true.
          do i=1,n-gap
!            if(len(i).gt.len(i+gap))then
            if(len(i)-len(i+gap).gt.1.d-10)then
              ltemp=len(i)
              vtemp(1:3)=v(1:3,i)
              len(i)=len(i+gap)
              v(1:3,i)=v(1:3,i+gap)
              len(i+gap)=ltemp
              v(1:3,i+gap)=vtemp(1:3)
              done=.false.
            endif
          enddo
        enddo
        gap=gap/2
      enddo
      return
      end subroutine shelsort
!EOC
