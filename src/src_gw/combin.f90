!BOP
!
! !ROUTINE: combin
!
! !INTERFACE:
      real(8) function combin(num,denom1,denom2)
      
! !DESCRIPTION:
!
! Calculates the factorial quotients $\frac{n_1!}{n_2!n_3!}$
!      
! !INPUT PARAMETERS:
      
      implicit none
      
      integer(4), intent(in) :: num,denom1,denom2
      
! !LOCAL VARIABLES:      

      integer(4) :: i
      integer(4) :: dmin,dmax
      real(8)    :: cm

 
! !EXTERNAL ROUTINES: 

      real(8), external :: factr

! !INTRINSIC ROUTINES: 

      intrinsic min
      intrinsic max
!
! !REVISION HISTORY:
! 
! Created March  2004 by RGA
! Last modified July 20th 2004 by RGA
!

!EOP
!BOC
!
!     initialize the value
!      
      cm=1.0d+0
!
!     select the larger and smaller denominator
!      
      dmin=min(denom1,denom2)
      dmax=max(denom1,denom2)
!
!     if one of the denominators is 0 or 1 then just call factr
!       
      if(dmin.le.1)then
        cm=factr(num,dmax)
      else
        if(num.ge.dmax)then
          do i=1,dmin
            cm=cm/dble(i)
          enddo
          do i=dmax+1,num
            cm = cm * dble(i)
          enddo
        elseif(num.ge.dmin)then
          do i=1,dmin
            cm=cm/dble(i)
          enddo
          do i=num+1,dmax
            cm = cm / dble(i)
          enddo
        else
          do i=1,num    
            cm=cm/dble(i)
          enddo
          do i=num+1,dmin    
            cm=cm/(dble(i)*dble(i))
          enddo
          do i=dmin+1,dmax
            cm=cm/dble(i)
          enddo
        endif
      endif
      
      combin=cm
      
      end function combin
      
!EOC      
          
            
