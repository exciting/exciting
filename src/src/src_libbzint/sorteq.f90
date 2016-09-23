!BOP
!
! !SUBROUTINE: sorteq
!
! !INTERFACE:
      subroutine sorteq(v,index,sigeq)
      
! !DESCRIPTION:
!
! This subroutine sorts the values of the vector v (integer) in an
! order that the first one is the one with the biggest possibility of 
! being equal to the others, then the second one, third one and the 
! forth one. Optionally, it also returns an integer vector that gives 
! the old indexes in vector v according to the new ones
!
! !INPUT PARAMETERS:
      use tetra_internal, only: ztol_sorteq
      implicit none
            
! !INPUT/OUTPUT PARAMETERS:

      real(8), dimension(4), intent(inout) :: v

! !OUTPUT PARAMETERS:      
      integer(4), dimension(4), intent(out) :: index  ! this index tells us 
!                       the relation between the sorted and unsorted arrays.

      integer(4), intent(out) :: sigeq    

! !LOCAL VARIABLES:

      integer(4) :: i,j,itmp
      integer(4) :: sig(5)   ! sig(1:4) number of items in the array equal to
                             !  this specific item. For example, if all are 
                             !  different sig(1:4)=1. sig(5) is the sum 
      real(8)    :: vtmp,dist
      
! !SYSTEM ROUTINES:

      intrinsic abs
            
! !REVISION HISTORY:
!
! Created 03rd. Nov. 2004 by XZL
!
!EOP
!BOC
      

      sig=0 
      do i=1,4
        do j=1,4
          if(f_rdist(v(i),v(j)).lt.ztol_sorteq) sig(i)=sig(i)+1
        enddo
      enddo
      
      sig(5)=sig(1)+sig(2)+sig(3)+sig(4)
      
      if((sig(5).eq.8).and.(maxval(sig(1:4)).eq.3))then
        
        do i=1,4
          if(sig(i).eq.2)sig(i)=sig(i)+1
        enddo
        sig(5)=sig(1)+sig(2)+sig(3)+sig(4)
      endif  

      if((sig(5).eq.10).and.(minval(sig(1:4)).eq.2))then
        do i=1,4
          sig(i)=4
        enddo
        sig(5)=sig(1)+sig(2)+sig(3)+sig(4)
      endif  
          
      if(sig(5).eq.12)sig(5)=16 
      if(sig(5).eq.14)sig(5)=16 
      do i=1,4
        index(i)=i
      enddo

      select case (sig(5))
      case (16)

        continue

      case (10)                  ! three of them are equal

        do i=1,3                  ! we want to make it as a=b=c while d not
          if(sig(i).eq.1) then
            vtmp=v(i)                            ! we make sig(4)=1
            v(i)=v(4)
            v(4)=vtmp
            itmp=index(i)
            index(i)=index(4)
            index(4)=itmp
            itmp=sig(4)
            sig(4)=sig(i)
            sig(i)=itmp
          endif
        enddo
    
      case (8)                   ! doubly equal but not all

        do i=3,4                                 ! make the first two equal
          if( f_rdist(v(i),v(1)).lt.f_rdist(v(2),v(1)) )then
            vtmp=v(i)
            v(i)=v(2)
            v(2)=vtmp
            itmp=index(i)
            index(i)=index(2)
            index(2)=itmp
            itmp=sig(i)
            sig(i)=sig(2)
            sig(2)=itmp
          endif
        enddo

      case (6)                   ! only two of them are equal

        j=1
        do i=1,4                         ! make the first one with sig(1)=2
          if((sig(i).eq.2)) then
            vtmp=v(j)
            v(j)=v(i)
            v(i)=vtmp
            itmp=index(j)
            index(j)=index(i)
            index(i)=itmp
            itmp=sig(j)
            sig(j)=sig(i)
            sig(i)=itmp
            j=j+1
           endif
        enddo

      case (4)      ! all different, nothing to do
        
        continue
      
      case default   ! None of the above... ERROR
         write(6,*)'ERROR in sorteq: case not found'
         write(6,*)'sig(i)=', sig
         write(6,'(a,4g16.6)')'v = ',v
         stop "ERROR in sorteq"
!                    sig(5)=4 and sig(5) corresponds to the all not equal 
!                    and the all equal case respectively
      end select
      sigeq = sig(5)

      contains 
!
!     this internal function measures the relative distance
!
      real(8) function f_rdist(v1,v2) 
      real(8),intent(in) :: v1,v2
      
      if((abs(v1)+abs(v2)).lt.1.0e-10) then 
        f_rdist = 0.d0 
      else 
        f_rdist = 2.d0*abs(v1-v2)/(abs(v1)+abs(v2))
      endif 
      end function 
      
      end subroutine sorteq
!EOC
 
