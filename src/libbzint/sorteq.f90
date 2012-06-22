
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !SUBROUTINE: sorteq
!
! !INTERFACE:
      subroutine sorteq(v,index,sig)
      
! !DESCRIPTION:
!
! This subroutine sorts the values of the vector v (integer) in an
! order that the first one is the one with the biggest possibility of 
! being equal to the others, then the second one, third one and the 
! forth one. Optionally, it also returns an integer vector that gives 
! the old indexes in vector v according to the new ones
!
! !INPUT PARAMETERS:

      implicit none
            
! !INPUT/OUTPUT PARAMETERS:

      real(8), dimension(4), intent(inout) :: v

! !OUTPUT PARAMETERS:      
      integer(4), dimension(4), intent(out) :: index  ! this index tells us 
!                       the relation between the sorted and unsorted arrays.

      integer(4), dimension(5), intent(out) :: sig    ! sig(1:n) means the 
!             number that how many items in the array equals to this specific
!             item. For example, if all of these four are not equal to each
!             other, sig(1:n)=1. 

! !LOCAL VARIABLES:

      integer(4) :: i,j,itmp
      real(8)    :: vtmp,ztol
      
! !SYSTEM ROUTINES:

      intrinsic abs
            
! !REVISION HISTORY:
!
! Created 03rd. Nov. 2004 by XZL
!
!EOP
!BOC
      ztol=1.0d-2
      
      do i=1,5
        sig(i)=0
      enddo

      do i=1,4
        do j=1,4
          if(abs(v(i))+abs(v(j)).gt.1.0d-05)then
          if(2.0d0*abs(v(i)-v(j)).le.(ztol*abs(v(i)+v(j))))sig(i)=sig(i)+1
          else
            sig(i)=sig(i)+1
          endif    
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
          if(2.0d0*abs(v(i)-v(1))/abs(v(i)+v(1)).lt.ztol) then
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
          if((sig(i).eq.2).and.(j.eq.1)) then
            vtmp=v(1)
            v(1)=v(i)
            v(i)=vtmp
            itmp=index(1)
            index(1)=index(i)
            index(i)=itmp
            itmp=sig(1)
            sig(1)=sig(i)
            sig(i)=itmp
            j=j+1
          endif
        enddo
        do i=3,4                ! make the last two with sig(i)=1, by
!                        exchanging it with the second item when sig(i)=2
          if(sig(i).eq.2) then
            vtmp=v(2)
            v(2)=v(i)
            v(i)=vtmp
            itmp=index(2)
            index(2)=index(i)
            index(i)=itmp
            itmp=sig(2)
            sig(2)=sig(i)
            sig(i)=itmp
          endif
        enddo

      case (4)      ! all different, nothing to do
        
        continue
      
      case default   ! None of the above... ERROR
         write(*,*)'sorteq: case not found'
         write(*,*)'sig(i)=', sig
         write(*,*)'v = ',v
         stop
!                    sig(5)=4 and sig(5) corresponds to the all not equal 
!                    and the all equal case respectively
      end select
  
      end subroutine sorteq
!EOC
 
