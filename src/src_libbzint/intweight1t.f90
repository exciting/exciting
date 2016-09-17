!BOP
!
! !ROUTINE: intweight1t
!
! !INTERFACE:
       subroutine intweight1t(e,ef,v,w)
!       
! !DESCRIPTION:
! 
!   This subroutine calculates the contribution to the integration weights at
!   the four corners of a tetrahedron. The energies have to be given in
!   increasing energy order.     
!       
! !USES:
             
       implicit none

! !INPUT PARAMETERS:
 
       real(8), intent(in) :: e(4) ! Eigenenergies at the corners of the
!                                    tetrahedron.
 
       real(8), intent(in) :: ef   ! Fermi energy.
       
       real(8), intent(in) :: v    ! Volume of the tetrahedron.

! !OUTPUT PARAMETERS:
       
       real(8), dimension(4) :: w ! Weight at each corner.
  
!   
! !REVISION HISTORY:
!
!   Created: 3th. March 2004, by RGA
!
! !LOCAL VARIABLES:
 
       
       integer(4) :: index
       
       integer(4) :: i
       
       real(8)    :: c1,c2,c3,vo4
       
       real(8)    :: f31,f41,f32,f42
 
       real(8)    :: e31,e41,e32,e42
 
 ! !SYSTEM ROUTINES:
 

!EOP
!BOC
      if(e(4).le.ef)then
        index=4
      else
        index=0
        do i=1,4
          if(e(i).le.ef)index=index+1
        enddo  
      endif
      vo4=v/4.0d0
      select case(index)
        case(0)          ! all states are unoccupied
          w(1:4) = 0.0d0
          
        case(1)          ! only the lowest energy is occupied.
          do i=2,4
            w(i) = (ef-e(1))/(e(i)-e(1))
          enddo
          c1 = vo4 * w(2) * w(3) * w(4)
          w(1) = 4.0d0 * c1
          do i=2,4
            w(i) = c1 * w(i) 
            w(1) = w(1) - w(i)
          enddo
          
        case(2)         ! the two lower energies are occupied
          f31=(ef-e(1))/(e(3)-e(1))
          f41=(ef-e(1))/(e(4)-e(1))
          f32=(ef-e(2))/(e(3)-e(2))
          f42=(ef-e(2))/(e(4)-e(2))
          e31=(e(3)-ef)/(e(3)-e(1))
          e41=(e(4)-ef)/(e(4)-e(1))
          e32=(e(3)-ef)/(e(3)-e(2))
          e42=(e(4)-ef)/(e(4)-e(2))
          c1 = vo4 * f41 * f31
          c2 = vo4 * f41 * f32 * e31
          c3 = vo4 * f42 * f32 * e41
          w(1) = c1 + (c1+c2)*e31 + (c1+c2+c3)*e41
          w(2) = c1+c2+c3 + (c2+c3)*e32 + c3*e42
          w(3) = (c1+c2)*f31 + (c2+c3)*f32
          w(4) = (c1+c2+c3)*f41 + c3*f31
          
        case(3)        ! Only the highest energy is unoccupied
          vo4=v/4.0d0
          do i=1,3
            w(i) = (e(4)-ef)/(e(4)-e(i))
          enddo
          c1 = vO4 * w(1) * w(2) * w(3) 
          w(4) = vo4 - 4.0d0 * c1
          do i=1,3
            w(4) = w(4) + c1 * w(i)
            w(i) = vo4 - c1 * w(i) 
          enddo

        case(4)        ! All states are occupied
          do i=1,4
            w(i) = vo4
          enddo
        
      end select

      end subroutine intweight1t

!EOC
