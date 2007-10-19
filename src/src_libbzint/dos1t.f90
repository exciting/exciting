!BOP
!
! !ROUTINE: dos1t 
!
! !INTERFACE:
       real(8) function dos1t(e,ei,v)
!
! !DESCRIPTION:
!  
!   This subroutine calculates the contribution to the density of states at
!   an energy $\epsilon = $ \verb"ei" of one tetrahedron, according to
!   Bl\"ochl {\it et al}, Phys. Rev. B, {\bf 49}, 16223 (1994). The energies at the 
!   corner of the tetrahedron have to be given in increasing order. 
!       
! !USES:
             
       implicit none
!
! !INPUT PARAMETERS:
 
       real(8), intent(in) :: e(4) ! Eigenenergies at the corners of the
!                                    tetrahedron.
 
       real(8), intent(in) :: ei   ! energy at which the contribution to
!                                     the dos1t is calculated
       
       real(8), intent(in) :: v    ! Volume of the tetrahedron.
       
      
!   
! !REVISION HISTORY:
!
!   Created: 3th. March 2004, by RGA
!
! !LOCAL VARIABLES:
 
       integer(4) :: index,i
       real(8)    :: br,denom,eij
       real(8)    :: e21,e31,e41,e32,e42,e43
 
 ! !SYSTEM ROUTINES:
 
 
 
 
!EOP
!BOC 

      if(e(4).le.ei)then
        index=4
      else
        index=0
        do i=1,4
          if(e(i).le.ei)index=index+1
        enddo  
      endif
      
      select case(index)
        case(0,4)          ! all states are either below or above ei
          dos1t = 0.0d0
          
        case(1)          ! only the lowest energy is occupied.
          e21=e(2)-e(1)
          e31=e(3)-e(1)
          e41=e(4)-e(1)
          eij=ei-e(1)
          denom=e21*e31*e41
          dos1t = 3.0d0*v*eij*eij/denom
        
        case(2)         ! e2 < ei < e3
          e21=e(2)-e(1)
          e31=e(3)-e(1)
          e41=e(4)-e(1)
          e32=e(3)-e(2)
          e42=e(4)-e(2)
          eij=ei-e(2)
          denom=e32*e42
          br = e21+2.0d0*eij-(e31+e42)*eij*eij/denom
          denom = e31*e41
          dos1t = 3.0d0*v*br/denom
          
        case(3)        ! e3 < ei < e4
          e43=e(4)-e(3)
          e42=e(4)-e(2)
          e41=e(4)-e(1)
          eij=e(4)-ei
          denom=e41*e42*e43
          dos1t = 3.0d0*v*eij*eij/denom
          
      end select
      
      end function dos1t

 
!EOC
