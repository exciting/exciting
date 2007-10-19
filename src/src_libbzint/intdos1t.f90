!BOP
!
! !ROUTINE: intdos1t
!
! !INTERFACE:
       real(8) function intdos1t(e,ei,v)
!
! !DESCRIPTION: 
!
!   This subroutine calculates the contribution to the integrated density of states at
!   an energy $\epsilon = $ \verb"ei" of one tetrahedron. The energies at the 
!   corner of the tetrahedron have to be given in increasing order.   
!
! !USES:
             
      implicit none

! !INPUT PARAMETERS:
 
      real(8), intent(in) :: e(4) ! Eigenenergies at the corners of the
!                                    tetrahedron.
 
      real(8), intent(in) :: ei   ! energy at which the contribution to
!                                     the DOS is calculated
       
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
        case(0)          ! all states are above ei
          intdos1t = 0.0d0
          
        case(1)          ! e1 < ei < e2.
          e21=e(2)-e(1)
          e31=e(3)-e(1)
          e41=e(4)-e(1)
          eij=ei-e(1)
          denom=e21*e31*e41
          intdos1t = v*eij*eij*eij/denom
        
        case(2)         ! e2 < ei < e3
          e21=e(2)-e(1)
          e31=e(3)-e(1)
          e41=e(4)-e(1)
          e32=e(3)-e(2)
          e42=e(4)-e(2)
          eij=ei-e(2)
          denom=e32*e42
          br = e21*e21+3.0d0*eij*e21+3.0d0*eij*eij-(e31+e42)*eij*eij*&
     &         eij/denom
          denom = e31*e41
          intdos1t = v*br/denom
          
        case(3)        ! e3 < ei < e4
          e43=e(4)-e(3)
          e42=e(4)-e(2)
          e41=e(4)-e(1)
          eij=e(4)-ei
          denom=e41*e42*e43
          intdos1t = v*(1.0d0-eij*eij*eij/denom)
          
        case(4)        ! e4 < ei
          intdos1t = v
          
      end select
      
      end function intdos1t
!EOC
