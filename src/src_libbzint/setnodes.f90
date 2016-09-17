!BOP
!
! !ROUTINE: setnodes
!
! !INTERFACE:
      subroutine setnodes
      
! !DESCRIPTION:
!
!This subroutine calculates the coordinates of the nodes of the
!integration polyhedron for convolutions within one tetrahedron.
!
! !USES:
      
      use polyhedron

! !LOCAL VARIABLES:

      implicit none
      integer(4) :: info
 
! !EXTERNAL ROUTINES:
      
      external surfnodes
      external unrepnodes
      external relnodes
      external sortnodes
      
! !REVISION HISTORY:
!      
! Created 1st. April 2004 by RGA
! Last Revision: 22st. April 2008 by RGA

!EOP
!BOC     
!     calculate the intersections between three planes that belongs to 
!     the surface of the tetrahedron
        
      call surfnodes(0)

!     Eliminate repetitions and asign the corresponding value to ntype
!
      call unrepnodes(0)

!     Select those nodes that surround the integration region (e<=ef and
!     f>=ef)
!
      call relnodes(0,info)
      if (info.ne.0)call setnod_debugerr
!      
!     Sort the nodes in the way needed by convw1t.      
!
      call sortnodes(info)
      if (info.ne.0)call setnod_debugerr
      
!EOC      
      contains
      
!BOP
!
! !ROUTINE: setnod_debugerr
!
! !INTERFACE:
      subroutine setnod_debugerr
      
! !DESCRIPTION:
!
! Internal subroutine, repeats the setnodes cycle with write option
! for debugging
!
! !REVISION HISTORY:
!      
! Created: 22st. April 2008 by RGA

!EOP
!BOC 
      
      
!      
! call all the subroutines again with the write option, for debuging
!
        call surfnodes(1)
        call unrepnodes(1)
        call relnodes(1,info)
        stop 'error in setnodes'
        
      end subroutine setnod_debugerr  
!EOC  
      end subroutine setnodes
      
