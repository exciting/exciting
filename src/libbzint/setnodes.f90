
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

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
 
! !EXTERNAL ROUTINES:
      
      external surfnodes
      external unrepnodes
      external relnodes
      external sortnodes
      
! !REVISION HISTORY:
!      
! Created 1st. April 2004 by RGA
! Last Revision: 22st. April 2004 by RGA

!EOP
!BOC     
!     calculate the intersections between three planes that belongs to 
!     the surface of the tetrahedron
        
      call surfnodes

!     Eliminate repetitions and asign the corresponding value to ntype
!
      call unrepnodes

!     Select those nodes that surround the integration region (e<=ef and
!     f>=ef)
!
      call relnodes
!     Sort the nodes in the way needed by convw1t.      
!
      call sortnodes
  
      end subroutine setnodes
      
!EOC
