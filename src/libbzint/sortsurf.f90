
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!BOP
!
! !ROUTINE: sortsurf
!
! !INTERFACE:
      subroutine sortsurf
      
! !DESCRIPTION:
! This subroutine sort the nodes of Fermi surface in the k-mesh, so that
! if it is a square region, it can be divided into two triangles that the
! nodes 1,2,3 form one triangle and nodes 2,3,4 form the other. It is used 
! for the surface integration.
!
! !USES:
      
      use polyhedron
      
! !LOCAL VARIABLES:

      implicit none

      integer(1) :: ibit, inod, nt
      integer(4) :: iind
      integer(1), dimension(4) :: sp
      
! !INTRINSIC ROUTINES:
      
      intrinsic ibits
      intrinsic btest  

! !REVISION HISTORY:
!
! Created 17th, Jan 2005 by XZL
!
!EOP
!BOC
      sp(1:4)=0
      allocate(index(4,1))
      nt=ntype(1)
      index(1:4,1)=0
      index(1,1)=1
      do inod=2,4 
       do ibit=0,3
         if((btest(nt,ibit)).and.(btest(ntype(inod),ibit)))sp(inod)=sp(inod)+1
       enddo
      enddo
      do inod=2,4
       if(sp(inod).eq.0)index(inod,1)=4
      enddo
      iind=2
      do inod=2,4
       if(index(inod,1).eq.0) then
         index(inod,1)=iind
         iind=iind+1
       else
         continue
       endif
      enddo
        
      end subroutine sortsurf
!EOC            
              
            
      
