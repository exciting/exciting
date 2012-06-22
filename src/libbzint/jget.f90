
! The code was developed at the Fritz Haber Institute, and
! the intellectual properties and copyright of this file
! are with the Max Planck Society. When you use it, please
! cite R. Gomez-Abal, X. Li, C. Ambrosch-Draxl, M. Scheffler,
! Extended linear tetrahedron method for the calculation of q-dependent
! dynamical response functions, to be published in Comp. Phys. Commun. (2010)

!bop
!
! !FUNCTION: jget
!
! !INTERFACE:
      integer(4) function jget(a,ind)

! !DESCRIPTION:
! 
! Reads the value of the bit ind of vector a
!
! !ARGUMENTS:
      implicit none

      integer(4), intent(in) :: a(*)

      integer(4), intent(in) :: ind
! !RETURN VALUE:
!     
!     integer(4) :: jget            

! !LOCAL VARIABLES:

      integer(4) :: i,off

! !SYSTEM ROUTINES:

      intrinsic ibits
      
!EOP
!BOC
      i=(ind-1)/32+1
      off=mod(ind-1,32)
      jget=ibits(a(i),off,1)
      return
      end function jget
!EOC      

