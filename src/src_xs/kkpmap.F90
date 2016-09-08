! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: kkpmap
! !INTERFACE:
subroutine kkpmap(ikkp, nkp, ik, ikp)
! !INPUT/OUTPUT PARAMETERS:
! IN
!   ikkp, integer: Index of unique k-point combination 
!   nkp, integer:  Number of k-points
! OUT
!   ik, integer:  Index of first k-point
!   ikp, integer: Index of second k-point
!
! !DESCRIPTION:
!   Calculates individual k-point indices for a combined
!   k-k' index. \\
!   Example: nkp=3 \\ 
!   Number of unique combinations 6\\
!            ikkp   ik    ikp \\
!             1     1     1  \\
!             2     1     2  \\
!             3     1     3  \\
!             4     2     2  \\
!             5     2     3  \\
!             6     3     3 
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC      
  implicit none

  ! Arguments
  integer, intent(in) :: ikkp, nkp
  integer, intent(out) :: ik, ikp

  ik = ceiling(0.5d0 + nkp - sqrt((0.5d0 + nkp)**2 - 2.d0*ikkp))
  ikp = ikkp + (ik * (ik-1) ) / 2 - nkp * (ik-1)
end subroutine kkpmap
!EOC

!BOP
! !ROUTINE: kkpmap_back
! !INTERFACE:
subroutine kkpmap_back(ikkp, nkp, ik, ikp)
! !INPUT/OUTPUT PARAMETERS:
! IN
!   nkp, integer:  Number of k-points
!   ik, integer:  Index of first k-point
!   ikp, integer: Index of second k-point
! OUT
!   ikkp, integer: Index of unique k-point combination 
!
! !DESCRIPTION:
!   Calculates the combined index given to k indices with
!   ikp>=ik.
!   Example: nkp=3 \\ 
!   Number of unique combinations 6\\
!            ikkp   ik    ikp \\
!             1     1     1  \\
!             2     1     2  \\
!             3     1     3  \\
!             4     2     2  \\
!             5     2     3  \\
!             6     3     3 
!
! !REVISION HISTORY:
!   Created as complement to kkpmap. (Aurich)
!EOP
!BOC      
  implicit none

  ! Arguments
  integer, intent(out) :: ikkp
  integer, intent(in) :: ik, ikp, nkp

  ikkp = (ik-1)*nkp - (ik-2)*(ik-1)/2 + (ikp-ik) + 1

end subroutine kkpmap_back
!EOC
