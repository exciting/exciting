! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl
! adapted for exciting 2014,Stefan Kontur
!
!
subroutine SORT(iidim) 
!
   use raman_ew
   use modinput
   implicit none
   integer :: i,j,ihelp,iidim
   double precision :: ehelp
   double precision, allocatable :: esort(:)
!
   allocate( esort(input%properties%raman%nstate) )
!
!  get nnumber lowest eigenvalues and remember their index
!
      do i = 1, input%properties%raman%nstate
         ehelp = 1.d9
         do j = 1, iidim
            if (i .gt. 1) then
               if (eigen(j) .le. esort(i-1)) cycle
            endif
            if (eigen(j) .lt. ehelp) then
               ehelp = eigen(j)
               ihelp = j
            endif
         enddo
         indexi(i) = ihelp
         esort(i) = ehelp
      enddo
      do i = 1, input%properties%raman%nstate
         eigen(i) = esort(i)
      enddo
   deallocate( esort )
   return
end subroutine sort
!
