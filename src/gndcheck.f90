
! Copyright (C) 2014 A. Gulans, C. Draxl
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine gndcheck
      Use modinput
!
! !DESCRIPTION:
! Tests whether selected Hamiltonian,methods and/or algorithms are compatible 
! among themselves. 
!
! !REVISION HISTORY:
!   Created July 2014 (Andris)
!EOP
!BOC

      Implicit None
      if ((input%groundstate%ValenceRelativity.ne.'none').and. &
          (input%groundstate%ValenceRelativity.ne.'zora').and. &
          (input%groundstate%ValenceRelativity.ne.'iora*')) then
        write(*,*) 'ValenceRelativity=', input%groundstate%ValenceRelativity,' is not supported beyond the spherical grid calculations'
        stop
      else
        if (input%groundstate%ValenceRelativity.eq.'iora*') then
          if (.not.(input%groundstate%SymmetricKineticEnergy)) then
            write(*,*) 'ValenceRelativity=', input%groundstate%ValenceRelativity,' is not supported in the non-symmetric kinetic energy mode'
            stop
          endif 
          If (( associated(input%groundstate%spin)) .And. (input%groundstate%ldapu.ne.'none')) Then
            write(*,*) 'ValenceRelativity=',input%groundstate%ValenceRelativity, 'is not supported for spin-resolved and DFT+U calculations'
            stop
          endif

        endif
      endif



      Return
End Subroutine
