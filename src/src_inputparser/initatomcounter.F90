
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine initatomcounters
      Use modinput
      Use mod_atoms
      nspecies=0
      if (allocated(natoms)) deallocate(natoms)
      allocate(natoms(nspecies))
      natoms(:)=0
      if (associated(input%structure%speciesarray)) then
        nspecies = size(input%structure%speciesarray)
        if (allocated(natoms)) deallocate(natoms)
        Allocate(natoms(nspecies))
        Do is = 1, nspecies
           natoms(is) = size(input%structure%speciesarray(is)%species%atomarray)
        End Do
      end if
End Subroutine
