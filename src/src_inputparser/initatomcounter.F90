!
!
!
!
Subroutine initatomcounters
!
      Use modinput
      Use mod_atoms
      nspecies = size (input%structure%speciesarray)
      Allocate (natoms(nspecies))
!
      Do is = 1, nspecies
         natoms (is) = size &
        & (input%structure%speciesarray(is)%species%atomarray)
      End Do
!
End Subroutine
!
