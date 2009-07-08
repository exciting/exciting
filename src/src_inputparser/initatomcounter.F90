



subroutine initatomcounters

  use modinput
  use mod_atoms
  nspecies=size(input%structure%speciesarray)
  allocate(natoms(nspecies))

  Do is=1,nspecies
    natoms(is)=size(input%structure%speciesarray(is)%species%atomarray)
  end do

end subroutine

