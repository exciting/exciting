
subroutine initatomcounters
  use modinput
  use modspacegroup
  implicit none
  integer::is
  nspecies=size(input%structure%symmetries%WyckoffPositions%wspeciesarray)
  Do is=1, nspecies
    nwpos(is)=size(input%structure%symmetries%WyckoffPositions%wspeciesarray(is)%wspecies%wposarray)
  end do
end subroutine
