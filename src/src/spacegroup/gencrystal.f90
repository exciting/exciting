
subroutine gencrystal
use modspacegroup
use modsymmetries
implicit none
! local variables
integer::is, ia, ip, i, j
integer::i1, i2, i3
integer::id(3), ngen, ngrp
real(8)::abr, acr, bcr
real(8)::sab, cab, cac, cbc
real(8)::v1(3), v2(3)
! space group generator Seitz matrices
real(8)::srgen(3, 3, 12), stgen(3, 12)
! space group Seitz matrices
real(8)::srgrp(3, 3, 192), stgrp(3, 192)
! external functions
real(8)::r3taxi
external r3taxi
! convert angles from degrees to radians
abr= symmetries%lattice%ab*(pi/180.d0)
acr= symmetries%lattice%ac*(pi/180.d0)
bcr= symmetries%lattice%bc*(pi/180.d0)
! setup lattice vectors
sab=sin(abr)
if (abs(sab).lt.symmetries%lattice%epslat) then
  write(*, *)
  write(*, '("Error(gencrystal): degenerate lattice vectors")')
  write(*, *)
  stop
end if
cab=cos(abr)
cac=cos(acr)
cbc=cos(bcr)
avecnew(1, 1)= symmetries%lattice%a
avecnew(2, 1)=0.d0
avecnew(3, 1)=0.d0
avecnew(1, 2)= symmetries%lattice%b*cab
avecnew(2, 2)= symmetries%lattice%b*sab
avecnew(3, 2)=0.d0
avecnew(1, 3)= symmetries%lattice%c*cac
avecnew(2, 3)= symmetries%lattice%c*(cbc-cab*cac)/sab
avecnew(3, 3)= symmetries%lattice%c*sqrt(sab**2-cac**2+2.d0*cab*cac*cbc-cbc**2)/sab
do i=1, 3
  do j=1, 3
    if (abs(avecnew(i, j)).lt.symmetries%lattice%epslat) avecnew(i, j)=0.d0
  end do
end do
! scale lattice vectors by the number of unit cells
do i=1, 3
  avecnew(:, i)=avecnew(:, i)*dble( symmetries%lattice%ncell(i))
end do
! determine the Hall symbol from the Hermann-Mauguin symbol
call sgsymb( symmetries%HermannMauguinSymbol, num, schn, hall)
! determine the space group generators
call seitzgen(hall, ngen, srgen, stgen)
! compute the space group operations
call gengroup(ngen, srgen, stgen, ngrp, srgrp, stgrp)
! compute the equivalent atomic positions
do is=1, nspecies
  natoms(is)=0
  do ip=1, nwpos(is)
    do j=1, ngrp
! apply the space group operation
      call r3mv(srgrp(:, 1, j),&
       symmetries%WyckoffPositions%wspeciesarray(is)%wspecies%wposarray(ip)%wpos%coord(:), v1)
      v1(:)=v1(:)+stgrp(:, j)
      do i1=0,  symmetries%lattice%ncell(1)-1
	do i2=0,  symmetries%lattice%ncell(2)-1
	  do i3=0,  symmetries%lattice%ncell(3)-1
	    v2(1)=(v1(1)+dble(i1))/dble( symmetries%lattice%ncell(1))
	    v2(2)=(v1(2)+dble(i2))/dble( symmetries%lattice%ncell(2))
	    v2(3)=(v1(3)+dble(i3))/dble( symmetries%lattice%ncell(3))
	    call r3frac(symmetries%lattice%epslat, v2, id)
! check if new position already exists
	    do ia=1, natoms(is)
	      if (r3taxi(v2, atposlnew(:, ia, is)).lt.symmetries%lattice%epslat) goto 30
	    end do
! add new position to list
	    natoms(is)=natoms(is)+1
	    if (natoms(is).gt.maxatoms) then
	      write(*, *)
	      write(*, '("Error(gencrystal): natoms too large")')
	      write(*, '(" for species ", I4)') is
	      write(*, '("Adjust maxatoms and recompile code")')
	      write(*, *)
	      stop
	    end if
	    atposlnew(:, natoms(is), is)=v2(:)
	  end do
	end do
      end do
30 continue
    end do
  end do
  natmtot=natmtot+natoms(is)
end do
! set magnetic fields to zero
! reduce conventional cell to primitive cell if required
if (symmetries%lattice%primcell) call findprim
! find the total number of atoms
natmtot=0
do is=1, nspecies
  natmtot=natmtot+natoms(is)
end do
! determine the Cartesian atomic coordinates
do is=1, nspecies
  do ia=1, natoms(is)
    call r3mv(avecnew, atposlnew(:, ia, is), atposcnew(:, ia, is))
  end do
end do
return
end subroutine
