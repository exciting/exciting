subroutine initlattice
 use modinput
 use mod_atoms
 use mod_lattice
implicit none
integer ::is, ia, js, ja, i
real(8)::t1,v(3)
! scale the lattice vectors (scaling not referenced again in code)
input%structure%crystal%basevect(:, 1) = input%structure%crystal%stretch(1) * &
input%structure%crystal%basevect(:, 1)
input%structure%crystal%basevect(:, 2) = input%structure%crystal%stretch(2) * &
input%structure%crystal%basevect(:, 2)
input%structure%crystal%basevect(:, 3) = input%structure%crystal%stretch(3) * &
input%structure%crystal%basevect(:, 3)
input%structure%crystal%basevect(:, :) = input%structure%crystal%scale * &
input%structure%crystal%basevect(:, :)
! check if system is an isolated molecule
if (input%structure%molecule) then
! set up cubic unit cell with vacuum region around molecule
  input%structure%crystal%basevect(:, :)=0.d0
  do is=1, nspecies
    do ia=1, natoms(is)
      do js=1, nspecies
	do ja=1, natoms(is)
	  do i=1, 3
	    t1 = abs(input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(i) -&
    &input%structure%speciesarray(js)%species%atomarray(ja)%atom%coord(i))
	    if (t1.gt.input%structure%crystal%basevect(i, i)) input%structure%crystal%basevect(i, i)=t1
	  end do
	end do
      end do
    end do
  end do
  do i=1, 3
    input%structure%crystal%basevect(i, i) = input%structure%crystal%basevect(i, i) + input%structure%vacuum
  end do
! convert atomic positions from Cartesian to lattice coordinates
  call r3minv(input%structure%crystal%basevect, ainv)
  do is=1, nspecies
    do ia=1, natoms(is)
      call r3mv(ainv, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), v)
      input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)=v(:)
    end do
  end do
end if
end subroutine
