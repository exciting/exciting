


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writegeom
! !INTERFACE:


subroutine writegeom(topt)
! !USES:
use modinput
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   topt : if .true. then the filename will be {\tt GEOMETRY_OPT.OUT}, otherwise
!          {\tt GEOMETRY.OUT} (in,logical)
! !DESCRIPTION:
!   Outputs the lattice vectors and atomic positions to file, in a format
!   which may be then used directly in {\tt exciting.in}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: topt
! local variables
integer::is, ia
real(8)::v(3)
if (topt) then
  open(50, file='GEOMETRY_OPT'//trim(filext), action='WRITE', form='FORMATTED')
else
  open(50, file='GEOMETRY'//trim(filext), action='WRITE', form='FORMATTED')
end if
write(50, *)
write(50, '("scale")')
write(50, '(" 1.0")')
write(50, *)
write(50, '("scale1")')
write(50, '(" 1.0")')
write(50, *)
write(50, '("scale2")')
write(50, '(" 1.0")')
write(50, *)
write(50, '("scale3")')
write(50, '(" 1.0")')
write(50, *)
write(50, '("avec")')
write(50, '(3G18.10)') input%structure%crystal%basevect(:, 1)
write(50, '(3G18.10)') input%structure%crystal%basevect(:, 2)
write(50, '(3G18.10)') input%structure%crystal%basevect(:, 3)
if (input%structure%molecule) then
  write(50, *)
  write(50, '("molecule")')
  write(50, '(" ", L1)') input%structure%molecule
end if
write(50, *)
write(50, '("atoms")')
write(50, '(I4, T40, " : nspecies")') nspecies
do is=1, nspecies
  write(50, '(" ''", A, "''", T40, " : spfname")') trim(input%structure%speciesarray(is)%species%speciesfile)
  write(50, '(I4, T40, " : natoms; atpos, bfcmt below")') natoms(is)
  do ia=1, natoms(is)
    if (input%structure%molecule) then
! write Cartesian coordinates for the molecular case
      call r3mv(input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:),&
   &v)
    else
! otherwise write lattice coordinates
      v(:)=input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
    end if
    write(50, '(3F14.8, "  ", 3F12.8)') v(:), input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
  end do
end do
close(50)
return
end subroutine
!EOC
