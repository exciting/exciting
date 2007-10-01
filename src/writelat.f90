
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writelat
use modmain
implicit none
open(50,file='LATTICE'//trim(filext),action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("+----------------------------+")')
write(50,'("| Real-space lattice vectors |")')
write(50,'("+----------------------------+")')
write(50,*)
write(50,'("vector a1 : ",3G18.10)') avec(:,1)
write(50,'("vector a2 : ",3G18.10)') avec(:,2)
write(50,'("vector a3 : ",3G18.10)') avec(:,3)
write(50,*)
write(50,'("Stored column-wise as a matrix :")')
write(50,'(3G18.10)') avec(1,:)
write(50,'(3G18.10)') avec(2,:)
write(50,'(3G18.10)') avec(3,:)
write(50,*)
write(50,'("Inverse of matrix :")')
write(50,'(3G18.10)') ainv(1,:)
write(50,'(3G18.10)') ainv(2,:)
write(50,'(3G18.10)') ainv(3,:)
write(50,*)
write(50,'("Unit cell volume : ",G18.10)') omega
write(50,*)
write(50,*)
write(50,'("+----------------------------------+")')
write(50,'("| Reciprocal-space lattice vectors |")')
write(50,'("+----------------------------------+")')
write(50,*)
write(50,'("vector b1 : ",3G18.10)') bvec(:,1)
write(50,'("vector b2 : ",3G18.10)') bvec(:,2)
write(50,'("vector b3 : ",3G18.10)') bvec(:,3)
write(50,*)
write(50,'("Stored column-wise as a matrix :")')
write(50,'(3G18.10)') bvec(1,:)
write(50,'(3G18.10)') bvec(2,:)
write(50,'(3G18.10)') bvec(3,:)
write(50,*)
write(50,'("Inverse of matrix :")')
write(50,'(3G18.10)') binv(1,:)
write(50,'(3G18.10)') binv(2,:)
write(50,'(3G18.10)') binv(3,:)
write(50,*)
write(50,'("Brillouin zone volume : ",G18.10)') (twopi**3)/omega
close(50)
return
end subroutine

