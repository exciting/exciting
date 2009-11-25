!
!
!
! Copyright (C) 2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writelat
      Use modmain
      Use modinput
      Implicit None
      Open (50, File='LATTICE'//trim(filext), Action='WRITE', Form='FOR&
     &MATTED')
      Write (50,*)
      Write (50, '("+----------------------------+")')
      Write (50, '("| Real-space lattice vectors |")')
      Write (50, '("+----------------------------+")')
      Write (50,*)
      Write (50, '("vector a1 : ", 3G18.10)') &
     & input%structure%crystal%basevect(:, 1)
      Write (50, '("vector a2 : ", 3G18.10)') &
     & input%structure%crystal%basevect(:, 2)
      Write (50, '("vector a3 : ", 3G18.10)') &
     & input%structure%crystal%basevect(:, 3)
      Write (50,*)
      Write (50, '("Stored column-wise as a matrix :")')
      Write (50, '(3G18.10)') input%structure%crystal%basevect(1, :)
      Write (50, '(3G18.10)') input%structure%crystal%basevect(2, :)
      Write (50, '(3G18.10)') input%structure%crystal%basevect(3, :)
      Write (50,*)
      Write (50, '("Inverse of matrix :")')
      Write (50, '(3G18.10)') ainv (1, :)
      Write (50, '(3G18.10)') ainv (2, :)
      Write (50, '(3G18.10)') ainv (3, :)
      Write (50,*)
      Write (50, '("Unit cell volume : ", G18.10)') omega
      Write (50,*)
      Write (50,*)
      Write (50, '("+----------------------------------+")')
      Write (50, '("| Reciprocal-space lattice vectors |")')
      Write (50, '("+----------------------------------+")')
      Write (50,*)
      Write (50, '("vector b1 : ", 3G18.10)') bvec (:, 1)
      Write (50, '("vector b2 : ", 3G18.10)') bvec (:, 2)
      Write (50, '("vector b3 : ", 3G18.10)') bvec (:, 3)
      Write (50,*)
      Write (50, '("Stored column-wise as a matrix :")')
      Write (50, '(3G18.10)') bvec (1, :)
      Write (50, '(3G18.10)') bvec (2, :)
      Write (50, '(3G18.10)') bvec (3, :)
      Write (50,*)
      Write (50, '("Inverse of matrix :")')
      Write (50, '(3G18.10)') binv (1, :)
      Write (50, '(3G18.10)') binv (2, :)
      Write (50, '(3G18.10)') binv (3, :)
      Write (50,*)
      Write (50, '("Brillouin zone volume : ", G18.10)') (twopi**3) / &
     & omega
      Close (50)
      Return
End Subroutine
