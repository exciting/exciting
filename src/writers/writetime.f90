!
!
!
! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writetime (t, fname)
      Use modmpi
      Implicit None
  ! arguments
      Real (8), Intent (In) :: t
      Character (*), Intent (In) :: fname
      Open (80, File=trim(fname), Form='formatted', Action='write', &
     & Status='replace')
      Write (80,*)
      Write (80, '("Total time elapsed")')
      Write (80, '(" Wall clock seconds : ", es18.10)') t
      Write (80, '(" Wall clock hours   : ", f14.2)') t / 3600.d0
      Write (80, '(" Wall clock days    : ", f14.2)') t / 3600.d0 / &
     & 24.d0
      Write (80,*)
#ifdef MPI
      Write (80, '("Total time elapsed on all processors")')
      Write (80, '(" Wall clock seconds : ", es18.10)') procs * t
      Write (80, '(" Wall clock hours   : ", f14.2)') procs * t / &
     & 3600.d0
      Write (80, '(" Wall clock days    : ", f14.2)') procs * t / &
     & 3600.d0 / 24.d0
      Write (80,*)
#endif
      Close (80)
End Subroutine writetime
