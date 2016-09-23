
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine setuphsvect (n, m, system, evecfv, ldv, h, s)
      Use diisinterfaces
      Use modfvsystem
      Use modmain, Only: nstfv, zone, zzero
      Implicit None
      Integer, Intent (In) :: n, m, ldv
      Complex (8), Intent (Inout) :: evecfv (ldv, nstfv)
      Type (evsystem) :: system
      Complex (8), Intent (Out) :: h (n, nstfv), s (n, nstfv)
      Integer :: i
      Complex (8) :: z
      Real (8) :: t
!
      Call zhemm ('L', 'U', n, m, zone, system%overlap%za(1, 1), n, &
     & evecfv(1, 1), ldv, zzero, s(1, 1), n)
      Do i = 1, m
         t = (dble(zdotc(n, evecfv(1, i), 1, s(1, i), 1)))
         z = 1.d0 / Sqrt (t)
   ! write(*,*)"z",z
         Call zscal (n, z, evecfv(1, i), 1)
         Call zscal (n, z, s(1, i), 1)
      End Do
      h = h
   !call zhemm('L','U',n,m,zone,system%hamilton%za(1,1),n,evecfv(1,1),ldv,&
   !  zzero,h(1,1),n)
!
End Subroutine setuphsvect
