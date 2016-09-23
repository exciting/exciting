
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!
!
Subroutine calcupdatevectors (n, iunconverged, P, w, r, evalfv, evecfv, &
& phi)
      Use modmain, Only: nstfv, nmatmax, zzero, zone
      Use diisinterfaces
      Implicit None
      Integer, Intent (In) :: n, iunconverged
      Complex (8), Intent (In) :: P (nmatmax, nmatmax)
      Complex (8), Intent (In) :: r (n, nstfv), evecfv (n, nstfv)
      Complex (8), Intent (Out) :: phi (n, nstfv)
      Real (8), Intent (In) :: w (nmatmax), evalfv (nstfv)
      Complex (8) :: z, v (n, nstfv), alpha
      Real (8) :: nrm
      Integer :: m, i, j
      m = iunconverged
      Call zgemm ('C', 'N', n, m, n, zone, P, nmatmax, r, n, zzero, v, &
     & n)
!
      Do j = 1, m
         Do i = 1, n
            z = dcmplx (w(i)-evalfv(j), 0.0)
     !   write(*,*)z,w(i),evalfv(j),i,j,iunconverged
            If (Abs(z) .Gt. 0.1e-6) Then
               v (i, j) = - v (i, j) / z
            Else
               v (i, j) = zzero
            End If
         End Do
!
      End Do
!
!
      Do i = 1, m
         Call zcopy (n, evecfv(1, i), 1, phi(1, i), 1)
      End Do
!
      Call zgemm ('N', 'N', n, m, n, zone, P, nmatmax, v, n, zone, phi, &
     & n)
  ! alpha=dcmplx(.6,0)
  ! call zscal(n,alpha,phi,1)
  ! call zaxpy(n,1-alpha,evecfv,1,phi,1)
!
      Do i = 1, m
         Call zcopy (n, phi(1, i), 1, evecfv(1, i), 1)
      End Do
!
End Subroutine calcupdatevectors
