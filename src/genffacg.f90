!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genffacg (is, ffacg)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: is
      Real (8), Intent (Out) :: ffacg (ngvec)
! local variables
      Integer :: ig
      Real (8) :: t1, t2, t3, t4
      t1 = fourpi / omega
      t2 = input%groundstate%cfdamp / input%groundstate%gmaxvr
      Do ig = 1, ngvec
         If (gc(ig) .Gt. input%structure%epslat) Then
            If (input%groundstate%cfdamp .Ne. 0.d0) Then
! use damping if required
               t3 = Exp (-(t2*gc(ig))**2)
            Else
               t3 = 1.d0
            End If
            t4 = gc (ig) * rmt (is)
            ffacg (ig) = t1 * t3 * (Sin(t4)-t4*Cos(t4)) / (gc(ig)**3)
         Else
            ffacg (ig) = (t1/3.d0) * rmt (is) ** 3
         End If
      End Do
      Return
End Subroutine
