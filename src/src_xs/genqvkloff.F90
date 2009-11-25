!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genqvkloff (vq, voff)
      Use modmain
      Use modinput
      Use modxs
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vq (3)
      Real (8), Intent (Out) :: voff (3)
  ! local variables
      Real (8) :: v1 (3)
      Integer :: iv (3)
      If &
     & (any(input%groundstate%vkloff/dble(input%groundstate%ngridk)+vq &
     & .Ge. 1.d0)) Then
     ! vector is outside Brillouine zone
         v1 = input%groundstate%vkloff / dble &
        & (input%groundstate%ngridk) + vq
         Call mapkto01 (v1)
         voff = v1 * dble (input%groundstate%ngridk)
         If (any(v1*dble(input%groundstate%ngridk) .Ge. 1.d0)) Then
            v1 = v1 * dble (input%groundstate%ngridk)
            Call mapkto01 (v1)
            voff = v1
         End If
      Else If &
     & (any(input%groundstate%vkloff+vq*dble(input%groundstate%ngridk) &
     & .Ge. 1.d0)) Then
     ! vector is inside Brillouine zone but outside k-point spacing
         v1 = input%groundstate%vkloff + vq * dble &
        & (input%groundstate%ngridk)
         Call mapkto01 (v1)
         voff = v1
      Else
     ! vector is inside k-point spacing
         voff = input%groundstate%vkloff + vq * &
        & input%groundstate%ngridk
      End If
  ! treatment of values close to zero or one
      Call r3frac (input%structure%epslat, voff, iv)
End Subroutine genqvkloff
!
!
Subroutine mapkto01 (v)
      Implicit None
  ! arguments
      Real (8), Intent (Inout) :: v (3)
  ! local variables
  !integer :: id(3)
      Integer (8) :: v2 (3), v3 (3)
      Real (8), Parameter :: fac = 1.d15
  !call r3frac(epslat,v,id)
      v2 = dint (v)
      v3 = dint (fac*v)
      v3 = v3 - v2 * dint (fac)
      v = dble (v3/dint(fac))
      Where (v .Lt. 0.d0) v = v + 1.d0
End Subroutine mapkto01
