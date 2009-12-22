!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: autoradmt
! !INTERFACE:
!
!
Subroutine autoradmt
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Automatically determines the muffin-tin radii from the formula
!   $$ R_i\propto 1+\zeta|Z_i|^{1/3}, $$
!   where $Z_i$ is the atomic number of the $i$th species, $\zeta$ is a
!   user-supplied constant ($\sim 0.625$). The parameter $\zeta$ is stored in
!   {\tt rmtapm(1)} and the value which governs the distance between the
!   muffin-tins is stored in {\tt rmtapm(2)}. When {\tt rmtapm(2)} $=1$, the
!   closest muffin-tins will touch.
!
! !REVISION HISTORY:
!   Created March 2005 (JKD)
!   Changed the formula, September 2006 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, js, ia, ja, i1, i2, i3
      Real (8) :: s, v1 (3), v2 (3), t1, t2, t3
! external functions
      Real (8) :: r3dist
      External r3dist
! initial muffin-tin radii
      Do is = 1, nspecies
         rmt (is) = 1.d0 + input%groundstate%rmtapm(1) * Abs (spzn(is)) &
        & ** (1.d0/3.d0)
      End Do
! determine scaling factor
      s = 1.d8
      Do i1 = - 1, 1
         Do i2 = - 1, 1
            Do i3 = - 1, 1
               v1 (:) = dble (i1) * input%structure%crystal%basevect(:, &
              & 1) + dble (i2) * input%structure%crystal%basevect(:, 2) &
              & + dble (i3) * input%structure%crystal%basevect(:, 3)
               Do is = 1, nspecies
                  Do ia = 1, natoms (is)
                     v2 (:) = v1 (:) + atposc (:, ia, is)
                     Do js = 1, nspecies
                        t1 = 1.d0 / (rmt(is)+rmt(js))
                        Do ja = 1, natoms (js)
                           If ((i1 .Ne. 0) .Or. (i2 .Ne. 0) .Or. (i3 &
                          & .Ne. 0) .Or. (is .Ne. js) .Or. (ia .Ne. &
                          & ja)) Then
                              t2 = r3dist (v2, atposc(:, ja, js))
                              t3 = t1 * t2
                              If (t3 .Lt. s) s = t3
                           End If
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do
      s = s * input%groundstate%rmtapm(2)
! scale all radii
      Do is = 1, nspecies
! limit number of decimal digits
         t1 = s * rmt (is) * 10000.d0
         t1 = dble (Int(t1)) / 10000.d0
         rmt (is) = t1
      End Do
      Return
End Subroutine
!EOC
