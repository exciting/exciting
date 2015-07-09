!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: checkmt
! !INTERFACE:
!
!
Subroutine checkmt
    ! !USES:
    Use modinput
    Use modmain
    ! !DESCRIPTION:
    !   Checks for overlapping muffin-tins. If any muffin-tins are found to
    !   intersect the program is terminated with error.
    !
    ! !REVISION HISTORY:
    !   Created May 2003 (JKD)
    !   Revised October 2014 (PAS)
    !EOP
    !BOC
    Implicit None
    ! local variables
    Integer :: ia, is, ja, js
    Integer :: i1, i2, i3
    Real (8) :: v1(3), v2(3), v3(3)
    Real (8) :: t1, t2
    !
    Do i1 = -1, 1
      Do i2 = -1, 1
        Do i3 = -1, 1
          v1(:) = dble(i1) * input % structure % crystal % basevect(:, 1) &
         &      + dble(i2) * input % structure % crystal % basevect(:, 2) &
         &      + dble(i3) * input % structure % crystal % basevect(:, 3)
          Do is = 1, nspecies
            Do ia = 1, natoms(is)
              v2(:) = v1(:) + atposc(:, ia, is)
              Do js = 1, nspecies
                t1 = (rmt(is) + rmt(js)) ** 2
                Do ja = 1, natoms(js)
                  If ( (i1 .ne.  0) .or. (i2 .ne.  0) .or. (i3 .ne. 0) .Or. &
                    &  (is .ne. js) .or. (ia .ne. ja) ) Then
                    v3(:) = v2(:) - atposc(:, ja, js)
                    t2 = v3(1)**2 + v3(2)**2 + v3(3)**2
                    If (t1 .gt. t2) Then
                      Write (*, *)
                      Write (*, '("   Error(checkmt): muffin-tin spheres overlap between")')
                      Write (*, '(T20, "species ", I4, " atom ", I4, " and")') is, ia
                      Write (*, '(T20, "species ", I4, " atom ", I4)') js, ja
                      Write (*, *)
                      Write (*, '("   Sum of muffin-tin radii : ", G13.5)') Sqrt (t1)
                      Write (*, '("   Distance between atoms  : ", G13.5)') Sqrt (t2)
                      Write (*, *)
                      Write (*, '("   Please check, and possibly reduce, the input (or default) values")')
                      Write (*, '("   of the attribute ''''rmt'''' inside the element ''''species'''' ")')
                      Write (*, '("   for species ", I4, " and ", I4)') is, js
                      Write (*, *)
                      Write (*, '("   See link: http://exciting-code.org/ref:structure#species ")')
                      Write (*, *)
                      Stop
                    End If
                  End If
                End Do
              End Do
            End Do
          End Do
        End Do
      End Do
    End Do
    Return
End Subroutine
!EOC
