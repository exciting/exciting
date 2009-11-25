!
!
!
!
Subroutine orthogonalise (n, nstfv, evecfv, ldv, s)
      Implicit None
      Integer, Intent (In) :: n, nstfv, ldv
      Complex (8), Intent (Inout) :: evecfv (ldv, nstfv)
      Complex (8), Intent (In) :: s (n, nstfv)
!
      Integer :: ist, jst, i
      Complex (8) :: zt1
      Complex (8), External :: zdotc
      Real (8) :: t1
! perform Gram-Schmidt orthonormalisation
      Do ist = 1, nstfv
         Do jst = 1, ist - 1
            zt1 = - zdotc (n, evecfv(1, jst), 1, s(1, ist), 1)
            Call zaxpy (n, zt1, evecfv(1, jst), 1, evecfv(1, ist), 1)
         End Do
      End Do
End Subroutine
