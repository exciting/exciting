!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine findgqmap (iq, iqr, nsc, sc, ivgsc, n, isc, isci, ivgu, &
& igqmap)
      Use modmain
      Use modinput
      Use modxs
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, iqr, nsc, sc (maxsymcrys), ivgsc (3, &
     & maxsymcrys), n
      Integer, Intent (Out) :: isc, isci, ivgu (3), igqmap (n)
  ! local variables
      Real (8) :: vqr (3), v2 (3), t1
      Integer :: iqrnr, j, isym, isymi, lspl, lspli, iv (3), ivg1 (3), &
     & igq1
  ! find map from G-vectors to rotated G-vectors
      iqrnr = iqmap (ivqr(1, iqr), ivqr(2, iqr), ivqr(3, iqr))
      vqr (:) = vqlr (:, iqr)
      Do j = 1, nsc
         isym = sc (j)
         lspl = lsplsymc (isym)
         isymi = scimap (isym)
         lspli = lsplsymc (isymi)
         Do igq1 = 1, n
            ivg1 (:) = ivg (:, igqig(igq1, iq))
        ! G1 = si^-1 * ( G + G_s ) , where si is the inverse of s
            iv = matmul (transpose(symlat(:, :, lspli)), ivg1+ivgsc(:, &
           & j))
        ! |G1 + q|
            v2 = matmul (bvec, iv+vqr)
            t1 = Sqrt (sum(v2**2))
            If ((n .Gt. 1) .And. (t1 .Gt. input%xs%gqmax)) Then
               Write (*,*)
               Write (*, '("Info(findgqmap): need one more symmetry operation")')
               Write (*,*)
               Go To 10
            End If
        ! locate G1 + q in G+q-vector set
            igqmap (igq1) = ivgigq (iv(1), iv(2), iv(3), iqrnr)
            If (igqmap(igq1) .Le. 0) Then
               Write (*,*)
               Write (*, '("Error(findgqmap): failed to map rotated G-v&
              &ector")')
               Write (*, '(" non-reduced q-point		       :", i8)') iq
               Write (*, '(" reduced q-point 		       :", i8)') iqr
               Write (*, '(" reduced q-point in non-reduced set     :",&
              & i8)') iqrnr
               Write (*, '(" G+q-vector index (non-reduced q-point) :",&
              & i8)') igq1
               Write (*, '(" rotated G-vector		       :", 3i8)') iv
               Write (*,*)
               Call terminate
            End If
        ! end loop over G
         End Do
     ! store G1 vector
         ivgu (:) = ivgsc (:, j)
         isc = isym
         isci = isymi
         Go To 20
10       Continue
     ! end loop over symmetry operations
      End Do
      Write (*,*)
      Write (*, '("Error(findgqmap): failed to reduce q-point: ", i8)') &
     & iq
      Write (*,*)
      Call terminate
20    Continue
End Subroutine findgqmap
