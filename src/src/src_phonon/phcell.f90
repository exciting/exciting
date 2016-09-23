!
!
!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine phcell (iph, dph, iq, is, ia, ip)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iph
      Real (8), Intent (In) :: dph
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
! local variables
      Integer :: js, ja, na, i, n, iv (3)
      Integer :: i1, i2, i3, m (3, 3)
      Real (8) :: v1 (3), v2 (3), dmin, t1
      If ((iph .Ne. 0) .And. (iph .Ne. 1)) Then
         Write (*,*)
         Write (*, '("Error(phcell): iph should be 0 or 1 : ", I8)') &
        & iph
         Write (*,*)
         Stop
      End If
! check for Gamma-point phonon
      If ((ivq(1, iq) .Eq. 0) .And. (ivq(2, iq) .Eq. 0) .And. (ivq(3, &
     & iq) .Eq. 0)) Then
         m (:, :) = 0
         m (1, 1) = 1
         m (2, 2) = 1
         m (3, 3) = 1
         nphcell = 1
         Go To 10
      End If
! find the first lattice vector
      dmin = 1.d8
      Do i1 = - ngridq (1), ngridq (1)
         Do i2 = - ngridq (2), ngridq (2)
            Do i3 = - ngridq (3), ngridq (3)
               t1 = dble (i1) * vql (1, iq) + dble (i2) * vql (2, iq) + &
              & dble (i3) * vql (3, iq)
               If (Abs(t1-Nint(t1)) .Lt. input%structure%epslat) Then
                  v1 (:) = dble (i1) * avec0 (:, 1) + dble (i2) * avec0 &
                 & (:, 2) + dble (i3) * avec0 (:, 3)
                  t1 = Sqrt (v1(1)**2+v1(2)**2+v1(3)**2)
                  If ((t1 .Lt. dmin) .And. (t1 .Gt. &
                 & input%structure%epslat)) Then
                     m (1, 1) = i1
                     m (2, 1) = i2
                     m (3, 1) = i3
                     dmin = t1
                  End If
               End If
            End Do
         End Do
      End Do
! find the second lattice vector
      dmin = 1.d8
      Do i1 = - ngridq (1), ngridq (1)
         Do i2 = - ngridq (2), ngridq (2)
            Do i3 = - ngridq (3), ngridq (3)
               t1 = dble (i1) * vql (1, iq) + dble (i2) * vql (2, iq) + &
              & dble (i3) * vql (3, iq)
               If (Abs(t1-Nint(t1)) .Lt. input%structure%epslat) Then
! area defined by first two lattice vectors
                  n = (i2*m(3, 1)-i3*m(2, 1)) ** 2 + (i3*m(1, &
                 & 1)-i1*m(3, 1)) ** 2 + (i1*m(2, 1)-i2*m(1, 1)) ** 2
                  If (n .Ne. 0) Then
                     v1 (:) = dble (i1) * avec0 (:, 1) + dble (i2) * &
                    & avec0 (:, 2) + dble (i3) * avec0 (:, 3)
                     t1 = v1 (1) ** 2 + v1 (2) ** 2 + v1 (3) ** 2
                     If (t1 .Lt. dmin) Then
                        m (1, 2) = i1
                        m (2, 2) = i2
                        m (3, 2) = i3
                        dmin = t1
                     End If
                  End If
               End If
            End Do
         End Do
      End Do
! find the third lattice vector
      nphcell = 0
      dmin = 1.d8
      Do i1 = - ngridq (1), ngridq (1)
         Do i2 = - ngridq (2), ngridq (2)
            Do i3 = - ngridq (3), ngridq (3)
               t1 = dble (i1) * vql (1, iq) + dble (i2) * vql (2, iq) + &
              & dble (i3) * vql (3, iq)
               If (Abs(t1-Nint(t1)) .Lt. input%structure%epslat) Then
! number of primitive unit cells in supercell
                  n = m (1, 2) * (i2*m(3, 1)-i3*m(2, 1)) + m (2, 2) * &
                 & (i3*m(1, 1)-i1*m(3, 1)) + m (3, 2) * (i1*m(2, &
                 & 1)-i2*m(1, 1))
                  If (n .Ne. 0) Then
                     v1 (:) = dble (i1) * avec0 (:, 1) + dble (i2) * &
                    & avec0 (:, 2) + dble (i3) * avec0 (:, 3)
                     t1 = v1 (1) ** 2 + v1 (2) ** 2 + v1 (3) ** 2
                     If (t1 .Lt. dmin) Then
                        nphcell = Abs (n)
                        m (1, 3) = i1
                        m (2, 3) = i2
                        m (3, 3) = i3
                        dmin = t1
                     End If
                  End If
               End If
            End Do
         End Do
      End Do
      If (nphcell .Eq. 0) Then
         Write (*,*)
         Write (*, '("Error(phcell): unable to generate supercell")')
         Write (*,*)
         Stop
      End If
10    Continue
! new lattice vectors
      Do i = 1, 3
         input%structure%crystal%basevect(:, i) = dble (m(1, i)) * &
        & avec0 (:, 1) + dble (m(2, i)) * avec0 (:, 2) + dble (m(3, i)) &
        & * avec0 (:, 3)
      End Do
! inverse of lattice vector matrix
      Call r3minv (input%structure%crystal%basevect, ainv)
! generate offset vectors for each primitive cell in the supercell
      n = 1
      vphcell (:, 1) = 0.d0
      Do i1 = - ngridq (1), ngridq (1)
         Do i2 = - ngridq (2), ngridq (2)
            Do i3 = - ngridq (3), ngridq (3)
               If (n .Eq. nphcell) Go To 30
               v1 (:) = dble (i1) * avec0 (:, 1) + dble (i2) * avec0 &
              & (:, 2) + dble (i3) * avec0 (:, 3)
               Call r3mv (ainv, v1, v2)
               Call r3frac (input%structure%epslat, v2, iv)
               Call r3mv (input%structure%crystal%basevect, v2, v1)
               Do i = 1, n
                  t1 = Abs (v1(1)-vphcell(1, i)) + Abs &
                 & (v1(2)-vphcell(2, i)) + Abs (v1(3)-vphcell(3, i))
                  If (t1 .Lt. input%structure%epslat) Go To 20
               End Do
               n = n + 1
               vphcell (:, n) = v1 (:)
20             Continue
            End Do
         End Do
      End Do
      Write (*,*)
      Write (*, '("Error(phcell): unable to generate supercell")')
      Write (*,*)
      Stop
30    Continue
! set up the supercell with a size of nphcell*natoms(js) atoms for each species js
      do js=1, nspecies
        do ja=1,natoms(js)
            deallocate(input%structure%speciesarray(js)%species%atomarray(ja)%atom)
        end do
        deallocate(input%structure%speciesarray(js)%species%atomarray)
        allocate(input%structure%speciesarray(js)%species%atomarray(nphcell*natoms(js)))
        do ja=1,nphcell*natoms(js)
            allocate(input%structure%speciesarray(js)%species%atomarray(ja)%atom)
        end do
      end do
! set up new atomic positions
      Do js = 1, nspecies
         na = 0
         Do ja = 1, natoms (js)
            Do i = 1, nphcell
               na = na + 1
               If (na .Gt. maxatoms) Then
                  Write (*,*)
                  Write (*, '("Error(phcell): too many atoms in superce&
                 &ll : ", I8)') na
                  Write (*, '(" for species ", I4)') js
                  Write (*, '("Adjust maxatoms in modmain and recompile&
                 & code")')
                  Write (*,*)
                  Stop
               End If
               v1 (:) = vphcell (:, i) + atposc0 (:, ja, js)
! add small periodic displacement
               If ((is .Eq. js) .And. (ia .Eq. ja)) Then
                  t1 = dot_product (vqc(:, iq), vphcell(:, i))
                  If (iph .Eq. 0) Then
                     v1 (ip) = v1 (ip) + dph * Cos (t1)
                  Else
                     v1 (ip) = v1 (ip) + dph * Sin (t1)
                  End If
               End If
! convert to new lattice coordinates
               Call r3mv (ainv, v1, input%structure%speciesarray(js)%species%atomarray(na)%atom%coord(:))
               Call r3frac (input%structure%epslat, input%structure%speciesarray(js)%species%atomarray(na)%atom%coord(:), iv)
            End Do
         End Do
         natoms (js) = na
      End Do
! muffin-tin magnetic fields should be zero
      Do js = 1, nspecies
         Do ja = 1, natoms (js)
            input%structure%speciesarray(js)%species%atomarray(ja)%atom%bfcmt(:) = 0.d0
         End Do
      End Do
!
      Return
End Subroutine
