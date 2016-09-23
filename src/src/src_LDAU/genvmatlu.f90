!
!
!
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genvmatlu
      Use modmain
! local variables
      Integer :: ias
      Real (8) :: t1, t2
! allocatable arrays
      Real (8), Allocatable :: enfll (:)
      Complex (8), Allocatable :: vmfll (:, :, :, :, :)
! fully localised limit (FLL) or around mean field (AFM)
      If ((ldapu .Eq. 1) .Or. (ldapu .Eq. 2)) Then
         Call genvmatlu_12
         Return
      End If
! interpolation between the two (PRB 67, 153106 (2003))
      If (ldapu .Eq. 3) Then
         Allocate (enfll(natmtot))
         Allocate (vmfll(lmmaxlu, lmmaxlu, nspinor, nspinor, natmtot))
         ldapu = 1
         Call genvmatlu_12
         vmfll (:, :, :, :, :) = vmatlu (:, :, :, :, :)
         enfll (:) = engyalu (:)
         ldapu = 2
         Call genvmatlu_12
! reset ldapu value
         ldapu = 3
         Do ias = 1, natmtot
            t1 = alphalu (ias)
            t2 = 1.d0 - t1
            vmatlu (:, :, :, :, ias) = t1 * vmfll (:, :, :, :, ias) + &
           & t2 * vmatlu (:, :, :, :, ias)
            engyalu (ias) = t1 * enfll (ias) + t2 * engyalu (ias)
         End Do
         Deallocate (enfll, vmfll)
         Return
      End If
      Write (*,*)
      Write (*, '("Error(genvmatlu): invalid ldapu : ", I8)') ldapu
      Write (*,*)
      Stop
End Subroutine
!
!
Subroutine genvmatlu_12
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ia, ias, ispn, jspn
      Integer :: l, m1, m2, m3, m4, nm
      Integer :: lm1, lm2, lm3, lm4
      Real (8) :: u, j, n, n0
      Real (8) :: mg (3), mg0 (3), mg2
      Real (8) :: edc, sum1, sum2
      Complex (8) zt1, zt2
! automatic arrays
      Real (8) :: vee &
     & (-lmaxlu:lmaxlu,-lmaxlu:lmaxlu,-lmaxlu:lmaxlu,-lmaxlu:lmaxlu)
      Complex (8) dm (lmmaxlu, lmmaxlu, nspinor, nspinor)
      Complex (8) dmt (nspinor, nspinor)
! zero the LDA+U potential for each atom
      vmatlu (:, :, :, :, :) = 0.d0
! zero the LDA+U energy for each atom
      engyalu (:) = 0.d0
! zero the interpolation constants
      alphalu (:) = 0.d0
! begin loop over species
      Do is = 1, nspecies
! define LDA+U parameters
         l = llu (is)
         If (l .Lt. 0) Go To 10
         nm = 2 * l + 1
         u = ujlu (1, is)
         j = ujlu (2, is)
         If ((Abs(u) .Lt. 1.d-10) .And. (Abs(j) .Lt. 1.d-10)) Go To 10
! calculate the Coulomb matrix elements
         Call genveelu (l, u, j, lmaxlu, vee)
! begin loop over atoms
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! copy the density matrix for this atom
            dm (:, :, :, :) = dmatlu (:, :, :, :, ias)
! trace of density matrix for each spin
            dmt (:, :) = 0.d0
            Do ispn = 1, nspinor
               Do jspn = 1, nspinor
                  Do m1 = - l, l
                     lm1 = idxlm (l, m1)
                     dmt (ispn, jspn) = dmt (ispn, jspn) + dm (lm1, &
                    & lm1, ispn, jspn)
                  End Do
               End Do
            End Do
! trace over spin
            n = dble (dmt(1, 1))
            If (associated(input%groundstate%spin)) n = n + dble &
           & (dmt(2, 2))
            n0 = n / dble (nspinor*nm)
! magnetisation
            If (associated(input%groundstate%spin)) Then
               mg (:) = 0.d0
               mg (3) = dble (dmt(1, 1)-dmt(2, 2))
! non-collinear terms
               If (ncmag) Then
                  mg (1) = dble (dmt(1, 2)+dmt(2, 1))
                  mg (2) = dble (zi*(dmt(1, 2)-dmt(2, 1)))
               End If
               mg0 (:) = mg (:) / dble (nspinor*nm)
            End If
! around mean field (AFM) approach
            If (ldapu .Eq. 2) Then
! modify density matrices
               Do m1 = - l, l
                  lm1 = idxlm (l, m1)
                  If (associated(input%groundstate%spin)) Then
                     dm (lm1, lm1, 1, 1) = dm (lm1, lm1, 1, 1) - &
                    & (n0+mg0(3))
                     dm (lm1, lm1, 2, 2) = dm (lm1, lm1, 2, 2) - &
                    & (n0-mg0(3))
! non-collinear terms
                     If (ncmag) Then
                        dm (lm1, lm1, 1, 2) = dm (lm1, lm1, 1, 2) - &
                       & (mg0(1)-zi*mg0(2))
                        dm (lm1, lm1, 2, 1) = dm (lm1, lm1, 2, 1) - &
                       & (mg0(1)+zi*mg0(2))
                     End If
                  Else
! spin-unpolarised case
                     dm (lm1, lm1, 1, 1) = dm (lm1, lm1, 1, 1) - n0
                  End If
               End Do
! determine alpha (PRB 67,153106 (2003))
               sum1 = 0.d0
               Do ispn = 1, nspinor
                  Do m1 = - l, l
                     lm1 = idxlm (l, m1)
                     Do jspn = 1, nspinor
                        Do m2 = - l, l
                           lm2 = idxlm (l, m2)
                           sum1 = sum1 + dble (dm(lm1, lm2, ispn, &
                          & jspn)*dm(lm2, lm1, jspn, ispn))
                        End Do
                     End Do
                  End Do
               End Do
               If (associated(input%groundstate%spin)) Then
                  mg2 = mg (3) ** 2
                  If (ncmag) mg2 = mg2 + mg (1) ** 2 + mg (2) ** 2
               Else
                  mg2 = 0.d0
               End If
               sum2 = n * (1.d0-0.5d0*n/dble(nm)) - 0.5d0 * mg2 / dble &
              & (nm)
               If (Abs(sum2) .Gt. 1.d-14) Then
                  alphalu (ias) = sum1 / sum2
               Else
                  alphalu (ias) = 0.d0
               End If
            End If
! calculation of LDA+U potential and energy
! begin loops over m1 and m2
            Do m1 = - l, l
               lm1 = idxlm (l, m1)
               Do m2 = - l, l
                  lm2 = idxlm (l, m2)
! begin loops over m3 and m4
                  Do m3 = - l, l
                     lm3 = idxlm (l, m3)
                     Do m4 = - l, l
                        lm4 = idxlm (l, m4)
                        Do ispn = 1, nspinor
                           Do jspn = 1, nspinor
                              zt1 = dm (lm2, lm1, ispn, ispn) * dm &
                             & (lm4, lm3, jspn, jspn)
                              zt2 = dm (lm4, lm1, jspn, ispn) * dm &
                             & (lm2, lm3, ispn, jspn)
                              engyalu (ias) = engyalu (ias) + dble &
                             & (zt1-zt2) * vee (m1, m3, m2, m4)
                              vmatlu (lm1, lm2, ispn, ispn, ias) = &
                             & vmatlu (lm1, lm2, ispn, ispn, ias) + dm &
                             & (lm4, lm3, jspn, jspn) * vee (m1, m3, &
                             & m2, m4)
                              vmatlu (lm1, lm4, ispn, jspn, ias) = &
                             & vmatlu (lm1, lm4, ispn, jspn, ias) - dm &
                             & (lm2, lm3, ispn, jspn) * vee (m1, m3, &
                             & m2, m4)
                           End Do
                        End Do
! end loops over m3 and m4
                     End Do
                  End Do
! end loops over m1 and m2
               End Do
            End Do
! multiply energy by factor 1/2
            engyalu (ias) = 0.5d0 * engyalu (ias)
! fully localised limit (FLL) approach: double counting corrections
            If (ldapu .Eq. 1) Then
               If (associated(input%groundstate%spin)) Then
! spin-polarised
                  If (ncmag) Then
! non-collinear case
! correction to the energy
                     edc = 0.5d0 * u * n * (n-1.d0)
                     edc = edc - 0.5d0 * j * dble (dmt(1, 1)*(dmt(1, &
                    & 1)-1.d0))
                     edc = edc - 0.5d0 * j * dble (dmt(2, 2)*(dmt(2, &
                    & 2)-1.d0))
                     edc = edc - 0.5d0 * j * dble (dmt(1, 2)*dmt(2, 1))
                     edc = edc - 0.5d0 * j * dble (dmt(2, 1)*dmt(1, 2))
! correction to the potential
                     Do m1 = - l, l
                        lm1 = idxlm (l, m1)
                        vmatlu (lm1, lm1, 1, 1, ias) = vmatlu (lm1, &
                       & lm1, 1, 1, ias) - u * (n-0.5d0) + j * (dmt(1, &
                       & 1)-0.5d0)
                        vmatlu (lm1, lm1, 2, 2, ias) = vmatlu (lm1, &
                       & lm1, 2, 2, ias) - u * (n-0.5d0) + j * (dmt(2, &
                       & 2)-0.5d0)
                        vmatlu (lm1, lm1, 1, 2, ias) = vmatlu (lm1, &
                       & lm1, 1, 2, ias) + j * dmt (1, 2)
                        vmatlu (lm1, lm1, 2, 1, ias) = vmatlu (lm1, &
                       & lm1, 2, 1, ias) + j * dmt (2, 1)
                     End Do
                  Else
! collinear case
! correction to the energy
                     edc = 0.5d0 * u * n * (n-1.d0)
                     edc = edc - 0.5d0 * j * dble (dmt(1, 1)*(dmt(1, &
                    & 1)-1.d0))
                     edc = edc - 0.5d0 * j * dble (dmt(2, 2)*(dmt(2, &
                    & 2)-1.d0))
! correction to the potential
                     Do m1 = - l, l
                        lm1 = idxlm (l, m1)
                        vmatlu (lm1, lm1, 1, 1, ias) = vmatlu (lm1, &
                       & lm1, 1, 1, ias) - u * (n-0.5d0) + j * (dmt(1, &
                       & 1)-0.5d0)
                        vmatlu (lm1, lm1, 2, 2, ias) = vmatlu (lm1, &
                       & lm1, 2, 2, ias) - u * (n-0.5d0) + j * (dmt(2, &
                       & 2)-0.5d0)
                     End Do
                  End If
               Else
!spin-unpolarised
! correction to the energy
                  edc = 0.5d0 * u * n * (n-1.d0)
                  edc = edc - 0.5d0 * j * n * (n-1.d0)
! correction to the potential
                  Do m1 = - l, l
                     lm1 = idxlm (l, m1)
                     vmatlu (lm1, lm1, 1, 1, ias) = vmatlu (lm1, lm1, &
                    & 1, 1, ias) - u * (n-0.5d0) + j * (n-0.5d0)
                  End Do
               End If
               engyalu (ias) = engyalu (ias) - edc
            End If
! trace of dmatlu times vmatlu
            sum1 = 0.d0
            Do ispn = 1, nspinor
               Do m1 = - l, l
                  lm1 = idxlm (l, m1)
                  Do jspn = 1, nspinor
                     Do m2 = - l, l
                        lm2 = idxlm (l, m2)
                        sum1 = sum1 + dble (dm(lm1, lm2, ispn, &
                       & jspn)*vmatlu(lm2, lm1, jspn, ispn, ias))
                     End Do
                  End Do
               End Do
            End Do
! subtract contribution to the energy of LDA+U potential
            engyalu (ias) = engyalu (ias) - sum1
! end loop over atoms
         End Do
10       Continue
! end loop over species
      End Do
! symmetrise the potential
      Call symdmat (lmaxlu, lmmaxlu, vmatlu)
      Return
End Subroutine
