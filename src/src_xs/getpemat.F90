!
!
!
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_getpemat
      Implicit None
Contains
!
!
      Subroutine getpemat (iq, ik, pfilnam, efilnam, m12, m34, p12, &
     & p34)
         Use modmain
         Use modinput
         Use modxs
         Use modtetra
         Use m_getpmat
         Use m_getemat
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, ik
         Character (*), Intent (In) :: pfilnam, efilnam
         Complex (8), Optional, Intent (Out) :: m12 (:, :, :), p34 (:, &
        & :, :)
         Complex (8), Optional, Intent (Out) :: p12 (:, :, :), m34 (:, &
        & :, :)
    ! local variables
         Character (*), Parameter :: thisnam = 'getpemat'
         Real (8) :: fourpisqt
         Integer :: n, igq, j, i1, i2
         Logical :: tq0
         Logical, External :: tqgamma
         tq0 = tqgamma (iq)
         n = ngq (iq)
         fourpisqt = Sqrt (fourpi)
         If (tq0 .And. ( .Not. present(p12))) Then
            Write (*,*)
            Write (*, '("Error(", a, "): Gamma q-point but momentum mat&
           &rix elements not requested.")') thisnam
            Write (*,*)
            Call terminate
         End If
    ! Gamma q-point
         If (tq0) Then
       ! read momentum matrix elements
            Call getpmat (ik, vkl0, istl1, istu1, istl2, istu2, .True., &
           & trim(pfilnam), p12)
            If (present(p34)) Call getpmat (ik, vkl0, istl3, istu3, &
           & istl4, istu4, .True., trim(pfilnam), p34)
       ! consider symmetric gauge wrt. Coulomb potential
       ! (multiply with v^(1/2))
       ! and normalize wrt. KS eigenvalues (no scissors correction!)
            Do j = 1, 3
               Do i1 = 1, nst1
                  Do i2 = 1, nst2
                     If (Abs(deou(i1, i2)) .Ge. input%xs%epsdfde) &
                    & Then
                        p12 (j, i1, i2) = - p12 (j, i1, i2) / deou (i1, &
                       & i2) * fourpisqt
                     Else
                        p12 (j, i1, i2) = zzero
                        If ((Abs(docc12(i1, i2)) .Gt. &
                       & input%groundstate%epsocc).And.(input%xs%dbglev .Gt. 0)) Then
                           Write (*, '("Warning(", a, "): divergent ene&
                          &rgy denominator: q-point, k-point, band indi&
                          &ces 1-2, delta E12, delta Occ:", 4i6, 2g18.10)') &
                          & thisnam, iq, ik, &
                          & i1 + istl1 - 1, i2 + istl2 - 1, deou (i1, &
                          & i2), docc12(i1,i2)
                        End If
                     End If
                     If (present(p34)) Then
                        If (Abs(deuo(i2, i1)) .Ge. &
                       & input%xs%epsdfde) Then
                           p34 (j, i2, i1) = - p34 (j, i2, i1) / deuo &
                          & (i2, i1) * fourpisqt
                        Else
                           p34 (j, i2, i1) = zzero
                           If ((Abs(docc21(i2, i1)) .Gt. &
                          & input%groundstate%epsocc).And.(input%xs%dbglev .Gt. 0)) Then
                              Write (*, '("Warning(", a, "): divergent &
                             &energy denominator: q-point, k-point, ban&
                             &d indices 3-4:", 4i6, g18.10)') thisnam, &
                             & iq, ik, i1 + istl1 - 1, i2 + istl2 - 1, &
                             & deuo (i2, i1)
                           End If
                        End If
                     End If
                  End Do
               End Do
            End Do
         End If
         If (( .Not. tq0) .Or. (n .Gt. 1)) Then
       ! for BSE(-kernel) matrix elements are calculated on the fly
            If (tscreen) Then
               m12 (:, :, :) = xiou (:, :, :)
               If (present(m34)) m34 (:, :, :) = xiuo (:, :, :)
            Else
          ! read matrix elemets of plane wave
               If (present(m34)) Then
                  Call getemat (iq, ik, .True., trim(efilnam), ngq(iq), &
                 & istl1, istu1, istl2, istu2, m12, istl3, istu3, &
                 & istl4, istu4, m34)
               Else
                  Call getemat (iq, ik, .True., trim(efilnam), ngq(iq), &
                 & istl1, istu1, istl2, istu2, m12)
               End If
            End If
       ! consider symmetric gauge wrt. Coulomb potential (multiply with v^(1/2))
            If ( .Not. tq0) Then
               m12 (:, :, 1) = m12 (:, :, 1) *sptclg(1,iq)
               If (present(m34)) m34 (:, :, 1) = m34 (:, :, 1) *sptclg(1,iq)
            End If
            If (n .Gt. 1) Then
               Forall (igq=2:n)
                  m12 (:, :, igq) = m12 (:, :, igq) *sptclg(igq,iq)
               End Forall
               If (present(m34)) Then
                  Forall (igq=2:n)
                     m34 (:, :, igq) = m34 (:, :, igq) *sptclg(igq,iq)
                  End Forall
               End If
            End If
         End If
      End Subroutine getpemat
!
End Module m_getpemat
