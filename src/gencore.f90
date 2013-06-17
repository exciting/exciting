!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: gencore
! !INTERFACE:
!
!
Subroutine gencore
! !USES:
      Use modinput
      Use modmain
! !DESCRIPTION:
!   Computes the core radial wavefunctions, eigenvalues and densities. The
!   radial Dirac equation is solved in the spherical part of the effective
!   potential to which the atomic potential has been appended for
!   $r>R^{\rm MT}$. In the case of spin-polarised calculations, the effective
!   magnetic field is ignored.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! local variables
      Integer :: is, ia, ja, ias, jas, ist, ir
      Real (8) :: t1
      logical dirac_eq
! automatic arrays
      Logical :: done (natmmax)
      Real (8) :: vr (spnrmax)
      

      dirac_eq=(input%groundstate%CoreRelativity.eq."dirac")
     
      Do is = 1, nspecies
         done (:) = .False.
         Do ia = 1, natoms (is)
            If ( .Not. done(ia)) Then
               ias = idxas (ia, is)
               vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
               If (input%groundstate%frozencore) Then
! use atomic potential for the frozen core approximation
                  vr (1:nrmt(is)) = spvr (1:nrmt(is), is)
               Else
! else use the spherical part of the crystal effective potential
                  vr (1:nrmt(is)) = veffmt (1, 1:nrmt(is), ias) * y00
               End If
! append the effective potential from the atomic calculation
               t1 = vr (nrmt(is)) - spvr (nrmt(is), is)
               Do ir = nrmt (is) + 1, spnr (is)
                  vr (ir) = spvr (ir, is) + t1
               End Do
               rhocr (:, ias) = 0.d0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ir,t1)
!$OMP DO
               Do ist = 1, spnst (is)
                  If (spcore(ist, is)) Then
! solve the Dirac equation
                     Call rdirac (0,spn(ist, is), spl(ist, is), spk(ist, &
                    & is), spnr(is), spr(:, is), vr, evalcr(ist, ias), rwfcr(:, 1, ist, ias), &
                    & rwfcr(:, 2, ist, ias),dirac_eq,.false.)
                     t1 = spocc (ist, is)
!$OMP CRITICAL
                     if (dirac_eq) then
                       Do ir = 1, spnr (is)
! add to the core density
                           rhocr (ir, ias) = rhocr (ir, ias) + t1 * &
                         & (rwfcr(ir, 1, ist, ias)**2+rwfcr(ir, 2, ist, &
                         & ias)**2)
                       End Do
                     else
                       Do ir = 1, spnr (is)
                           rhocr (ir, ias) = rhocr (ir, ias) + t1 * &
                         & rwfcr(ir, 1, ist, ias)**2
                       Enddo
                     endif
!$OMP END CRITICAL
                  End If
               End Do
!$OMP END DO
!$OMP END PARALLEL
               Do ir = 1, spnr (is)
                  rhocr (ir, ias) = rhocr (ir, ias) / (fourpi*spr(ir, &
                 & is)**2)
               End Do
               done (ia) = .True.
! copy to equivalent atoms
               Do ja = 1, natoms (is)
                  If (( .Not. done(ja)) .And. (eqatoms(ia, ja, is))) &
                 & Then
                     jas = idxas (ja, is)
                     Do ist = 1, spnst (is)
                        If (spcore(ist, is)) Then
                           evalcr (ist, jas) = evalcr (ist, ias)
                           rwfcr (1:spnr(is), :, ist, jas) = rwfcr &
                          & (1:spnr(is), :, ist, ias)
                        End If
                     End Do
                     rhocr (1:spnr(is), jas) = rhocr (1:spnr(is), ias)
                     done (ja) = .True.
                  End If
               End Do
            End If
         End Do
      End Do
      Return
End Subroutine
!EOC
