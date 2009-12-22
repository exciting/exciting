!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: hmlaa
! !INTERFACE:
!
!
Subroutine hmlaa (tapp, is, ia, ngp, apwalm, v, h)
! !USES:
      Use modinput
      Use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tapp   : .true. if the Hamiltonian is to be applied to the input vector,
!            .false. if the full matrix is to be calculated (in,logical)
!   is     : species number (in,integer)
!   ia     : atom number (in,integer)
!   ngp    : number of G+p-vectors (in,integer)
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   v      : input vector to which H is applied if tapp is .true., otherwise
!            not referenced (in,complex(nmatmax))
!   h      : H applied to v if tapp is .true., otherwise it is the Hamiltonian
!            matrix in packed form (inout,complex(*))
! !DESCRIPTION:
!   Calculates the APW-APW contribution to the Hamiltonian matrix.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Logical, Intent (In) :: tapp
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: v (nmatmax)
      Complex (8), Intent (Inout) :: h (*)
! local variables
      Integer :: ias, io1, io2
      Integer :: l1, l2, l3, m1, m2, m3, lm1, lm2, lm3
      Real (8) :: t1
      Complex (8) zt1, zsum
! automatic arrays
      Complex (8) zv (ngp)
      ias = idxas (ia, is)
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
               zv (:) = 0.d0
               Do l3 = 0, input%groundstate%lmaxapw
                  Do m3 = - l3, l3
                     lm3 = idxlm (l3, m3)
                     If (lm1 .Ge. lm3) Then
                        Do io2 = 1, apword (l3, is)
                           zsum = 0.d0
                           Do l2 = 0, input%groundstate%lmaxvr
                              If (Mod(l1+l2+l3, 2) .Eq. 0) Then
                                 Do m2 = - l2, l2
                                    lm2 = idxlm (l2, m2)
                                    If ((l2 .Eq. 0) .Or. (l1 .Ge. l3)) &
                                   & Then
                                       zt1 = gntyry (lm1, lm2, lm3) * &
                                      & haa (io1, l1, io2, l3, lm2, &
                                      & ias)
                                    Else
                                       zt1 = gntyry (lm1, lm2, lm3) * &
                                      & haa (io1, l3, io2, l1, lm2, &
                                      & ias)
                                    End If
                                    zsum = zsum + zt1
                                 End Do
                              End If
                           End Do
                           If (lm1 .Eq. lm3) zsum = zsum * 0.5d0
                           If (Abs(dble(zsum))+Abs(aimag(zsum)) .Gt. &
                          & 1.d-14) Then
                              Call zaxpy (ngp, zsum, apwalm(:, io2, &
                             & lm3, ias), 1, zv, 1)
                           End If
                        End Do
                     End If
                  End Do
               End Do
               Call zmatinp (tapp, ngp, zone, apwalm(:, io1, lm1, ias), &
              & zv, v, h)
            End Do
         End Do
      End Do
! kinetic surface contribution
      t1 = 0.25d0 * rmt (is) ** 2
      Do l1 = 0, input%groundstate%lmaxmat
         Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword (l1, is)
               Do io2 = 1, apword (l1, is)
                  zt1 = t1 * apwfr (nrmt(is), 1, io1, l1, ias) * apwdfr &
                 & (io2, l1, ias)
                  Call zmatinp (tapp, ngp, zt1, apwalm(:, io1, lm1, &
                 & ias), apwalm(:, io2, lm1, ias), v, h)
               End Do
            End Do
         End Do
      End Do
      Return
End Subroutine
!EOC
