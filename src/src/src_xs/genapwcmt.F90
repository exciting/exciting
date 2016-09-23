!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genapwcmt (lmax, ngp, isti, istf, apwalm, evecfv, wfcmt)
      Use modmain
      Use modinput
      Implicit None
  ! arguments
      Integer, Intent (In) :: lmax, isti, istf
      Integer, Intent (In) :: ngp
      Complex (8), Intent (In) :: apwalm (ngkmax, apwordmax, lmmaxapw, &
     & natmtot)
      Complex (8), Intent (In) :: evecfv (nmatmax, nstfv, nspnfv)
      Complex (8), Intent (Out) :: wfcmt (istf-isti+1, apwordmax, &
     & (lmax+1)**2, natmtot)
  ! local variables
      Integer :: ist, istc, is, ia, ias
      If (lmax .Gt. input%groundstate%lmaxapw) Then
         Write (*,*)
         Write (*, '("Error(genapwcmt): lmax > lmaxapw : ", I8)') lmax
         Write (*,*)
         Stop
      End If
      Do istc = isti, istf
         ist = istc - isti + 1
         Do is = 1, nspecies
            Do ia = 1, natoms (is)
               ias = idxas (ia, is)
               Call genapwcmt_part (lmax, ngp, ia, is, apwalm, &
              & evecfv(:, istc, 1), wfcmt(ist, :, :, ias))
            End Do
         End Do
      End Do
Contains
!
!
      Subroutine genapwcmt_part (lmax, ngp, ia, is, apwalm, evecfv, &
     & fcmt)
         Use modmain
         Implicit None
    ! arguments
         Integer, Intent (In) :: lmax, ia, is
         Integer, Intent (In) :: ngp
         Complex (8), Intent (In) :: apwalm (:, :, :, :)
         Complex (8), Intent (In) :: evecfv (:)
         Complex (8), Intent (Out) :: fcmt (:, :)
    ! external functions
         Complex (8) zdotu
         External zdotu
    ! local variables
         Integer :: ias, l, m, lm, io
         ias = idxas (ia, is)
         Do l = 0, lmax
            Do m = - l, l
               lm = idxlm (l, m)
               Do io = 1, apword (l, is)
                  fcmt (io, lm) = zdotu (ngp, evecfv, 1, apwalm(1, io, &
                 & lm, ias), 1)
               End Do
            End Do
         End Do
      End Subroutine genapwcmt_part
End Subroutine genapwcmt
