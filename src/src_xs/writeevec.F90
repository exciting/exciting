!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writeevec (vq, voff, filxt)
      Use modmain
      Use modinput
      Use modmpi
      Use modxs
      Use m_gndstateq
      Use m_filedel
      Implicit None
  ! arguments
      Real (8), Intent (In) :: vq (3), voff (3)
      Character (*), Intent (In) :: filxt
  ! local variables
      Integer :: ik, j
      Complex (8), Allocatable :: apwalm (:, :, :, :)
  ! read from STATE.OUT exclusively
      isreadstate0 = .True.
  ! SCF calculation with one cycle
      Call gndstateq (voff, filxt)
      If (allocated(evecfv)) deallocate (evecfv)
      Allocate (evecfv(nmatmax, nstfv, nspnfv))
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (apwcmt(nstfv, apwordmax, lmmaxapw, natmtot))
      Allocate (locmt(nstfv, nlomax,-lolmax:lolmax, natmtot))
  ! delete existing coefficients files
      If (rank .Eq. 0) Call filedel ('APWCMT'//trim(filxt))
      If (rank .Eq. 0) Call filedel ('LOCMT'//trim(filxt))
      Call genparidxran ('k', nkpt)
      Do ik = kpari, kparf
         apwcmt (:, :, :, :) = zzero
         locmt (:, :, :, :) = zzero
         Call getevecfv (vkl(1, ik), vgkl(1, 1, 1, ik), evecfv)
         Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
        & sfacgk(1, 1, 1, ik), apwalm)
         Call genapwcmt (input%groundstate%lmaxapw, ngk(1, ik), 1, &
        & nstfv, apwalm, evecfv, apwcmt)
         Call genlocmt (ngk(1, ik), 1, nstfv, evecfv, locmt)
         Do j = 0, procs - 1
            If (rank .Eq. j) Then
               Call putapwcmt ('APWCMT'//trim(filxt), ik, vkl(1, ik), &
              & vq, apwcmt)
               Call putlocmt ('LOCMT'//trim(filxt), ik, vkl(1, ik), vq, &
              & locmt)
            End If
            Call barrier
         End Do
      End Do
      Call endloopbarrier (nkpt, procs)
      isreadstate0 = .False.
      Deallocate (evecfv, apwalm, apwcmt, locmt)
End Subroutine writeevec
