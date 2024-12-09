!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! !ROUTINE: seceqn
!
!
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!   Removed a call to arpack July 2022 (Andris)
!   Adapted call to seceqnfv and changed to FORD documentation, Oct 2024 (Ronaldo)
Subroutine seceqn (ik, evalfv, evecfv, evecsv)
      Use modinput
      Use modmain
      Use modmpi

      !> k-point index
      Integer, Intent (In) :: ik
      !> first-variational eigenvalues
      Real (8), Intent (Out) :: evalfv (nstfv, nspnfv)
      !> first-variational eigenvectors
      Complex (8), Intent (Out) :: evecfv (nmatmax, nstfv, nspnfv)
      !> second-variational eigenvectors
      Complex (8), Intent (Out) :: evecsv (nstsv, nstsv)
  ! local variables
      Integer :: ispn!,ib
  ! time
      Real (8) :: ts0,ts1
!
  ! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)

      
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      apwalm=zzero
  ! loop over first-variational spins (nspnfv=2 for spin-spirals only)

  !
  !-IMPORTANT: the first-variational spinor index and the k-point index have been
  ! swapped in the following arrays: ngk, igkig, vgkl, vgkc, gkc, tpgkc, sfacgk
  !
      Do ispn = 1, nspnfv
         current_igkig => igkig(:,ispn,ik)
         current_vgkc => vgkc(:,:,ispn,ik)
     ! find the matching coefficients
         Call timesec(ts0)
         Call match (ngk(ispn, ik), gkc(:, ispn, ik), tpgkc(:, :, ispn, &
        & ik), sfacgk(:, :, ispn, ik), apwalm(:, :, :, :, ispn))
         Call timesec(ts1)
         timematch=ts1-ts0+timematch
     ! solve the first-variational secular equation
         Call seceqnfv(ik, nmat(ispn,ik), ngk(ispn,ik), &
        &  igkig(:,ispn,ik), vgkc(:,:,ispn,ik), apwalm(:,:,:,:,ispn), &
        & evalfv(:,ispn), evecfv(:,:,ispn))
      End Do
      If (isspinspiral()) Then
     ! solve the spin-spiral second-variational secular equation
         Call seceqnss (ik, apwalm, evalfv, evecfv, evecsv)
      Else
     ! solve the second-variational secular equation
        if (input%groundstate%modifiedSV) then
          Call seceqnsv2 (ik, apwalm, evalfv, evecfv, evecsv)
        else
          Call seceqnsv (ik, apwalm, evalfv, evecfv, evecsv)
        endif
      End If
!
      Deallocate (apwalm)
End Subroutine seceqn
