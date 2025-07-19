!
!
!
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
module secular_equation
  implicit none

  private

  public :: seceqn

contains

! !REVISION HISTORY:
!   Created March 2004 (JKD)
!   Removed a call to arpack July 2022 (Andris)
!   Introduced an optional argument, changed to FORD documentation, Oct 2024 (Ronaldo)
!> Solve the first- and second-variational secular equations. See routines
!> `match`, `seceqnfv`, `seceqnss`, and `seceqnsv`.
Subroutine seceqn (ik, evalfv, evecfv, evecsv, cdft_maximum_overlap)
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
      !> If `.true.`, the maximum overlap method is employed within a constrained DFT calculation
      logical, optional, intent(in) :: cdft_maximum_overlap

  ! local variables
      Integer :: ispn!,ib
  ! time
      Real (8) :: ts0,ts1
      logical  :: is_maximum_overlap_method_used
!
  ! allocatable arrays
      Complex (8), Allocatable :: apwalm (:, :, :, :, :)

      
      Allocate (apwalm(ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      apwalm=zzero
      is_maximum_overlap_method_used = .false.
      if( present(cdft_maximum_overlap) ) is_maximum_overlap_method_used = cdft_maximum_overlap
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
         Call seceqnfv(ik, ispn, nmat(ispn,ik), ngk(ispn,ik), &
        &  igkig(:,ispn,ik), vgkc(:,:,ispn,ik), apwalm(:,:,:,:,ispn), sfacgk(:, :, ispn, ik), tpgkc(:, :, ispn, ik), &
        &  is_maximum_overlap_method_used, &
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
!EOC
end module
