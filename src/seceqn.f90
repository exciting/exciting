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
      use precision, only: dp
      use ghost_band_filter, only: detect_ghost_bands, filter_ghost_bands, reshape_sv_arrays, restore_shape_sv_arrays, &
                                  & set_default_ghost_band_parameters
#ifdef USEOMP
      use omp_lib
#endif

      !> k-point index
      Integer, Intent (In) :: ik
      !> first-variational eigenvalues
      Real (dp), Intent (Out) :: evalfv (nstfv, nspnfv)
      !> first-variational eigenvectors
      Complex (dp), Intent (Out) :: evecfv (nmatmax, nstfv, nspnfv)
      !> second-variational eigenvectors
      Complex (dp), Intent (Out) :: evecsv (nstsv, nstsv)
      !> If `.true.`, the maximum overlap method is employed within a constrained DFT calculation
      logical, optional, intent(in) :: cdft_maximum_overlap
  ! local variables
      Integer :: ispn!,ib
      ! number of ghost states
      Integer :: n_ghost_states
      ! temporary second-variational eigenvectors expanded to match shape of first-variational eigenvectors
      Complex (dp) :: evecsv_temp(1, nstsv, nstsv)
      Real (dp) :: evalsv_temp (nstsv, 1)
  ! time
      Real (dp) :: ts0,ts1
      logical  :: is_maximum_overlap_method_used
!
  ! allocatable arrays
      Complex (dp), Allocatable :: apwalm (:, :, :, :, :)

#ifdef USEOMP
      ! It is really a bad idea from a programer point of view to execute this from more than
      ! one thread as the diagonalization is performed using multithreading. Moreover, in some
      ! cases it can lead to hanging executions.
      if (omp_in_parallel()) call terminate("Error(seceqn): cannot be called inside a parallel region (OpenMP).")
#endif
      
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

      call set_default_ghost_band_parameters

      call detect_ghost_bands(lorbe0, apwe0, nstfv, evalfv(1:nstfv, 1:nspnfv), &
                              & input%groundstate%GhostBands%toleranceSmallestAllowedEval, n_ghost_states)

      if ( (input%groundstate%GhostBands%filterGhostBands) .and. (n_ghost_states > 0) ) then
        call filter_ghost_bands(nstfv, n_ghost_states, evalfv(1:nstfv, 1:nspnfv), evecfv(1:nmatmax, 1:nstfv, 1:nspnfv))
      end if

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

      ! reshape arrays, so that they fit the input of detect_ghost_bands
      call reshape_sv_arrays(evalsv, evecsv, nstsv, ik, evalsv_temp, evecsv_temp)

      call detect_ghost_bands(lorbe0, apwe0, nstsv, evalsv_temp(1:nstsv, 1:1), &
                              & input%groundstate%GhostBands%toleranceSmallestAllowedEval, n_ghost_states)

      if ( (input%groundstate%GhostBands%filterGhostBands) .and. (n_ghost_states > 0) ) then
        call filter_ghost_bands(nstsv, n_ghost_states, evalsv_temp(1:nstsv, 1:1), evecsv_temp(1:1, 1:nstsv, 1:nstsv))
        ! restore the original shape of the arrays
        call restore_shape_sv_arrays(evalsv_temp, evecsv_temp, nstsv, ik, evalsv, evecsv)
      end if

      Deallocate (apwalm)
End Subroutine seceqn
!EOC
end module
