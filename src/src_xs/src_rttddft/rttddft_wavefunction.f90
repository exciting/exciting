! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created May 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module that contains the subroutines envolved in the update of KS WFs
module rttddft_Wavefunction
  use rttddft_GlobalVariables, only: ham_time, ham_past, overlap, &
    & evecfv_time, tstep, method
  use mod_kpoint, only: nkpt
  use mod_eigensystem, only: nmat
  use mod_eigenvalue_occupancy, only: nstfv
  use modinput, only: input
  use constants, only: zone, zzero, zi
  use precision, only: dp
  use modmpi
  use asserts, only: assert

  implicit none

  private

  public :: UpdateWavefunction

contains
  !> This subroutine updates KS wavefunctions.  
  !> Here, we employ a propagator to evolve the Kohn-Sham wavefunctions.  
  !> Extrapolation scheme for the hamiltonian (predcorr .False.)
  !> \[ \hat{H}(t+f\Delta t) = (1+f)\hat{H}(t) - f\hat{H}(t-\Delta t). \]
  !> \(\hat{H}(t)\) is stored in `ham_time`, whereas \(\hat{H}(t -\Delta t)\),
  !> in `ham_past`
  !> Extrapolation scheme for the hamiltonian (predcorr .True.)
  !> \[ \hat{H}(t+f\Delta t) = f\hat{H}(t + \Delta t) + (1-f)\hat{H}(t). \]
  !> \(\hat{H}(t)\) is stored in `ham_past`, whereas \(\hat{H}(t +\Delta t)\),
  !> in `ham_time` (which comes from a previous iteration in the predictor corrector
  !> loop
  subroutine UpdateWavefunction( predcorr )
    use matrix_exp, only: &
      & exp_propagator => exp_hermitianoperator_times_wavefunctions, & 
      & exphouston_propagator => exphouston_hermitianoperator_times_wavefunctions
    use normalize, only: normalizeWF => normalize_vectors
    use integration, only: rk4 => ODESolver_RungeKutta4thOrder

    implicit none
    !> tells if we are in the loop of the predictor-Corrector scheme
    logical, intent(in)       :: predcorr
    !> Counter for loops with k-points
    integer                   :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer                   :: nmatp
    !> indexes of the first and the last k-points
    integer                   :: first_kpt, last_kpt
    !> Factors that multiply the hamiltonian in the following propagator:
    !> Commutator-Free Magnus expansion of 4th order
    real(dp)                  :: f1, f2, a1, a2
    !> auxiliary variable to store the overlap and hamiltonian matrices
    complex(dp), allocatable  :: overl(:,:), ham(:,:), hamold(:,:)

    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik,nmatp,overl,ham,hamold), &
!$OMP& SHARED(first_kpt, last_kpt,a1,a2,f1,f2,nstfv,nkpt,method,predcorr), &
!$OMP& SHARED(input,nmat,tstep,ham_time,ham_past,overlap,evecfv_time)
!$OMP DO
#endif
    do ik = first_kpt, last_kpt
      ! Dimension of the Hamiltonian and Overlap matrices for the current k-point
      nmatp = nmat(1,ik)

      allocate(overl(nmatp,nmatp))
      allocate(ham(nmatp,nmatp))
      overl(1:nmatp,1:nmatp) = overlap(1:nmatp,1:nmatp,ik)

      select case(method)
        ! SE (simple exponential)
        ! CN (Crank-Nicolson)
        ! EMR (Exponential at midpoint rule)
        ! AETRS (approximate enforced time-reversal symmetry)
        ! CFM4 (Commutator-Free Magnus expansion of 4th order)
        ! EH (exponential using a basis of the hamiltonian-eigenvectors)
        ! EHM (same as before, but uses the hamiltonian at midpoint)
        ! RK4 (Runge-Kutta of 4th order)
        case ('SE')
          ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)
          call exp_propagator( &
            & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
            & alpha=-zi*tstep, H=ham, S=overl, &
            & vectors=evecfv_time(1:nmatp, :, ik) )
        case ('EMR')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = 1.5_dp*ham_time(1:nmatp,1:nmatp,ik) -0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik) +0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          end if
          call exp_propagator( &
            & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
            & alpha=-zi*tstep, H=ham, S=overl, &
            & vectors=evecfv_time(1:nmatp, :, ik) )
        case ('AETRS')
          if ( .not. predcorr ) then
            ! 1/2*H(t)
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik)
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ! extrapolated 1/2*H(t+\Delta t) as H(t) - 1/2*H(t-\Delta t)
            ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)-0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik)
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          end if
        case ('CFM4')
          f1 =  0.21132486540518713_dp ! 1/2 - sqrt(3)/6
          f2 =  0.78867513459481290_dp ! 1/2 + sqrt(3)/6
          a1 = -0.03867513459481287_dp ! 1/4 - sqrt(3)/6
          a2 =  0.53867513459481290_dp ! 1/4 + sqrt(3)/6
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = a1*((1+f2)*ham_time(1:nmatp,1:nmatp,ik)-f2*ham_past(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1+f1)*ham_time(1:nmatp,1:nmatp,ik)-f1*ham_past(1:nmatp,1:nmatp,ik))
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ham(1:nmatp,1:nmatp) = a1*((1+f1)*ham_time(1:nmatp,1:nmatp,ik)-f1*ham_past(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1+f2)*ham_time(1:nmatp,1:nmatp,ik)-f2*ham_past(1:nmatp,1:nmatp,ik))
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          else
            ham(1:nmatp,1:nmatp) = a1*((1-f2)*ham_past(1:nmatp,1:nmatp,ik)+f2*ham_time(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1-f1)*ham_past(1:nmatp,1:nmatp,ik)+f1*ham_time(1:nmatp,1:nmatp,ik))
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep,  &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ham(1:nmatp,1:nmatp) = a1*((1-f1)*ham_past(1:nmatp,1:nmatp,ik)+f1*ham_time(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1-f2)*ham_past(1:nmatp,1:nmatp,ik)+f2*ham_time(1:nmatp,1:nmatp,ik))
            call exp_propagator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep,  &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          end if
        case ('EH')
          ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)
          call exphouston_propagator( alpha=-zi*tstep, &
            & H=ham, S=overl, &
            & vectors=evecfv_time(1:nmatp, :, ik), &
            & tol=input%groundstate%solver%evaltol)
        case ('EHM')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = 1.5_dp*ham_time(1:nmatp,1:nmatp,ik) -0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik) +0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          end if
          call exphouston_propagator( alpha=-zi*tstep, H=ham, S=overl, &
            & vectors=evecfv_time(1:nmatp, :, ik), &
            & tol=input%groundstate%solver%evaltol)
        case ('RK4')
          ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)
          allocate(hamold(1:nmatp,1:nmatp))
          hamold(1:nmatp,1:nmatp) = ham_past(1:nmatp,1:nmatp,ik)
          if( .not. predcorr ) then
            call rk4( time_step=tstep, alpha=zi, &
              & H=ham, H_past=hamold, S=overl, &
              & x=evecfv_time(1:nmatp, :, ik))
          else
            ! Trick: H(t-dt) = 2*H(t)-H(t+dt), where H(t) = hamold, H(t+dt)=ham
            call rk4( time_step=tstep, alpha=zi, &
              & H=hamold, H_past=2_dp*hamold-ham, S=overl, &
              & x=evecfv_time(1:nmatp, :, ik))
          end if
          deallocate(hamold)
      end select

      ! Normalize WFs, if this is the case
      ! The propagator operator should be unitary, so this step would be unnecessary
      ! However, numerically this is almost never possible
      ! This normalization may help to avoid numerical issues
      if ( input%xs%realTimeTDDFT%normalizeWF ) then
        call normalizeWF( S=overl, &
          & vectors=evecfv_time(1:nmatp, :, ik) )
      end if

      deallocate(overl)
      deallocate(ham)
    end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

  end subroutine UpdateWavefunction

end module rttddft_Wavefunction
