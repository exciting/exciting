! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created May 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module that contains the subroutines envolved in the update of KS WFs
module rttddft_Wavefunction
  use asserts, only: assert
  use constants, only: zone, zzero, zi
  use integration, only: rk4 => ODESolver_RungeKutta4thOrder
  use matrix_exp, only: &
      & exp_hermitian => exp_hermitian_matrix_times_vectors, & 
      & exp_general => exp_general_matrix_times_vectors, &
      & exphouston_propagator => exphouston_hermitian_matrix_times_vectors
  use mod_kpoint, only: nkpt
  use mod_eigensystem, only: nmat
  use mod_eigenvalue_occupancy, only: nstfv
  use modinput, only: input
  use modmpi
  use normalize, only: normalizeWF => normalize_vectors
  use precision, only: dp
  use rttddft_GlobalVariables, only: ham_time, ham_past, overlap, &
    & evecfv_time, tstep, method, B_time, B_past, mathcalB

  implicit none

  private

  public :: UpdateWavefunction, Update_basis_derivative

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
  subroutine UpdateWavefunction( predcorr, atoms_velocities )
    !> tells if we are in the loop of the predictor-Corrector scheme
    logical, intent(in)       :: predcorr
    !> if present, we need to add corrections due to the nuclei movement (which changes the basis set)
    real(dp), intent(in), optional  :: atoms_velocities(:, :)
    !> Counter for loops with k-points
    integer                   :: ik
    !> Dimension of the Hamiltonian and Overlap matrices (given a k-point)
    integer                   :: nmatp
    !> indexes of the first and the last k-points
    integer                   :: first_kpt, last_kpt
    ! Factors that multiply the hamiltonian in the following propagator:
    ! Commutator-Free Magnus expansion of 4th order
    real(dp), parameter       :: f1 =  0.21132486540518713_dp ! 1/2 - sqrt(3)/6
    real(dp), parameter       :: f2 =  0.78867513459481290_dp ! 1/2 + sqrt(3)/6
    real(dp), parameter       :: a1 = -0.03867513459481287_dp ! 1/4 - sqrt(3)/6
    real(dp), parameter       :: a2 =  0.53867513459481290_dp ! 1/4 + sqrt(3)/6f1, f2, a1, a2
    !> auxiliary variable to store the overlap and hamiltonian matrices
    complex(dp), allocatable  :: overl(:,:), ham(:,:), hamold(:,:)
    logical                   :: atoms_movement
    procedure(exp_general), pointer      :: exp_operator




    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)

    atoms_movement = present( atoms_velocities )
    if ( atoms_movement ) then
      call Update_basis_derivative( atoms_velocities, mathcalB, B_time, B_past )
      exp_operator => exp_general
      call assert( trim(method)/='RK4', 'atoms_movement does not support RK4')
      call assert( trim(method)/='EH' .and. trim(method)/='EHM', 'atoms_movement does not support EH and EHM')
      call assert( .not. predcorr, 'atoms_movement does not support predictor-corrector')
    else
      exp_operator => exp_hermitian
    end if

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik,nmatp,overl,ham,hamold), &
!$OMP& SHARED(first_kpt, last_kpt,nstfv,nkpt,method,predcorr), &
!$OMP& SHARED(input,nmat,tstep,ham_time,ham_past,overlap,evecfv_time) &
!$OMP& SHARED(B_time, B_past, atoms_movement,exp_operator)
!$OMP DO
#endif
    do ik = first_kpt, last_kpt
      ! Dimension of the Hamiltonian and Overlap matrices for the current k-point
      nmatp = nmat(1,ik)

      allocate(overl(nmatp,nmatp), source=overlap(1:nmatp,1:nmatp,ik))
      allocate(ham(nmatp,nmatp))

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
          if( atoms_movement ) ham = ham - zi*B_time(1:nmatp,1:nmatp,ik)
          call exp_operator( &
              order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              alpha=-zi*tstep, H=ham, S=overl, &
              vectors=evecfv_time(1:nmatp, :, ik) )
        case ('EMR')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = 1.5_dp*ham_time(1:nmatp,1:nmatp,ik) -0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            if( atoms_movement ) ham = ham - zi*( 1.5_dp*B_time(1:nmatp,1:nmatp,ik)-0.5_dp*B_past(1:nmatp,1:nmatp,ik) )
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik) +0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
          end if
          call exp_operator( &
            & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
            & alpha=-zi*tstep, H=ham, S=overl, &
            & vectors=evecfv_time(1:nmatp, :, ik) )
        case ('AETRS')
          if ( .not. predcorr ) then
            ! 1/2*H(t)
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik)
            if( atoms_movement ) ham = ham - zi*0.5_dp*B_time(1:nmatp,1:nmatp,ik)
            call exp_operator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ! extrapolated 1/2*H(t+\Delta t) as H(t) - 1/2*H(t-\Delta t)
            ham(1:nmatp,1:nmatp) = ham_time(1:nmatp,1:nmatp,ik)-0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            if( atoms_movement ) ham = ham - zi*( B_time(1:nmatp,1:nmatp,ik)-0.5_dp*B_past(1:nmatp,1:nmatp,ik) )
            call exp_operator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          else
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_past(1:nmatp,1:nmatp,ik)
            call exp_hermitian( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ham(1:nmatp,1:nmatp) = 0.5_dp*ham_time(1:nmatp,1:nmatp,ik)
            call exp_hermitian( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          end if
        case ('CFM4')
          if ( .not. predcorr ) then
            ham(1:nmatp,1:nmatp) = (a1*(1+f2) + a2*(1+f1))*ham_time(1:nmatp,1:nmatp,ik) &
              -(a1*f2 + a2*f1)*ham_past(1:nmatp,1:nmatp,ik)
            if( atoms_movement ) ham = ham - zi*( (a1*(1+f2) + a2*(1+f1))*B_time(1:nmatp,1:nmatp,ik)&
              -(a1*f2 + a2*f1)*B_past(1:nmatp,1:nmatp,ik) )
            call exp_operator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ham(1:nmatp,1:nmatp) = (a1*(1+f1) + a2*(1+f2))*ham_time(1:nmatp,1:nmatp,ik) &
              -(a1*f1 + a2*f2)*ham_past(1:nmatp,1:nmatp,ik)
            if( atoms_movement ) ham = ham - zi*( (a1*(1+f1) + a2*(1+f2))*B_time(1:nmatp,1:nmatp,ik)&
              -(a1*f1 + a2*f2)*B_past(1:nmatp,1:nmatp,ik) )
            call exp_operator( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep, &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
          else
            ham(1:nmatp,1:nmatp) = a1*((1-f2)*ham_past(1:nmatp,1:nmatp,ik)+f2*ham_time(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1-f1)*ham_past(1:nmatp,1:nmatp,ik)+f1*ham_time(1:nmatp,1:nmatp,ik))
            call exp_hermitian( &
              & order_taylor=input%xs%realTimeTDDFT%TaylorOrder, &
              & alpha=-zi*tstep,  &
              & H=ham, S=overl, vectors=evecfv_time(1:nmatp, :, ik) )
            ham(1:nmatp,1:nmatp) = a1*((1-f1)*ham_past(1:nmatp,1:nmatp,ik)+f1*ham_time(1:nmatp,1:nmatp,ik)) + &
                                   a2*((1-f2)*ham_past(1:nmatp,1:nmatp,ik)+f2*ham_time(1:nmatp,1:nmatp,ik))
            call exp_hermitian( &
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

  !> Update \(B_k\) as
  !> \[ B_\mathbf{k}(t) = \sum_J \dot{\mathbf{R}_J}\cdot 
  !> \mathcal{B}_{J\mathbf{k}}(t) \]
  !> where \(J\) indexes the atoms
  subroutine Update_basis_derivative( atoms_velocities, mathcal_B, B_now, B_old )
    !> the velocities (in cartesian coordinates) of all atoms
    real(dp), intent(in)       :: atoms_velocities(:, :)
    !> `mathcalB` measures how the ions displacements affect overlap elements
    !> \[ \mathcal{B}_{J\mu'\mu}^{\mathbf{k}} = \left \langle
    !> \phi_{\mu'}^{\mathbf{k}}\bigg| \frac{\partial}{\partial \mathbf{R}_J}
    !> \phi_{\mu}^{\mathbf{k}} \right\rangle \]
    complex(dp), intent(in)    :: mathcal_B(:, :, :, :, :)
    !> on entry: \(B\) at time \(t-\Delta t\), on exit: \(B\) at time \(t\)
    complex(dp), intent(inout) :: B_now(:, :, :)
    !> on exit: \(B\) at time \(t-\Delta t\)
    complex(dp), intent(out)   :: B_old(:, :, :)
    
    integer :: ias, ik, n_atoms, n_kpt

    call assert( size(atoms_velocities,1)==3, 'atoms_velocities must have size = 3 along dim = 1' )
    call assert( size(atoms_velocities,2)==size(mathcal_B,4), &
      'size(atoms_velocities,2) and size(mathcal_B,4) must be equal' )
    call assert( size(atoms_velocities,2)==size(mathcal_B,4), &
      'size(atoms_velocities,2) and size(mathcal_B,4) must be equal' )
    call assert( size(atoms_velocities,2)==size(mathcal_B,4), &
      'size(atoms_velocities,2) and size(mathcal_B,4) must be equal' )

    n_kpt = size( mathcal_B, 5)
    n_atoms = size( atoms_velocities, 2 )
    B_old = B_now
    B_now = zzero
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik,ias), &
!$OMP& SHARED(n_kpt,n_atoms,B_now,atoms_velocities,mathcal_B)
!$OMP DO
#endif
    do ik = 1, n_kpt
      do ias = 1, n_atoms
        B_now(:,:,ik) = B_now(:,:,ik) + &
          & atoms_velocities(1,ias)*mathcal_B(:,:,1,ias,ik) + &
          & atoms_velocities(2,ias)*mathcal_B(:,:,2,ias,ik) + &
          & atoms_velocities(3,ias)*mathcal_B(:,:,3,ias,ik)
      end do
    end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
  end subroutine

end module rttddft_Wavefunction
