! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created: July 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Refactored: January 2025 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module to obtain the number of excitations. 
module rttddft_NumberExcitations
  use asserts, only: assert
  use exciting_mpi, only: mpiinfo, xmpi_allreduce
  use precision, only: dp, i32
  use rttddft_Wavefunction, only: obtain_occupations, obtain_projection_coefficients
  
  implicit none

  private

  public :: Obtain_number_excitations

contains

  !> Within this subroutine, we obtain the number of excitations, as described
  !> below.  
  !> The number of excited electrons after the interaction with a laser pulse
  !> In RT-TDDFT, the occupation number \( f_{j\mathbf{k}} \) of a KS state is
  !> kept fixed to its initial value. As the wavefunctions evolve, they are not
  !> any longer eigenstates of \( \hat{H}(t) \). It is possible to describe
  !> the number of excitations by projecting \( | \psi_{i\mathbf{k}}(t)\rangle \)
  !> onto the reference ground state at \( t=0 \).
  !> For a given k-point, we define the number of electrons that have
  !> been excited to an unoccupied KS state, labeled  \( j \), as
  !> \[
  !> 	m_{j\mathbf{k}}(t)= \sum_{i} f_{i\mathbf{k}}| \langle \psi_{j\mathbf{k}}(0)
  !>	         | \psi_{i\mathbf{k}}(t)\rangle |^2.
  !> \]
  !> Similarly, the number of holes created in an occupied KS \( j' \) state can
  !> specified as
  !> 	\[
  !> 	m_{j'\mathbf{k}}(t)= f_{j'\mathbf{k}} - \sum_{i}
  !> 	f_{i\mathbf{k}}	| \langle \psi_{j'\mathbf{k}}(0)| \psi_{i\mathbf{k}}(t)\rangle |^2.
  !> 	\]
  !> Thus, the total number of excited electrons in a unit cell can be
  !> obtained by considering all the unoccupied states
  !> \[
  !> 	N_{exc}(t)=
  !> 	\sum_{j\mathbf{k}}^{j\, unocc}
  !> 	w_\mathbf{k} m_{j\mathbf{k}}(t) = \sum_{j'\mathbf{k}}^{j'\, occ}
  !> 	w_\mathbf{k} m_{j'\mathbf{k}}(t) .
  !> 	\]
  subroutine Obtain_number_excitations( psi_gnd, psi, overlap, eps_occ, occ_gnd, wkpt, mpi_env, &
      & n_exc, n_gs )
    !> Basis-expansion coefficients of the KS-wavefunctions at \( t=0 \).
    complex(dp), contiguous, intent(in)   :: psi_gnd(:, :, :)
    !> Basis-expansion coefficients of the KS-wavefunctions at current time \(t\).
    complex(dp), contiguous, intent(in)   :: psi(:, :, :)
    !> Overlap matrices
    complex(dp), contiguous, intent(in)   :: overlap(:, :, :)
    !> Occupation threshold above which a state is considered occupied
    real(dp), intent(in) :: eps_occ
    !> List of occupations at \(t=0\)
    real(dp), contiguous, intent(in) :: occ_gnd(:, :)
    !> k-point integration weights
    real(dp), contiguous, intent(in) :: wkpt(:)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> number of excited electrons
    real(dp), intent(out)     :: n_exc
    !> number of electrons on the groundstate state
    real(dp), intent(out)     :: n_gs

    integer(i32) :: ik, n_kpt
    real(dp) :: buffer(2)
    real(dp), allocatable :: occ(:, :), aux_tot(:), aux_exc(:)
    complex(dp), allocatable  :: proj(:, :, :)
    
    n_kpt = size( psi, 3 )
    call assert( size(wkpt) == n_kpt, 'wkpt must have n_kpt elements')

    allocate( aux_tot(n_kpt), aux_exc(n_kpt) )
    call obtain_projection_coefficients( psi_gnd, overlap, psi, proj )
    call obtain_occupations( proj, occ_gnd, occ )
    do concurrent (ik = 1:n_kpt)
      aux_tot(ik) = sum( occ(:, ik) )
      aux_exc(ik) = sum( occ(:, ik), occ_gnd(:, ik) <= eps_occ )
    end do
    n_exc = dot_product(wkpt, aux_exc)
    n_gs = dot_product(wkpt, aux_tot) - n_exc
    buffer = [ n_exc, n_gs ]
    call xmpi_allreduce( buffer, mpi_env )
    n_exc = buffer(1); n_gs = buffer(2)

  end subroutine Obtain_number_excitations
end module rttddft_NumberExcitations
