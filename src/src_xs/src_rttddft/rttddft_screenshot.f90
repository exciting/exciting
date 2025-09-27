! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! History
! Created by Ronaldo Rodrigues Pela, July 2019
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module to manage "screenshots" of desired properties during a RT-TDDFT
!> propagation
module rttddft_screenshot
  use asserts, only: assert
  use constants, only: real_zero, zzero
  use exciting_mpi, only: mpiinfo, xmpi_gatherv
  use precision, only: dp, i32
  use rttddft_Hamiltonian, only: hamiltonian_set
  use rttddft_input, only: screenshot_keys
  use rttddft_io, only: out_dens => write_density_to_file, &
                        out_eigs => write_eigenvalues, &
                        out_occs => write_occupations, &
                        out_proj => write_projection_coefficients
  use rttddft_Overlap, only: overlap_set
  use rttddft_Wavefunction, only: obtain_occupations, obtain_projection_coefficients, wavefunction_set
  use xlapack, only: solve_generalized_hermitian_eigenproblem


  implicit none

  private

  public :: screenshot

contains

  !> subroutine that takes "screenshots" during a RT-TDDFT propagation
  subroutine screenshot( it, input_keys, overlap, psi, H, &
      rho_MT, rho_interstitial, rho_MT_0, rho_interstitial_0, mpi_env )
    !> number of the current iteration (to name output files)
    integer(i32), intent(in) :: it
    !> Type that encapsulates the elements/attributes defined inside `screenshots` (in the input file)
    type(screenshot_keys), intent(in) :: input_keys
    !> overlap matrix
    class(overlap_set), intent(in) :: overlap
    !> Basis-expansion coefficients of the KS-wavefunctions.
    class(wavefunction_set), intent(in) :: psi
    !> Objtect that encapsulates the Hamiltonian matrix
    class(hamiltonian_set), intent(in) :: H
    !> electron density inside MT spheres
    real(dp), contiguous, intent(in) :: rho_MT(:, :, :)
    !> electron density in the interstitial region
    real(dp), contiguous, intent(in) :: rho_interstitial(:)
    !> electron density inside MT spheres at \( t=0 \)
    real(dp), contiguous, optional, intent(in) :: rho_MT_0(:, :, :)
    !> electron density in the interstitial region at \( t=0 \)
    real(dp), contiguous, optional, intent(in) :: rho_interstitial_0(:)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: dim_k, m
    integer(i32), allocatable :: dimensions_buffer(:)
    real(dp), allocatable :: eigenvalues(:, :), occupations(:, :), buffer(:, :)
    complex(dp), allocatable :: proj_time(:, :, :), proj_buffer(:, :, :)
    complex(dp), allocatable :: complete_filled_set(:, :, :)
    logical :: my_rank_writes

    my_rank_writes = mpi_env%is_root
    dim_k = size( H%H_t%array, 3 )

    associate( p => input_keys%projection_coefficients, occ => input_keys%occupations )
      if( p%on .or. occ%on ) then

        allocate( complete_filled_set( psi%n_basis(), psi%n_occupied(), psi%n_kpts() ) )
        complete_filled_set(:, psi%first_active(): psi%n_occupied(), :) = psi%active
        if ( psi%has_frozen() ) complete_filled_set(:, 1: psi%n_frozen() , :) = psi%frozen

        ! Project the current WFs onto the ground-state ones
        call obtain_projection_coefficients( psi%groundstate, overlap%array, complete_filled_set, proj_time )
        if( p%on ) then
          ! Send results to root rank, storing in the buffer
          call xmpi_gatherv( mpi_env, proj_time, proj_buffer )
          ! Write to output
          if( my_rank_writes ) call out_proj( it, p%print_absolute_value, p%output_format, proj_buffer )
        end if
        if( occ%on ) then
          call obtain_occupations( proj_time, psi%occupations(1 : psi%n_occupied(), : ), occupations )
          ! Send results to root rank, storing in the buffer
          call xmpi_gatherv( mpi_env, occupations, buffer )
          ! Write to output
          if( my_rank_writes ) call out_occs( it, buffer, &
            occ%output_text_format, occ%output_binary_format, occ%output_format )
        end if
      end if
    end associate

    if( input_keys%eigenvalues%on ) then
      associate( n_eigs => input_keys%eigenvalues%n_eigenvalues, tol => input_keys%eigenvalues%tol )
        m = merge( size(overlap%array, 1), n_eigs, n_eigs <= 0 )
        allocate( eigenvalues(m, dim_k), source=real_zero )
        call obtain_eigenvalues( H%H_t%array, overlap%array, H%dims, n_eigs, tol, eigenvalues )
        ! Send results to root rank
        call xmpi_gatherv( mpi_env, eigenvalues, buffer )
        if( n_eigs <= 0 ) then
          call xmpi_gatherv( mpi_env, H%dims, dimensions_buffer )
        else
          allocate( dimensions_buffer(size(buffer, 2)), source=n_eigs )
        end if
      end associate
      if( my_rank_writes ) call out_eigs( it, buffer, dimensions_buffer )
    end if

    if( input_keys%density%on ) then
      associate( plot3d => input_keys%density%plot3d, rho_I => rho_interstitial, rho_I0 => rho_interstitial_0 )
        if( it == 0 ) then
          call out_dens( it, rho_MT, rho_I, delta_rho=.false., plot3d=plot3d, my_rank_writes=my_rank_writes )
        else
          call assert( present(rho_MT_0) .and. present(rho_interstitial_0), "arguments must be present")
          call assert( all( shape(rho_MT) == shape(rho_MT_0) ), "rho_MT and rho_MT_0 must have same shape" )
          call assert( size(rho_I) == size(rho_I0), "rho_interstitial and rho_interstitial_0 must have same size" )
          call out_dens( it, rho_MT-rho_MT_0, rho_I-rho_I0, delta_rho=.true., plot3d=plot3d, my_rank_writes=my_rank_writes )
        end if
      end associate
    end if

  end subroutine screenshot

  !> Obtain the eigenvalues of a list of matrices `H` by solving the problem
  !> \[ H_k X = \lambda S_k X \]
  !> where `S` is a list of overlap matrices
  subroutine obtain_eigenvalues( H, S, dimensions, n_eigenvalues, tol, eigenvalues )
    !> List of Hamiltonian matrices (each of them must be hermitian)
    complex(dp), contiguous, intent(in) :: H(:, :, :)
    !> List of overlap matrices (each of them must be positive definite)
    complex(dp), contiguous, intent(in) :: S(:, :, :)
    !> Dimension of each `H_k` and `S_k` in the list
    integer(i32), contiguous, intent(in) :: dimensions(:)
    !> Number of eigenvalues to be evaluated
    integer(i32), intent(in) :: n_eigenvalues
    !> Tolerance for solving the generalized eigenvalue problem
    real(dp), intent(in) :: tol 
    !> List of eigenvalues obtained
    real(dp), contiguous, intent(out) :: eigenvalues(:, :)

    integer(i32) :: ik, dim, m, n
    complex(dp), allocatable :: H_copy(:, :), S_copy(:, :)

    m = size( H, 1 )
    do ik = 1, size( H, 3 )
      dim = dimensions(ik)
      n = merge( dim, n_eigenvalues, n_eigenvalues <= 0 )
      allocate( S_copy, source=S(1:dim, 1:dim, ik) )
      allocate( H_copy, source=H(1:dim, 1:dim, ik) )
      call solve_generalized_hermitian_eigenproblem( H_copy, S_copy, tol, eigenvalues(1:n, ik) )
      deallocate( H_copy, S_copy )
    end do
  end subroutine

end module rttddft_screenshot
