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
  use normalize, only: normalize_vectors
  use precision, only: dp, i32
  use projection, only: project_y_onto_x
  use to_char_conversion, only: to_char
  use xlapack, only: hermitian_matrix_multiply, matrix_multiply

  implicit none

  private

  public :: normalize_wavefunctions, obtain_occupations, obtain_projection_coefficients, Update_basis_derivative

  

contains

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
        B_now(:,:, ik) = B_now(:,:, ik) + &
          & atoms_velocities(1,ias)*mathcal_B(:,:,1,ias, ik) + &
          & atoms_velocities(2,ias)*mathcal_B(:,:,2,ias, ik) + &
          & atoms_velocities(3,ias)*mathcal_B(:,:,3,ias, ik)
      end do
    end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
  end subroutine

  !> Normalize the wavefunctions \(|\Psi_{i\mathbf{k}}\rangle\)
  !> It is essentially a wrapper to the subroutine [[normalize_vectors]]
  subroutine normalize_wavefunctions( overlap_matrices, wavefunctions )
    !> List of overlap matrices. The 1st and 2nd indexes refers to the matrix elements,
    !> the 3rd index refers to the k-points
    complex(dp), intent(in)    :: overlap_matrices(:, :, :)
    !> List of wavefunctions \(|\Psi_{i\mathbf{k}}\rangle\). The 1st index refers to LAPW basis, 
    !> the 2nd index, to the number of states, and 
    !> the 3rd index, to the k-points
    complex(dp), intent(inout) :: wavefunctions(:, :, :)

    integer(i32) :: ik

    do ik = 1, size( wavefunctions, 3 )
      call normalize_vectors( S=overlap_matrices(:, :, ik), vectors=wavefunctions(:, :, ik) )
    end do
  end subroutine

  !> Project the wavefunctions `y` onto `x` and store the projection coefficients.   
  !> For each `k-point` (3rd dimension), the projection `p` is calculated as
  !> \[ p_k = x_k^\dagger S_k y_k \]
  subroutine obtain_projection_coefficients( x, S, y, proj_coeff )
    !> Wavefunctions onto which the projection is carried out
    complex(dp), contiguous, intent(in) :: x(:, :, :)
    !> Overlap matrix
    complex(dp), contiguous, intent(in) :: S(:, :, :)
    !> Wavefunctions to be projected
    complex(dp), contiguous, intent(in) :: y(:, :, :)
    !> Projection coefficients
    complex(dp), allocatable, intent(out) :: proj_coeff(:, :, :)

    integer(i32) :: ik
    complex(dp), allocatable :: aux(:, :)

    associate( mx => size( x, 2 ), my => size( y, 1 ), n => size( y, 2 ), k => size( y, 3 ))
      call assert( size(S, 3) == k, 'S and y must have same size along 3rd dim.')
      call assert( size(x, 3) == k, 'x and y must have same size along 3rd dim.')

      allocate( aux(my, n) )
      allocate( proj_coeff(mx, n, k) )
      do ik = 1, k
        call project_y_onto_x( y(:, :, ik), x(:, :, ik), S(:, :, ik), proj_coeff(:, :, ik), aux )
      end do
    end associate
  end subroutine


  !> Obtain the occupation factors given the projections onto a reference basis set
  !> Given the projection coefficients \(p_{ijk}\) of \(|\Psi_{jk}\rangle\) onto
  !> \(|\phi^0_{ik}\rangle\) as
  !> \[ |\Psi_{jk}\rangle = \sum_{i=1}^m p_{ijk} |\phi^0_{ik}\rangle, \quad j = 1, \ldots, n. \]
  !> The occupation factors \(f_{ik}\) are obtained as
  !> \[ f_{ik} = \sum_{j=1}^n f^0_{jk}|p_{ijk}|^2, \quad i = 1, \ldots, m, \]
  !> where \(f^0_{jk}\) are the original occupation factors of \(|\phi^0_{ik}\rangle\) usually taken for \(t=0\)
  subroutine obtain_occupations( proj, occ_gnd, occ )
    !> List of projection coefficients. Each set of projection coefficients is a rank-2 array
    complex(dp), contiguous, intent(in) :: proj(:, :, :)
    !> List of occupations at \(t=0\). Each set of occupations is a rank-1 array
    real(dp), contiguous, intent(in) :: occ_gnd(:, :)
    !> List of new occupations. Each set of occupations is a rank-1 array
    real(dp), allocatable, intent(out) :: occ(:, :)

    integer(i32) :: ik
    real(dp), parameter :: tol = 1.e-8_dp
    
    associate( m => size(proj, 1), n => size(proj, 2), dim_k => size(proj, 3))
      call assert( size( occ_gnd, 2) == dim_k , 'occ_gnd and proj must have compatible dimensions' )
      call assert( size( occ_gnd, 1) == m , 'occ_gnd and proj must have compatible dimensions' )
      call assert( n <= m , 'n must be <= m' )
      do ik = 1, dim_k
        ! \sum_{i=1}^m |p_{ijk}|^2 must be <= 1 (is equal to 1 only if the basis |\phi^0_{ik}\rangle is complete)
        call assert( maxval( sum(abs(proj(:, :, ik))**2, dim=1) ) <= 1._dp + tol , &
          'proj cannot represent projection factors along ik = ' // to_char(ik) )
      end do

      allocate( occ(m, dim_k) )
      
      do ik = 1, dim_k
        call matrix_multiply( abs(proj(:, :, ik))**2, occ_gnd(1:n, ik), occ(:, ik) )
      end do
    end associate
  end subroutine
end module rttddft_Wavefunction
