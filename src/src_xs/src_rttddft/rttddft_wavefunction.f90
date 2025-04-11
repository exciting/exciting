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

  public :: initialize_wavefunction_set, obtain_occupations, obtain_projection_coefficients, update_basis_derivative
  
  !> Type for the set of wavefunctions expanded on a basis set
  type, public, abstract :: wavefunction_set
    !> Basis expansion coefficients of the frozen states
    complex(dp), allocatable :: frozen(:, :, :)
    !> Initial basis expansion coefficients of all states
    complex(dp), allocatable :: groundstate(:, :, :)
    !> Basis expansion coefficients of the active states at time \(t\)
    complex(dp), allocatable :: active(:, :, :)
    !> Storage for the expansion coefficients of the active states
    complex(dp), allocatable :: active_save(:, :, :)

    contains
      procedure(initialize_), public, deferred :: initialize
      procedure :: normalize => normalize_wavefunctions
      procedure :: has_frozen
      procedure :: save => save_wavefunctions
      procedure :: restore => restore_wavefunctions
      procedure :: n_kpts
      procedure :: n_active
      procedure :: first_active
      procedure :: n_frozen
      procedure :: n_basis
      procedure :: n_empty
      procedure :: n_occupied
      procedure :: expanded_in_lapwlo
  end type

  !> Type to encapsulate the wavefunction set expanded in the LAPW+lo basis
  type, extends (wavefunction_set) :: wavefunction_set_lapwlo_basis
  contains
    procedure :: initialize => initialize_set_in_lapwlo_basis
  end type

  !> Type to encapsulate the wavefunction set expanded in the KS basis
  type, extends (wavefunction_set) :: wavefunction_set_ks_basis
  contains
    procedure :: initialize => initialize_set_in_ks_basis
  end type

  abstract interface
    subroutine initialize_( this, save_needed, n_frozen_, complete_gnd_set_lapwlo, &
        occupations, occs_tol )
      import :: wavefunction_set, i32, dp
      class(wavefunction_set), intent(inout) :: this
      !> If `.true.`, active_save component should be allocated
      logical, intent(in) :: save_needed
      !> Number of the frozen states
      integer(i32), intent(in) :: n_frozen_
      !> Initial wavefunction set in the LAPW+lo basis, (n_basis, n_states, n_kpt)
      complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, :)
      !> State occupations array (n_states, n_kpt)
      real(dp), contiguous, intent(in) :: occupations(:, :)
      !> Minimal value of occupation for the state to be 'occupied'
      real(dp), intent(in) :: occs_tol
    end subroutine
  end interface

contains

  !> Initializes the wavefunction set class ([[wavefunction_set]])
  !> with a concrete type, depending on the basis set
  subroutine initialize_wavefunction_set( psi, use_lapwlo_basis, save_needed, n_frozen_, &
      complete_gnd_set_lapwlo, occupations, occs_tol )
    class(wavefunction_set), allocatable, intent(out) :: psi
    !> Whether the LAPW+lo basis should be used
    logical, intent(in) :: use_lapwlo_basis
    !> If `.true.`, active_save component should be allocated
    logical, intent(in) :: save_needed
    !> Number of the frozen states
    integer(i32), intent(in) :: n_frozen_
    !> Initial wavefunction set in the LAPW+lo basis, (n_basis, n_states, n_kpt)
    complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, :)
    !> State occupations array (n_states, n_kpt)
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Minimal value of occupation for the state to be 'occupied'
    real(dp), intent(in) :: occs_tol

    if ( use_lapwlo_basis ) then
      allocate( wavefunction_set_lapwlo_basis :: psi )
    else
      allocate( wavefunction_set_ks_basis :: psi )      
    end if
    call psi%initialize( save_needed, n_frozen_, complete_gnd_set_lapwlo, occupations, occs_tol )
  end subroutine

  !> Initialize the set from the ground state WFs expanded in LAWP+lo basis
  subroutine initialize_set_in_lapwlo_basis( this, save_needed, n_frozen_, complete_gnd_set_lapwlo, &
      occupations, occs_tol )
    class(wavefunction_set_lapwlo_basis), intent(inout) :: this
    !> If `.true.`, active_save component should be allocated
    logical, intent(in) :: save_needed
    !> Number of the frozen states
    integer(i32), intent(in) :: n_frozen_
    !> Initial wavefunction set in the LAPW+lo basis, (n_basis, n_states, n_kpt)
    complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, :)
    !> State occupations array (n_states, n_kpt)
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Minimal value of occupation for the state to be 'occupied'
    real(dp), intent(in) :: occs_tol

    call assert( n_frozen_ <= size( complete_gnd_set_lapwlo, 2 ), &
      'n_frozen_ > n_states')
    call assert( size( complete_gnd_set_lapwlo, 2 ) == size( occupations, 1 ), &
      'complete_gnd_set_lapwlo and occupations have different n_states')
    call assert( size( complete_gnd_set_lapwlo, 3 ) == size( occupations, 2 ), &
      'complete_gnd_set_lapwlo and occupations have different n_kpts')

    allocate( this%groundstate, source = complete_gnd_set_lapwlo )
    allocate( this%active, source = complete_gnd_set_lapwlo(:, n_frozen_ + 1 : &
      last_occupied_for_current_rank( occupations, occs_tol ), :) )
    if ( save_needed ) allocate( this%active_save, source = this%active )
    if ( n_frozen_ > 0 ) allocate( this%frozen, source = complete_gnd_set_lapwlo(:, 1 : n_frozen_, :) )

  end subroutine

  !> Initialize the set from the ground state WFs expanded in LAWP+lo basis
  subroutine initialize_set_in_ks_basis( this, save_needed, n_frozen_, complete_gnd_set_lapwlo, &
      occupations, occs_tol )
    class(wavefunction_set_ks_basis), intent(inout) :: this
    !> If `.true.`, active_save component should be allocated
    logical, intent(in) :: save_needed
    !> Number of the frozen states
    integer(i32), intent(in) :: n_frozen_
    !> Initial wavefunction set in the LAPW+lo basis, (n_basis, n_states, n_kpt)
    complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, :)
    !> State occupations array (n_states, n_kpt)
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Minimal value of occupation for the state to be 'occupied'
    real(dp), intent(in) :: occs_tol

    integer :: n_gnd_states, i

    n_gnd_states = size( complete_gnd_set_lapwlo, 2 )

    call assert( n_frozen_ <= n_gnd_states, 'n_frozen_ > n_states')
    call assert( n_gnd_states == size( occupations, 1 ), &
      'complete_gnd_set_lapwlo and occupations have different n_states')
    call assert( size( complete_gnd_set_lapwlo, 3 ) == size( occupations, 2 ), &
      'complete_gnd_set_lapwlo and occupations have different n_kpts')

    allocate( this%groundstate( n_gnd_states, n_gnd_states, &
      size( complete_gnd_set_lapwlo, 3 ) ), source = zzero )
    do i = 1, n_gnd_states
      this%groundstate(i, i, :) = zone
    end do

    allocate( this%active, source = this%groundstate(:, n_frozen_ + 1 : &
      last_occupied_for_current_rank( occupations, occs_tol ), :) )
    if ( save_needed ) allocate( this%active_save, source = this%active )
    if ( n_frozen_ > 0 ) allocate( this%frozen, source = this%groundstate(:, 1 : n_frozen_, :) )

  end subroutine

  !> Save current active component into the save component
  pure subroutine save_wavefunctions( this )
    class(wavefunction_set), intent(inout) :: this

    this%active_save = this%active
  end subroutine

  !> Restore active component from the save component
  pure subroutine restore_wavefunctions( this )
    class(wavefunction_set), intent(inout) :: this

    this%active = this%active_save
  end subroutine

  !> Normalize the wavefunctions \(|\Psi_{i\mathbf{k}}\rangle\)
  !> It is essentially a wrapper to the subroutine [[normalize_vectors]]
  subroutine normalize_wavefunctions( this, overlap_matrices, normalize_all )
    class(wavefunction_set), intent(inout) :: this
    !> Overlap matrices of the basis functions
    complex(dp), intent(in) :: overlap_matrices(:, :, :)
    !> If `.true.`, frozen, save, and ground components are also normalized
    logical, optional, intent(in) :: normalize_all

    integer(i32) :: ik, nkpts
    logical :: normalize_all_

    normalize_all_ = .false.
    if ( present( normalize_all ) ) normalize_all_ = normalize_all
    nkpts = this%n_kpts()

    do ik = 1, nkpts
      call normalize_vectors( S=overlap_matrices(:, :, ik), vectors=this%active(:, :, ik) )
    end do

    if ( normalize_all_ ) then
      
      do ik = 1, nkpts
        call normalize_vectors( S=overlap_matrices(:, :, ik), vectors=this%groundstate(:, :, ik) )
      end do

      if ( this%has_frozen() ) then
        do ik = 1, nkpts
          call normalize_vectors( S=overlap_matrices(:, :, ik), vectors=this%frozen(:, :, ik) )
        end do
      end if

      if ( allocated( this%active_save ) ) then
        do ik = 1, nkpts
          call normalize_vectors( S=overlap_matrices(:, :, ik), vectors=this%active_save(:, :, ik) )
        end do
      end if

    end if
  end subroutine

  !> Tells whether the set is expanded over the LAPW+lo basis
  logical function expanded_in_lapwlo( this )
    class(wavefunction_set), intent(in) :: this
    
    select type( this )
    type is( wavefunction_set_lapwlo_basis )
    expanded_in_lapwlo = .true.
    type is( wavefunction_set_ks_basis )
    expanded_in_lapwlo = .false.
    class default
      call assert( .false., 'unrecognized type passed to expanded_in_lapwlo' )
    end select
  end function


  !> Tells whether there are frozen states
  pure logical function has_frozen( this )
    class(wavefunction_set), intent(in) :: this

    has_frozen = allocated( this%frozen )
  end function

  !> Returns the number of frozen states
  pure integer function n_frozen( this )
    class(wavefunction_set), intent(in) :: this

    n_frozen = 0
    if ( this%has_frozen() ) n_frozen = size( this%frozen, 2 )
  end function

  !> Returns the number of empty states
  pure integer function n_empty( this )
    class(wavefunction_set), intent(in) :: this
    
    n_empty = size( this%groundstate, 2 ) - this%n_occupied()
  end function

  !> Returns the number of occupied states
  pure integer function n_occupied( this )
    class(wavefunction_set), intent(in) :: this
    
    n_occupied = this%n_frozen() + this%n_active()
  end function

  !> Returns the number of active states
  pure integer function n_active( this )
    class(wavefunction_set), intent(in) :: this

    n_active = size( this%active, 2 )
  end function

  !> Checks the k dimensions and returns the number of k points
  integer function n_kpts( this )
    class(wavefunction_set), intent(in) :: this

    call assert( size( this%active, 3 ) == size( this%groundstate, 3 ), &
    "active and groundstate must have the same number of elements along 3rd dim." )

    call assert( size( this%active, 3 ) == size( this%active_save, 3 ), &
    "active and active_save must have the same number of elements along 3rd dim." )

    if ( this%has_frozen() ) call assert( size( this%active, 3 ) == size( this%frozen, 3 ), &
    "active and frozen must have the same number of elements along 3rd dim." )

    n_kpts = size( this%active, 3 )
  end function

  !> Returns basis size
  pure integer function n_basis( this )
    class(wavefunction_set), intent(in) :: this

    n_basis = size( this%active, 1 )
  end function

  !> Returns the position of the first active state
  pure integer function first_active( this )
    class(wavefunction_set), intent(in) :: this
    
    first_active = this%n_frozen() + 1
  end function

  !> Update \(B_k\) as
  !> \[ B_\mathbf{k}(t) = \sum_J \dot{\mathbf{R}_J}\cdot 
  !> \mathcal{B}_{J\mathbf{k}}(t) \]
  !> where \(J\) indexes the atoms
  subroutine update_basis_derivative( atoms_velocities, mathcal_B, B_now, B_old )
    !> the velocities (in cartesian coordinates) of all atoms
    real(dp), intent(in) :: atoms_velocities(:, :)
    !> `mathcalB` measures how the ions displacements affect overlap elements
    !> \[ \mathcal{B}_{J\mu'\mu}^{\mathbf{k}} = \left \langle
    !> \phi_{\mu'}^{\mathbf{k}}\bigg| \frac{\partial}{\partial \mathbf{R}_J}
    !> \phi_{\mu}^{\mathbf{k}} \right\rangle \]
    complex(dp), intent(in) :: mathcal_B(:, :, :, :, :)
    !> on entry: \(B\) at time \(t-\Delta t\), on exit: \(B\) at time \(t\)
    complex(dp), intent(inout) :: B_now(:, :, :)
    !> on exit: \(B\) at time \(t-\Delta t\)
    complex(dp), intent(out) :: B_old(:, :, :)
    
    integer :: ias, ik, n_atoms, n_kpt

    call assert( size( atoms_velocities, 1 ) == 3, 'atoms_velocities must have size = 3 along dim = 1' )
    call assert( size( atoms_velocities, 2 ) == size( mathcal_B, 4 ), &
      'size(atoms_velocities,2) and size(mathcal_B,4) must be equal' )
    call assert( size( atoms_velocities, 2 ) == size( mathcal_B, 4 ), &
      'size(atoms_velocities,2) and size(mathcal_B,4) must be equal' )
    call assert( size( atoms_velocities, 2 ) == size( mathcal_B, 4 ), &
      'size(atoms_velocities,2) and size(mathcal_B,4) must be equal' )

    n_kpt = size( mathcal_B, 5)
    n_atoms = size( atoms_velocities, 2 )
    B_old = B_now
    B_now = zzero
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(ik,ias), &
    !$OMP& SHARED(n_kpt,n_atoms,B_now,atoms_velocities,mathcal_B)
    !$OMP DO
    do ik = 1, n_kpt
      do ias = 1, n_atoms
        B_now(:, :, ik) = B_now(:, :, ik) + &
          & atoms_velocities(1, ias) * mathcal_B(:, :, 1,ias, ik) + &
          & atoms_velocities(2, ias) * mathcal_B(:, :, 2,ias, ik) + &
          & atoms_velocities(3, ias) * mathcal_B(:, :, 3,ias, ik)
      end do
    end do
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
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
  !> where \(f^0_{jk}\) are the original occupation factors of \(|\Psi_{jk}\rangle\) usually taken for \(t=0\)
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
      call assert( size( occ_gnd, 1) == n , 'occ_gnd and proj must have compatible dimensions' )
      call assert( n <= m , 'n must be <= m' )
      do ik = 1, dim_k
        ! \sum_{i=1}^m |p_{ijk}|^2 must be <= 1 (is equal to 1 only if the basis |\phi^0_{ik}\rangle is complete)
        call assert( maxval( sum(abs(proj(:, :, ik))**2, dim=1) ) <= 1._dp + tol , &
          'proj cannot represent projection factors along ik = ' // to_char(ik) )
      end do

      allocate( occ(m, dim_k) )
      
      do ik = 1, dim_k
        call matrix_multiply( abs(proj(:, :, ik))**2, occ_gnd(:, ik), occ(:, ik) )
      end do
    end associate
  end subroutine

  !> Returns the index of the last occupied state for the current rank
  pure integer function last_occupied_for_current_rank( occupations, occs_tol )
    !> State occupations array (n_states, n_kpt)
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Minimal value of occupation for the state to be 'occupied'
    real(dp), intent(in) :: occs_tol

    integer(i32) :: ik, i

    last_occupied_for_current_rank = -1
    do ik = 1, size( occupations, 2 )
      do i = size( occupations, 1 ), 1, -1
        if ( occupations(i, ik) > occs_tol ) exit
      end do
      if ( i > last_occupied_for_current_rank ) last_occupied_for_current_rank = i
    end do

  end function
end module rttddft_Wavefunction
