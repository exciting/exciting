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
  use exciting_mpi, only: mpiinfo, xmpi_allgather, xmpi_allreduce
  use modmpi, only: mpiglobal
  use normalize, only: normalize_vectors
  use precision, only: dp, i32
  use projection, only: project_y_onto_x
  use rttddft_Overlap, only: overlap_set
  use to_char_conversion, only: to_char
  use xlapack, only: hermitian_matrix_multiply, matrix_multiply

  implicit none

  private

  public :: initialize_wavefunction_set, obtain_occupations, obtain_projection_coefficients
  
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
      procedure :: obtain_number_excitations => wavefunction_set_obtain_number_excitations
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

    integer(i32) :: n_active_states
    integer(i32), allocatable :: buffer(:)

    call assert( n_frozen_ <= size( complete_gnd_set_lapwlo, 2 ), 'n_frozen_ > n_states')
    call assert( size( complete_gnd_set_lapwlo, 2 ) == size( occupations, 1 ), &
      'complete_gnd_set_lapwlo and occupations have different n_states')
    call assert( size( complete_gnd_set_lapwlo, 3 ) == size( occupations, 2 ), &
      'complete_gnd_set_lapwlo and occupations have different n_kpts')

    allocate( this%groundstate, source = complete_gnd_set_lapwlo )
    n_active_states = last_occupied_for_current_rank( occupations, occs_tol )
    ! Force the same number of active states over all MPI ranks
    call xmpi_allgather( mpiglobal, n_active_states, buffer )
    n_active_states = maxval( buffer )
    allocate( this%active, source = complete_gnd_set_lapwlo(:, n_frozen_ + 1 : n_active_states, :) )
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

    integer(i32) :: n_gnd_states, i, n_active_states
    integer(i32), allocatable :: buffer(:)

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

    n_active_states = last_occupied_for_current_rank( occupations, occs_tol )
    ! Force the same number of active states over all MPI ranks
    call xmpi_allgather( mpiglobal, n_active_states, buffer )
    n_active_states = maxval( buffer )
    allocate( this%active, source = this%groundstate(:, n_frozen_ + 1 : n_active_states, :) )
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

  !> Normalize the wavefunctions \(|\Psi_{i\mathbf{k}}\rangle\).   
  !> It is essentially a wrapper to the subroutine [[normalize_vectors]]
  subroutine normalize_wavefunctions( this, S, normalize_all )
    class(wavefunction_set), intent(inout) :: this
    !> Object that encapsulates the overlap matrices
    class(overlap_set), intent(in) :: S
    !> If `.true.`, frozen, save, and ground components are also normalized
    logical, optional, intent(in) :: normalize_all

    logical :: local_normalize_all

    local_normalize_all = .false.
    if ( present( normalize_all ) ) local_normalize_all = normalize_all

    call wrapper_normalize_vectors( S%array, this%active )
    if ( local_normalize_all ) then
      call wrapper_normalize_vectors( S%array, this%groundstate )
      if ( this%has_frozen() ) call wrapper_normalize_vectors( S%array, this%frozen )
      if ( allocated( this%active_save ) ) call wrapper_normalize_vectors( S%array, this%active_save )
    end if
    contains 
      subroutine wrapper_normalize_vectors(overlap, vectors)
        complex(dp), contiguous, intent(in) :: overlap(:, :, :)
        complex(dp), contiguous, intent(inout) :: vectors(:, :, :)
        
        integer(i32) :: ik

        call assert( size( overlap, 3 ) == size( vectors, 3 ), "Incompatible size" )
        do ik = 1, size( vectors, 3 )
          call normalize_vectors( S=overlap(:, :, ik), vectors=vectors(:, :, ik) )
        end do
      end subroutine
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

    if( allocated( this%active_save ) ) call assert( size( this%active, 3 ) == size( this%active_save, 3 ), &
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
  subroutine wavefunction_set_obtain_number_excitations( this, overlap, eps_occ, occ_gnd, wkpt, mpi_env, &
    & n_exc, n_gs )
    !> Basis-expansion coefficients of the KS-wavefunctions.
    class(wavefunction_set), intent(in) :: this
    !> Overlap matrices
    class(overlap_set), intent(in) :: overlap
    !> Occupation threshold above which a state is considered occupied
    real(dp), intent(in) :: eps_occ
    !> List of occupations at \(t=0\)
    real(dp), contiguous, intent(in) :: occ_gnd(:, :)
    !> k-point integration weights
    real(dp), contiguous, intent(in) :: wkpt(:)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> number of excited electrons
    real(dp), intent(out) :: n_exc
    !> number of electrons on the groundstate state
    real(dp), intent(out) :: n_gs

    integer(i32) :: ik, n_kpt
    real(dp) :: buffer(2)
    real(dp), allocatable :: occ(:, :), aux_tot(:), aux_exc(:)
    complex(dp), allocatable :: proj(:, :, :)
    
    n_kpt = this%n_kpts()
    call assert( size( wkpt ) == n_kpt, 'wkpt must have n_kpt elements')

    allocate( aux_tot(n_kpt), aux_exc(n_kpt) )
    call obtain_projection_coefficients( this%groundstate, overlap%array, this%active, proj )
    call obtain_occupations( proj, occ_gnd(this%first_active(): this%n_occupied(), :), occ )
    if ( this%has_frozen() ) occ(1 : this%n_frozen(), :) = occ(1 : this%n_frozen(), :) + occ_gnd(1 : this%n_frozen(), :)
    do concurrent (ik = 1:n_kpt)
      aux_tot(ik) = sum( occ(:, ik) )
      aux_exc(ik) = sum( occ(:, ik), occ_gnd(:, ik) <= eps_occ )
    end do
    n_exc = dot_product( wkpt, aux_exc )
    n_gs = dot_product( wkpt, aux_tot ) - n_exc
    buffer = [ n_exc, n_gs ]
    call xmpi_allreduce( buffer, mpi_env )
    n_exc = buffer(1); n_gs = buffer(2)

  end subroutine wavefunction_set_obtain_number_excitations

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

