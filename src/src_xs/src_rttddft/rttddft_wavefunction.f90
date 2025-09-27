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
  use mod_kpointset, only: k_set
  use modmpi, only: mpiglobal
  use normalize, only: normalize_vectors
  use precision, only: dp, i32
  use projection, only: project_y_onto_x
  use rttddft_file_formats, only: file_handler
  use rttddft_io_unformatted, only: read_wavefunction, t, t_minus_dt, write_wavefunction
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
    !> Initial occupations of KS states
    real(dp), allocatable :: occupations(:, :)
    !> Threshold above which a state is considered occupied
    real(dp) :: eps_occ
    !> Object that encapsulates the \( \mathbf{k} \)-points used by all MPI ranks
    type(k_set) :: kset

    contains
      procedure :: expanded_in_lapwlo => wavefunction_set_expanded_in_lapwlo
      procedure :: first_active => wavefunction_set_first_active
      procedure :: first_kpt => wavefunction_set_first_kpt
      procedure :: has_frozen => wavefunction_set_has_frozen
      procedure(initialize_), private, deferred :: initialize
      procedure :: last_kpt => wavefunction_set_last_kpt
      procedure :: n_active => wavefunction_set_n_active
      procedure :: n_basis => wavefunction_set_n_basis
      procedure :: n_empty => wavefunction_set_n_empty
      procedure :: n_frozen => wavefunction_set_n_frozen
      procedure :: n_kpts => wavefunction_set_n_kpts
      procedure :: n_occupied => wavefunction_set_n_occupied
      procedure :: normalize => normalize_wavefunctions
      procedure :: obtain_number_excitations => wavefunction_set_obtain_number_excitations
      procedure :: read_from_file => wavefunction_set_read
      procedure :: restore => restore_wavefunctions
      procedure :: save => save_wavefunctions
      procedure :: write_to_file => wavefunction_set_write
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
    subroutine initialize_( this, first_kpt, kset, save_needed, n_frozen, complete_gnd_set_lapwlo, &
        occupations, occs_tol )
      import :: dp, i32, k_set, wavefunction_set
      class(wavefunction_set), intent(inout) :: this
      !> Index of the first \( \mathbf{k} \)-point assigned to this MPI rank
      integer(i32), intent(in) :: first_kpt
      !> Object that encapsulates the \( \mathbf{k} \)-points used by all MPI ranks
      type(k_set), intent(in) :: kset
      !> If `.true.`, active_save component should be allocated
      logical, intent(in) :: save_needed
      !> Number of the frozen states
      integer(i32), intent(in) :: n_frozen
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
  subroutine initialize_wavefunction_set( psi, use_lapwlo_basis, first_kpt, kset, &
      save_needed, n_frozen, complete_gnd_set_lapwlo, occupations, occs_tol )
    class(wavefunction_set), allocatable, intent(out) :: psi
    !> Whether the LAPW+lo basis should be used
    logical, intent(in) :: use_lapwlo_basis
    !> Index of the first \( \mathbf{k} \)-point assigned to this MPI rank
    integer(i32), intent(in) :: first_kpt
    !> Object that encapsulates the \( \mathbf{k} \)-points used by all MPI ranks
    type(k_set), intent(in) :: kset
    !> If `.true.`, active_save component should be allocated
    logical, intent(in) :: save_needed
    !> Number of the frozen states
    integer(i32), intent(in) :: n_frozen
    !> Initial wavefunction set in the LAPW+lo basis, (n_basis, n_states, n_kpt)
    complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, first_kpt:)
    !> State occupations array (n_states, n_kpt)
    real(dp), contiguous, intent(in) :: occupations(:, first_kpt:)
    !> Minimal value of occupation for the state to be 'occupied'
    real(dp), intent(in) :: occs_tol

    if ( use_lapwlo_basis ) then
      allocate( wavefunction_set_lapwlo_basis :: psi )
    else
      allocate( wavefunction_set_ks_basis :: psi )      
    end if
    call psi%initialize( first_kpt, kset, save_needed, n_frozen, complete_gnd_set_lapwlo, occupations, occs_tol )
  end subroutine

  !> Initialize the set from the ground state WFs expanded in LAWP+lo basis. See [[initialize_]] for documentation.
  subroutine initialize_set_in_lapwlo_basis( this, first_kpt, kset, save_needed, n_frozen, &
      complete_gnd_set_lapwlo, occupations, occs_tol )
    class(wavefunction_set_lapwlo_basis), intent(inout) :: this
    integer(i32), intent(in) :: first_kpt
    type(k_set), intent(in) :: kset
    logical, intent(in) :: save_needed
    integer(i32), intent(in) :: n_frozen
    complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, first_kpt:)
    real(dp), contiguous, intent(in) :: occupations(:, first_kpt:)
    real(dp), intent(in) :: occs_tol

    integer(i32) :: last_kpt, n_active_states, n_basis

    last_kpt = ubound( complete_gnd_set_lapwlo, 3 )
    n_basis = size( complete_gnd_set_lapwlo, 1 )
    call assert( n_frozen <= size( complete_gnd_set_lapwlo, 2 ), 'n_frozen > n_states')
    call assert( size( complete_gnd_set_lapwlo, 2 ) == size( occupations, 1 ), &
      'complete_gnd_set_lapwlo and occupations have different n_states')
    call assert( size( complete_gnd_set_lapwlo, 3 ) == size( occupations, 2 ), &
      'complete_gnd_set_lapwlo and occupations have different n_kpts')

    this%eps_occ = occs_tol
    allocate( this%occupations, source = occupations )
    allocate( this%groundstate, source = complete_gnd_set_lapwlo )
    this%kset = kset
    ! TODO: n_active_states to be defined by the user (issue #248)
    n_active_states = last_occupied_for_all_ranks( occupations, occs_tol, mpiglobal )
    allocate( this%active(n_basis, n_active_states-n_frozen, first_kpt:last_kpt), &
      source = complete_gnd_set_lapwlo(:, n_frozen + 1 : n_active_states, :) )
    if ( save_needed ) allocate( this%active_save, source = this%active )
    if ( n_frozen > 0 ) allocate( this%frozen(n_basis, n_frozen, first_kpt:last_kpt), &
      source = complete_gnd_set_lapwlo(:, 1 : n_frozen, :) )
  end subroutine

  !> Initialize the set from the ground state WFs expanded in KS basis. See [[initialize_]] for documentation.
  subroutine initialize_set_in_ks_basis( this, first_kpt, kset, save_needed, n_frozen, complete_gnd_set_lapwlo, &
      occupations, occs_tol )
    class(wavefunction_set_ks_basis), intent(inout) :: this
    integer(i32), intent(in) :: first_kpt
    type(k_set), intent(in) :: kset
    logical, intent(in) :: save_needed
    integer(i32), intent(in) :: n_frozen
    complex(dp), contiguous, intent(in) :: complete_gnd_set_lapwlo(:, :, first_kpt:)
    real(dp), contiguous, intent(in) :: occupations(:, first_kpt:)
    real(dp), intent(in) :: occs_tol

    integer(i32) :: i, last_kpt, n_active_states, n_gnd_states

    n_gnd_states = size( complete_gnd_set_lapwlo, 2 )
    last_kpt = ubound( complete_gnd_set_lapwlo, 3 )

    call assert( n_frozen <= n_gnd_states, 'n_frozen > n_gnd_states')
    call assert( n_gnd_states == size( occupations, 1 ), &
      'complete_gnd_set_lapwlo and occupations have different n_states')
    call assert( size( complete_gnd_set_lapwlo, 3 ) == size( occupations, 2 ), &
      'complete_gnd_set_lapwlo and occupations have different n_kpts')

    this%eps_occ = occs_tol
    allocate( this%occupations, source = occupations )
    allocate( this%groundstate(n_gnd_states, n_gnd_states, first_kpt:last_kpt), source = zzero )
    do i = 1, n_gnd_states
      this%groundstate(i, i, :) = zone
    end do
    this%kset = kset

    ! TODO: issue a warning if n_active_states reduced nempty given in the input.xml
    n_active_states = last_occupied_for_all_ranks( occupations, occs_tol, mpiglobal )
    allocate( this%active(n_gnd_states, n_active_states-n_frozen, first_kpt:last_kpt), &
      source = this%groundstate(:, n_frozen + 1 : n_active_states, :) )
    if ( save_needed ) allocate( this%active_save, source = this%active )
    if ( n_frozen > 0 ) allocate( this%frozen(n_gnd_states, n_frozen, first_kpt:last_kpt), &
      source = this%groundstate(:, 1 : n_frozen, :) )
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
  logical function wavefunction_set_expanded_in_lapwlo( this ) result( expanded_in_lapwlo )
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
  pure logical function wavefunction_set_has_frozen( this ) result( has_frozen )
    class(wavefunction_set), intent(in) :: this

    has_frozen = allocated( this%frozen )
  end function

  !> Returns the number of frozen states
  pure integer(i32) function wavefunction_set_n_frozen( this ) result( n_frozen )
    class(wavefunction_set), intent(in) :: this

    n_frozen = 0
    if ( this%has_frozen() ) n_frozen = size( this%frozen, 2 )
  end function

  !> Returns the number of empty states
  pure integer(i32) function wavefunction_set_n_empty( this ) result( n_empty )
    class(wavefunction_set), intent(in) :: this
    
    n_empty = size( this%groundstate, 2 ) - this%n_occupied()
  end function

  !> Returns the number of occupied states
  pure integer(i32) function wavefunction_set_n_occupied( this ) result( n_occupied )
    class(wavefunction_set), intent(in) :: this
    
    n_occupied = this%n_frozen() + this%n_active()
  end function

  !> Returns the number of active states
  pure integer(i32) function wavefunction_set_n_active( this ) result( n_active)
    class(wavefunction_set), intent(in) :: this

    n_active = size( this%active, 2 )
  end function

  !> Checks the k dimensions and returns the number of k points
  integer(i32) function wavefunction_set_n_kpts( this ) result( n_kpts )
    class(wavefunction_set), intent(in) :: this

    call assert( size( this%active, 3 ) == size( this%groundstate, 3 ), &
      "active and groundstate must have the same number of elements along 3rd dim." )

    if( allocated( this%active_save ) ) call assert( size( this%active, 3 ) == size( this%active_save, 3 ), &
      "active and active_save must have the same number of elements along 3rd dim." )

    if ( this%has_frozen() ) call assert( size( this%active, 3 ) == size( this%frozen, 3 ), &
      "active and frozen must have the same number of elements along 3rd dim." )

    n_kpts = size( this%active, 3 )
  end function

  !> Return the index of the first \( \mathbf{k} \)-point assigned to this MPI rank.
  pure integer(i32) function wavefunction_set_first_kpt( this ) result( first_kpt )
    class(wavefunction_set), intent(in) :: this

    first_kpt = lbound( this%active, 3 )
  end function

  !> Return the index of the last \( \mathbf{k} \)-point assigned to this MPI rank.
  pure integer(i32) function wavefunction_set_last_kpt( this ) result( last_kpt )
    class(wavefunction_set), intent(in) :: this

    last_kpt = ubound( this%active, 3 )
  end function

  !> Returns basis size
  pure integer(i32) function wavefunction_set_n_basis( this ) result( n_basis)
    class(wavefunction_set), intent(in) :: this

    n_basis = size( this%active, 1 )
  end function

  !> Returns the position of the first active state
  pure integer(i32) function wavefunction_set_first_active( this ) result( first_active )
    class(wavefunction_set), intent(in) :: this
    
    first_active = this%n_frozen() + 1
  end function

  !> Wrapper for calling [[rttddft_io_unformatted(module):write_wavefunction]]
  subroutine wavefunction_set_write( this, mpi_env, handler )
    class(wavefunction_set), intent(inout) :: this
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    associate( ki => this%first_kpt(), kf => this%last_kpt() )
    call write_wavefunction( t, ki, this%kset%vkl(:, ki:kf), &
      this%active, mpi_env, handler, this%kset%nkpt )
    if( allocated( this%active_save ) ) call write_wavefunction( t_minus_dt, &
      ki, this%kset%vkl(:, ki:kf), this%active_save, mpi_env, handler, this%kset%nkpt )
    end associate
  end subroutine

  !> Wrapper for calling [[rttddft_io_unformatted(module):read_wavefunction]]
  subroutine wavefunction_set_read( this, mpi_env, handler )
    class(wavefunction_set), intent(inout) :: this
    !> MPI environment (needed to write in parallel over MPI procs.)
    type(mpiinfo), intent(in) :: mpi_env
    !> File handler
    type(file_handler), optional, intent(in) :: handler

    associate( ki => this%first_kpt(), kf => this%last_kpt() )
    call read_wavefunction( t, ki, this%kset%vkl(:, ki:kf), this%active, &
      mpi_env, handler )
    if( allocated( this%active_save ) ) call read_wavefunction( t_minus_dt, &
      ki, this%kset%vkl(:, ki:kf), this%active_save, mpi_env, handler )
    end associate
  end subroutine

  !> Compute the number of excitations in RT-TDDFT by projecting the
  !> time-evolved wavefunctions onto the ground state at \(t=0\). 
  !> First, the number of electrons electrons in the ground state is obtained as
  !> \[
  !>  N_{gs}(t)= \sum_{j\mathbf{k}}^{j\, occ} w_\mathbf{k} f_{j\mathbf{k}}
  !>          | \langle \psi_{j\mathbf{k}}(0) | \psi_{i\mathbf{k}}(t)\rangle |^2
  !> \]
  !> Then, the total number of electrons is calculated as
  !> \[
  !>  N_{tot}(t) = \sum_{j\mathbf{k}}^{j\, occ} w_\mathbf{k} f_{j\mathbf{k}}
  !>          | \langle \psi_{j\mathbf{k}}(t) | \psi_{i\mathbf{k}}(t)\rangle |^2
  !> \]
  !> Finally, the number of excitations is then given by
  !> \[
  !>  N_{exc}(t) = N_{tot}(t) - N_{gs}(t)
  !> \]
  subroutine wavefunction_set_obtain_number_excitations( this, overlap, mpi_env, &
    & n_exc, n_gs )
    !> Basis-expansion coefficients of the KS-wavefunctions.
    class(wavefunction_set), intent(in) :: this
    !> Overlap matrices
    class(overlap_set), intent(in) :: overlap
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> number of excited electrons
    real(dp), intent(out) :: n_exc
    !> number of electrons on the groundstate state
    real(dp), intent(out) :: n_gs

    integer(i32) :: ik, n_kpt, shift
    real(dp) :: buffer(2)
    real(dp), allocatable :: occ(:, :), aux_tot(:), aux_exc(:)
    complex(dp), allocatable :: proj(:, :, :)
    
    n_kpt = this%n_kpts()
    shift = this%first_kpt() - 1

    allocate( aux_tot(n_kpt), source = 0._dp )
    allocate( aux_exc(n_kpt), source = 0._dp )
    call obtain_projection_coefficients( this%groundstate, overlap%array, this%active, proj )
    call obtain_occupations( proj, this%occupations(this%first_active(): this%n_occupied(), :), occ )
    if ( this%has_frozen() ) &
      occ(1 : this%n_frozen(), :) = occ(1 : this%n_frozen(), :) + this%occupations(1 : this%n_frozen(), :)
    do concurrent (ik = 1:n_kpt)
      aux_tot(ik) = sum( occ(:, ik) )
      aux_exc(ik) = sum( occ(:, ik), this%occupations(:, ik+shift) <= this%eps_occ )
    end do
    associate( wkpt => this%kset%wkpt, ki => this%first_kpt(), kf => this%last_kpt() )
      n_exc = dot_product( wkpt(ki:kf), aux_exc )
      n_gs = dot_product( wkpt(ki:kf), aux_tot ) - n_exc
    end associate
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
  pure integer(i32) function last_occupied_for_current_rank( occupations, occs_tol )
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

  !> (private) Returns the index of the last occupied state for all MPI ranks
  integer(i32) function last_occupied_for_all_ranks( occupations, occs_tol, mpi_env ) result( last_occupied )
    !> State occupations array (n_states, n_kpt)
    real(dp), contiguous, intent(in) :: occupations(:, :)
    !> Minimal value of occupation for the state to be 'occupied'
    real(dp), intent(in) :: occs_tol
    !> MPI environment
    type(mpiinfo), intent(inout) :: mpi_env

    integer(i32), allocatable :: buffer(:)

    last_occupied = last_occupied_for_current_rank( occupations, occs_tol )
    ! Force the same number of active states over all MPI ranks
    call xmpi_allgather( mpi_env, last_occupied, buffer )
    last_occupied = maxval( buffer )
  end function

end module rttddft_Wavefunction

