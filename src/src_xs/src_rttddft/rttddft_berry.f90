!> This module handles all Berry-phase-related calculations required for RT-TDDFT 
!> using the dynamical Berry phase approach to describe the interaction with the external field.
module rttddft_berry
  use asserts, only: assert
  use constants, only: zone, zzero, zi, fourpi
  use determinant, only: determinant_LU
  use exciting_mpi, only: xmpi_allgatherv
  use inverse, only: inverse_LU
  use mod_kpointset, only: k_set
  use mod_lattice, only: avec
  use modmpi, only: mpi_env_k
  use precision, only: dp, i32
  use rttddft_electric_field, only: Electric_Field
  use rttddft_Wavefunction, only: wavefunction_set
  use xlapack, only: matrix_multiply

  implicit none

  private
  public :: get_td_overlap_det_and_berry_coupling_term

  !> Number of cartesian directions
  integer(i32), parameter :: n_cartesian = 3
  !> Treshold for the field to be considered non-zero
  real(dp), parameter :: eps_field = 1.e-14_dp

contains

  !> Wrapper for calling the private `rttddft_berry` routines.
  !> First, active states are gathered from all MPI processes.
  !> Then, [[get_jumps_array]] is called to determine how many Cartesian directions
  !> and jumps are active, based on the external field and the \( \mathbf{k} \)-grid dimensions.
  !> Next, the time-dependent overlap is computed in [[build_td_overlap_and_det]].
  !> Finally, the interaction term is obtained via [[get_berry_coupling_term]].
  !> For further details, please refer to the documentation of the respective routines.
  subroutine get_td_overlap_det_and_berry_coupling_term( first_kpt, e_vec, pws_for_berry_phase, psi, &
      kset, k_ptrs, td_overlap_det, berry_coupling_term, use_save )
    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Electric field \( \mathbf{E} (t) \)
    type(Electric_Field), intent(in) :: e_vec
    !> Time-independent matrix \( W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pr} \)
    complex(dp), contiguous, intent(in) :: pws_for_berry_phase(:, :, :, :, :)
    !> Basis-expansion coefficients of the KS-WFs
    class(wavefunction_set), intent(in) :: psi
    !> Set of \( \mathbf{k} \) points used throughout the module
    type(k_set), intent(in) :: kset
    !> Pointers to the neighbouring \( \mathbf{k} \) points (nkpt, 3, 4)
    integer(i32), contiguous, intent(in) :: k_ptrs(:, :, :)
    !> Determinants of the time-dependent overlap matrices gathered from all MPI procs, 
    !> det \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{+}} \) (3, nkpt)
    complex(dp), contiguous, intent(out) :: td_overlap_det(:, :)
    !> Resulting interaction term
    complex(dp), contiguous, intent(out) :: berry_coupling_term(:, :, :)
    !> If `.true.`, `active_save` wavefunction component should be used instead of `active` one
    logical, optional, intent(in) :: use_save

    complex(dp), allocatable :: all_active_states(:, :, :), td_overlap(:, :, :, :, :)
    integer(i32) :: max_jumps(n_cartesian), last_kpt
    logical :: use_save_local

    last_kpt = first_kpt + psi%n_kpts() - 1
    use_save_local = .false.
    if ( present( use_save ) ) use_save_local = use_save

    call assert( psi%n_basis() == size( berry_coupling_term, 1 ), "n_basis is different for psi and berry_coupling_term" )
    call assert( psi%n_basis() == size( pws_for_berry_phase, 1 ), "n_basis is different for psi and pws_for_berry_phase" )
    call assert( size( k_ptrs, 1) == size( td_overlap_det, 2), "n_kpt is different for k_ptrs and td_overlap_det" )
    call assert( psi%n_kpts() == size( berry_coupling_term, 3 ), "n_kpt is different for psi and berry_coupling_term" )

    allocate( td_overlap(psi%n_occupied(), psi%n_occupied(), first_kpt : last_kpt, n_cartesian, 4) )
    allocate( all_active_states(psi%n_basis(), psi%n_active(), kset%nkpt), source = zzero )

    ! gather active states across MPI processes
    if ( use_save_local ) then
      all_active_states(:, :, first_kpt : last_kpt) = psi%active_save
    else
      all_active_states(:, :, first_kpt : last_kpt) = psi%active
    end if

    call xmpi_allgatherv( mpi_env_k, all_active_states, size( psi%active ) )

    call get_jumps_array( kset%ngridk, e_vec, max_jumps )
    call build_td_overlap_and_det( first_kpt, all_active_states, &
      pws_for_berry_phase, k_ptrs, max_jumps, td_overlap, td_overlap_det )
    call get_berry_coupling_term( first_kpt, all_active_states, pws_for_berry_phase, e_vec, &
      k_ptrs, kset%ngridk, max_jumps, td_overlap, berry_coupling_term )

  end subroutine


  !> Get the maximum number of jumps along the \( \mathbf{k} \) grid 
  !> \( \sigma = +, -, +2, -2 \) needed for each cartesian direction \( \alpha \) 
  !> before evaluating the time-dependent overlap 
  !> \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} \) (see [[build_td_overlap_and_det]]).
  !> 1 jump along a direction is needed to evaluate polarization; 2 -- when the field is on, but 
  !> only 3 \( \mathbf{k} \) points are available for numerical derivative evaluation, 
  !> 4 -- when the field is on, and atleast 5 \( \mathbf{k} \) points are available.
  subroutine get_jumps_array( k_grid_size, e_vec, max_jumps )
    !> Dimensions of the \( \mathbf{k} \) grid (3)
    integer(i32), intent(in) :: k_grid_size(:)
    !> Electric field \( \mathbf{E} (t) \)
    type(Electric_Field), intent(in) :: e_vec
    !> Max jumps array (3)
    integer(i32), intent(out) :: max_jumps(:)

    integer(i32) :: k_dir

    max_jumps = 0
    do k_dir = 1, n_cartesian
      if ( k_grid_size(k_dir) < 2 ) cycle
      
      max_jumps(k_dir) = 1
      if ( abs( dot_product( e_vec%components, avec(:, k_dir) ) ) > eps_field ) then
        ! field is on, numerical k derivative will be needed
        max_jumps(k_dir) = 2
        if ( k_grid_size(k_dir) > 4 ) max_jumps(k_dir) = 4 ! 5 k points for the 2nd order derivative
      end if
    end do

  end subroutine

  !> Evaluate the time-dependent \( n_{\rm occupied} \times n_{\rm occupied} \) overlap 
  !> matrices \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} \) ( \( \alpha = 1, 2, 3 \), \( \sigma = +, -, +2, -2 \) ) 
  !> for each \( \mathbf{k} \) point of the current MPI process. The set of used directions \( \alpha \) and jumps 
  !> \( \sigma \) is precalculated in [[get_jumps_array]]. There are \( n_{\rm frozen} \) 
  !> frozen states (out of \( n_{\rm occupied} \)), therefore overlap consists of four blocks: frozen, top, bottom, active:
  !> \[
  !>     S^{\mathbf{k}, \mathbf{k}_{\alpha}^{\sigma}, {\rm frozen}}_{ln} (t) = 
  !> W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{ln}, \; \; l,n = 1, \dots, n_{\rm frozen},
  !> \]
  !> \[
  !>     S^{\mathbf{k}, \mathbf{k}_{\alpha}^{\sigma}, {\rm top}}_{ln} (t) = \sum_r W^{\mathbf{k} 
  !> \mathbf{k}_{\alpha}^{\sigma}}_{lr} c_{\mathbf{k}_{\alpha}^{\sigma} rn}(t), \; l = 1, \dots, n_{\rm frozen}, 
  !> \; n = n_{\rm frozen} + 1, \dots, n_{\rm occupied},
  !> \]
  !> \[
  !>     S^{\mathbf{k}, \mathbf{k}_{\alpha}^{\sigma}, {\rm bottom}}_{ln} (t) = 
  !> \sum_p c^*_{\mathbf{k} pl} (t) W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pn}, 
  !> \; l = n_{\rm frozen} + 1, \dots, n_{\rm occupied}, \; n = 1, \dots, n_{\rm frozen},
  !> \]
  !> \[
  !>     S^{\mathbf{k}, \mathbf{k}_{\alpha}^{\sigma}, {\rm active}}_{ln} (t) = 
  !> \sum_p c^*_{\mathbf{k} pl} (t) \left( \sum_r W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pr} 
  !> c_{\mathbf{k}_{\alpha}^{\sigma} rn}(t) \right ) \; l = n_{\rm frozen} + 1, 
  !> \dots, n_{\rm occupied}, \; n = n_{\rm frozen} + 1, \dots, n_{\rm occupied},
  !> \]
  !> where 
  !> \[
  !> W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pr} = \left \langle 
  !> u_{\mathbf{k} p} | u_{\mathbf{k}_{\alpha}^{\sigma} r} \right \rangle
  !> \]
  !> is the matrix pre-evaluated in the [[calc_pw_mes]] from periodic parts 
  !> \( u_{\mathbf{k} p} \) of the unperturbed KS states. After construction of the 
  !> overlap, evaluate det \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{+}} \).
  subroutine build_td_overlap_and_det( first_kpt, active_states, pws_for_berry_phase, &
      k_ptrs, max_jumps, td_overlap, td_overlap_det )
    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Active states expansion coefficients \( c_{\mathbf{k}} \) from all MPI procs
    complex(dp), contiguous, intent(in) :: active_states(:, :, :)
    !> Time-independent matrix \( W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pr} \)
    complex(dp), contiguous, intent(in) :: pws_for_berry_phase(:, :, first_kpt:, :, :)
    !> Pointers to the neighbouring \( \mathbf{k} \) points (nkpt, 3, 4)
    integer(i32), contiguous, intent(in) :: k_ptrs(:, :, :)
    !> Max jumps array (3)
    integer(i32), contiguous, intent(in) :: max_jumps(:)
    !> Time-dependent overlap matrix \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} \)
    !> (n_occupied, n_occupied, first_kpt:last_kpt, 3, 4)
    complex(dp), contiguous, intent(out) :: td_overlap(:, :, first_kpt: ,:, :)
    !> Determinants of the time-dependent overlap matrices gathered from all MPI procs, 
    !> det \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{+}} \) (3, nkpt)
    complex(dp), contiguous, intent(out) :: td_overlap_det(:, :)

    integer(i32) :: k_dir, k_jump, k_left, k_right, last_kpt, n_frozen
    complex(dp), allocatable :: g_matrix(:, :, :)

    last_kpt = ubound( pws_for_berry_phase, 3 )
    ! n_frozen = n_occupied - n_active:
    n_frozen = size( td_overlap, 1 ) - size( active_states, 2 )
    
    allocate( g_matrix(size( pws_for_berry_phase, 1 ), size( active_states, 2 ), 4) )
    td_overlap = zzero
    td_overlap_det(:, first_kpt : last_kpt) = zone
    
    !$omp parallel default(none) private(k_left, k_dir, g_matrix, k_jump, k_right), &
    !$omp shared(first_kpt, last_kpt, max_jumps, n_frozen, td_overlap_det, &
    !$omp active_states, pws_for_berry_phase, td_overlap, k_ptrs)
    !$omp do
    do k_left = first_kpt, last_kpt
      do k_dir = 1, n_cartesian

        do k_jump = 1, max_jumps(k_dir) ! +1 -1 +2 -2
          k_right = k_ptrs(k_left, k_dir, k_jump)

          if ( n_frozen > 0 ) then
            ! time-independent frozen part (n_frozen x n_frozen):
            td_overlap(1 : n_frozen, 1 : n_frozen, k_left, k_dir, k_jump) = &
              pws_for_berry_phase(1 : n_frozen, 1 : n_frozen, k_left, k_dir, k_jump)
            ! top part (n_frozen x n_active):
            call matrix_multiply( pws_for_berry_phase(1 : n_frozen , :, &
              k_left, k_dir, k_jump), active_states(:, :, k_right), &
              td_overlap(1 : n_frozen, n_frozen + 1 :, k_left, k_dir, k_jump) )
            ! bot part (n_active x n_frozen):
            call matrix_multiply( active_states(:, :, k_left), &
              pws_for_berry_phase(:, 1 : n_frozen, k_left, k_dir, k_jump), &
              td_overlap(n_frozen + 1 :, 1 : n_frozen, k_left, k_dir, k_jump), trans_A = 'C' )
          end if

          ! g_matrix(n_basis x n_active):
          call matrix_multiply( pws_for_berry_phase(:, :, k_left, k_dir, k_jump), &
            active_states(:, :, k_right), g_matrix(:, :, k_jump) )
          ! active part (n_active x n_active):
          call matrix_multiply( active_states(:, :, k_left), g_matrix(:, :, k_jump), &
            td_overlap(n_frozen + 1 :, n_frozen + 1 :, k_left, k_dir, k_jump), trans_A = 'C' )
        end do ! cycle over 2 or 4 derivative steps, k_jump

        if ( max_jumps(k_dir) > 0 ) td_overlap_det(k_dir, k_left) = &
          determinant_LU( td_overlap(:, :, k_left, k_dir, 1) )

      end do ! cycle over 3 lattice k directions, k_dir
    end do ! cycle over k points, k_left
    !$omp end parallel

    call xmpi_allgatherv( mpi_env_k, td_overlap_det, size( td_overlap_det, 1 ) * (last_kpt - first_kpt + 1) )
  end subroutine

  !> Evaluate the time-dependent field coupling matrix
  !> \[
  !> \frac{i f}{4 \pi} \sum_{\alpha = 1}^3 (\mathbf{a}_{\alpha} \cdot \mathbf{E} (t) ) N_{\mathbf{k}_{\alpha}} 
  !> \sum_{\sigma = \pm} \sigma \Bigg\{ \sum_{n = n_{\rm frozen} + 1}^{n_{\rm occupied}} c^*_{\mathbf{k} mn} (t) 
  !> \left[ \sum_{l = 1}^{n_{\rm frozen}} W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pl} 
  !> (S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}})^{-1}_{ln} (t) + \sum_{l = n_{\rm frozen} + 1}^{n_{\rm occupied}} 
  !> G^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pl}(t)  (S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}})^{-1}_{ln} (t) \right]
  !> + \sum_{n = 1}^{n_{\rm frozen}} \delta_{mn} \left[ \sum_{l = 1}^{n_{\rm frozen}} 
  !> W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pl} (S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}})^{-1}_{ln} (t) + 
  !> \sum_{l = n_{\rm frozen} + 1}^{n_{\rm occupied}} G^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pl}(t) 
  !> (S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}})^{-1}_{ln} (t) \right] \Bigg\} + {\rm H.a.}
  !> \]
  !> where \( \mathbf{a}_{\alpha} \) are the lattice vectors, \( \mathbf{E} (t) \) 
  !> is the field strength, \( N_{\mathbf{k}_{\alpha}} \) is the number of \( \mathbf{k} \)
  !> points along direction \( \alpha \), \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} \) is 
  !> the time-dependent overlap built in [[build_td_overlap_and_det]], 
  !> \[
  !> G^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pl}(t) = \sum_{r} W^{\mathbf{k} 
  !> \mathbf{k}_{\alpha}^{\sigma}}_{pr} c_{\mathbf{k}_{\alpha}^{\sigma} rl}(t),
  !> \]
  !> and
  !> \[
  !> W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}}_{pr} = \left \langle 
  !> u_{\mathbf{k} p} | u_{\mathbf{k}_{\alpha}^{\sigma} r} \right \rangle
  !> \]
  !> is the matrix pre-evaluated in the [[calc_pw_mes]] from periodic parts 
  !> \( u_{\mathbf{k} p} \) of the unperturbed KS states.
  subroutine get_berry_coupling_term( first_kpt, active_states, pws_for_berry_phase, e_vec, &
      k_ptrs, k_grid_size, max_jumps, td_overlap, berry_coupling_term )

    !> The first \( \mathbf{k} \) point
    integer(i32), intent(in) :: first_kpt
    !> Active states expansion coefficients \( c_{\mathbf{k}} \) from all MPI procs
    complex(dp), contiguous, intent(in) :: active_states(:, :, :)
    !> Time-independent initial PW MEs \( W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} \)
    complex(dp), contiguous, intent(in) :: pws_for_berry_phase(:, :, first_kpt :, :, :)
    !> Electric field \( \mathbf{E} (t) \)
    type(Electric_Field), intent(in) :: e_vec
    !> Pointers to the neighbouring \( \mathbf{k} \) points (nkpt, 3, 4)
    integer(i32), contiguous, intent(in) :: k_ptrs(:, :, :)
    !> Dimensions of the \( \mathbf{k} \) grid (3)
    integer(i32), contiguous, intent(in) :: k_grid_size(:)
    !> Max jumps array (3)
    integer(i32), contiguous, intent(in) :: max_jumps(:)
    !> Time-dependent overlap matrix \( S^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} \) (n_occupied, n_occupied, :, 3, 4)
    complex(dp), contiguous, intent(in) :: td_overlap(:, :, first_kpt : ,:, :)
    !> Resulting interaction term
    complex(dp), contiguous, intent(out) :: berry_coupling_term(:, :, first_kpt :)

    real(dp), parameter :: four_over_three = 4._dp/3._dp, one_over_six = 1._dp/6._dp
    integer(i32) :: k_dir, k_jump, i, k_left, k_right, n_basis, &
      n_active, last_kpt, n_frozen, n_occupied
    complex(dp), allocatable :: td_overlap_inv(:, :), g_matrix(:, :), dir_term(:, :), &
      dir_sum(:, :), sigma_term(:, :, :), tmp_n_basis_n_fr(:, :), aux_n_basis_n_ac(:, :), &
      tmp_n_basis_n_ac(:, :)
    real(dp) :: factor(n_cartesian)

    berry_coupling_term = zzero
    if ( dot_product( e_vec%components, e_vec%components ) < eps_field ) return

    last_kpt = ubound( pws_for_berry_phase, 3 )
    n_basis = size( pws_for_berry_phase, 1 )
    n_active = size( active_states, 2 )
    n_occupied = size( td_overlap, 1 )
    n_frozen = n_occupied - n_active
    
    allocate( g_matrix(n_basis, n_active) )
    allocate( dir_sum(n_basis, n_basis), dir_term(n_basis, n_basis) )
    allocate( td_overlap_inv(n_occupied, n_occupied) )

    allocate( sigma_term(n_basis, n_basis, 4) )
    allocate( tmp_n_basis_n_fr(n_basis, n_frozen) )
    allocate( aux_n_basis_n_ac(n_basis, n_active) )
    allocate( tmp_n_basis_n_ac(n_basis, n_active) )

    do k_dir = 1, n_cartesian
      ! factor = N_k_alpha x (F(t) \cdot a_alpha)
      factor(k_dir) = real( k_grid_size(k_dir), dp ) * dot_product( e_vec%components, avec(:, k_dir) )
    end do

    !$omp parallel default(none) private(k_left, dir_term, dir_sum, k_dir, &
    !$omp g_matrix, k_jump, k_right, td_overlap_inv, sigma_term, tmp_n_basis_n_fr, &
    !$omp aux_n_basis_n_ac, tmp_n_basis_n_ac), &
    !$omp shared(first_kpt, last_kpt, max_jumps, k_grid_size, n_frozen, active_states, &
    !$omp pws_for_berry_phase, berry_coupling_term, k_ptrs, td_overlap, factor)
    !$omp do
    do k_left = first_kpt, last_kpt

      dir_sum = zzero
      do k_dir = 1, n_cartesian
        sigma_term = zzero
        do k_jump = 1, max_jumps(k_dir) ! +1 -1 +2 -2
          
          td_overlap_inv = inverse_LU( td_overlap(:, :, k_left, k_dir, k_jump) )
          k_right = k_ptrs(k_left, k_dir, k_jump)

          ! W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} x C = G^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} (n_basis x n_active):
          call matrix_multiply( pws_for_berry_phase(:, :, k_left, k_dir, k_jump), &
            active_states(:, :, k_right), g_matrix )
          
          ! G^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} x S^-1[active-active]
          call matrix_multiply( g_matrix, td_overlap_inv(n_frozen + 1 :, &
            n_frozen + 1 :), aux_n_basis_n_ac )

          if ( n_frozen > 0 ) then
            ! W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} x S^-1[frozen-active]
            call matrix_multiply( pws_for_berry_phase(:, 1 : n_frozen, k_left, k_dir, k_jump), &
              td_overlap_inv(1 : n_frozen, n_frozen + 1 :), tmp_n_basis_n_ac )
            aux_n_basis_n_ac = aux_n_basis_n_ac + tmp_n_basis_n_ac
          end if
          ! C^* x aux_n_basis_n_ac
          call matrix_multiply( aux_n_basis_n_ac, active_states(:, :, k_left), &
            sigma_term(:, :, k_jump), trans_B = 'C' )

          if ( n_frozen > 0 ) then
            ! W^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} x S^-1[frozen-frozen]
            call matrix_multiply( pws_for_berry_phase(:, 1 : n_frozen, k_left, k_dir, k_jump), &
              td_overlap_inv(1 : n_frozen, 1 : n_frozen), tmp_n_basis_n_fr )
            sigma_term(:, 1 : n_frozen, k_jump) = sigma_term(:, 1 : n_frozen, k_jump) + tmp_n_basis_n_fr

            ! G^{\mathbf{k} \mathbf{k}_{\alpha}^{\sigma}} x S^-1[active-frozen]
            call matrix_multiply( g_matrix, td_overlap_inv(n_frozen + 1 :, &
              1 : n_frozen), tmp_n_basis_n_fr )
            sigma_term(:, 1 : n_frozen, k_jump) = sigma_term(:, 1 : n_frozen, k_jump) + tmp_n_basis_n_fr
          end if
        end do ! cycle over 2 or 4 derivative steps, k_jump

        if ( abs( factor(k_dir) ) > eps_field ) then
          ! first order numerical derivative D(dk)
          ! see Eq. (96) of [PRB 63, 155107 (2001)]
          dir_term = sigma_term(:, :, 1) - sigma_term(:, :, 2)

          if ( max_jumps(k_dir) == 4 ) then
            ! second order numerical derivative, ( 4*D(dk) - D(2*dk) ) / 3
            ! see Eq. (97) of [PRB 63, 155107 (2001)]
            dir_term = four_over_three * dir_term - one_over_six * & ! 1/3 * 1/2 from N_kpt
              (sigma_term(:, :, 3) - sigma_term(:, :, 4))
          end if
          dir_sum = dir_sum + factor(k_dir) * dir_term
        end if
      end do ! cycle over 3 lattice k directions, k_dir
      
      berry_coupling_term(:, :, k_left) = zi / fourpi * dir_sum
      ! make the matrix hermitian
      berry_coupling_term(:, :, k_left) = berry_coupling_term(:, :, k_left) + &
        conjg( transpose( berry_coupling_term(:, :, k_left) ) )

    end do ! cycle over k points, k_left
    !$omp end parallel

  end subroutine

end module
