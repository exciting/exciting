module MD
  use asserts, only: assert
  use constants, only: y00
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  use modinput, only: input
  use precision, only: dp, i32
  implicit none

  private 
  ! Types
  public :: force, MD_input_keys, MD_timing
  ! Procedures
  public :: obtain_core_corrections, &
            obtain_force_ext, &
            obtain_Hellmann_Feynman_force, &
            obtain_valence_corrections_part1, &
            obtain_valence_corrections_part2

  type force
    !> External force
    real(dp), allocatable     :: EXT(:,:)
    !> Hellmann-Feynman force
    real(dp), allocatable     :: HF(:,:)
    !> Core corrections to the force
    real(dp), allocatable     :: core(:,:)
    !> Valence corrections
    real(dp), allocatable     :: val(:,:)
    !> Total force = sum of individual components
    real(dp), allocatable     :: total(:,:)
    !> Total force in the previous time step
    real(dp), allocatable     :: total_save(:,:)
  contains
    procedure, public  :: allocate_arrays => force_allocate_arrays
    procedure, private :: deallocate_arrays => force_deallocate_arrays
    procedure, public  :: evaluate_total_force => force_evaluate_total_force
    procedure, public  :: save_total_force => force_save_total_force
    final              :: force_destructor
  end type

  type trajectory
    integer(i32), private     :: n_atoms
    !> Positions of the nuclei in cartesian coordinates
    real(dp), allocatable     :: positions(:,:)
    !> Velocities in cartesian coordinates
    real(dp), allocatable     :: velocities(:,:)
  contains
    procedure, private :: deallocate_arrays => trajectory_deallocate_arrays
    final              :: trajectory_destructor
  end type

  type MD_input_keys
    logical           :: on
    logical           :: print_all_force_components
    logical           :: update_overlap
    logical           :: update_pmat
    logical           :: basis_derivative
    logical           :: valence_corrections
    logical           :: core_corrections
    character(len=80) :: MD_type
    character(len=80) :: integration_algorithm
    real(dp)          :: time_step
  contains
    procedure         :: parse_input => MD_input_keys_parse_input
  end type

  !> Type for store timings spent in different subroutines of a MD
  type MD_timing
    !> Overall time spent in an MD step
    real(dp) :: t_MD_step
    !> Time spent to evaluate the 1st part of the total force
    real(dp) :: t_MD_1st
    !> Time spent to evaluate the 2nd part of valence corrections to the total force
    real(dp) :: t_MD_2nd
    !> Time spent to sum all contributions to the total force
    real(dp) :: t_MD_sumforces
    !> Time spent to update the positions of the ions
    real(dp) :: t_MD_moveions
    !> Time spent to update the basis after moving the ions
    real(dp) :: t_MD_updateBasis
  end type 

  contains 
  subroutine force_allocate_arrays( this, n_atoms )
    class(force), intent(inout) :: this
    integer, intent(in) :: n_atoms

    allocate( this%EXT(3, n_atoms), source = 0._dp )
    allocate( this%HF(3, n_atoms), source = 0._dp )
    allocate( this%core(3, n_atoms), source = 0._dp )
    allocate( this%val(3, n_atoms), source = 0._dp )
    allocate( this%total(3, n_atoms), source = 0._dp )
    allocate( this%total_save(3, n_atoms), source = 0._dp )
  end subroutine

  subroutine force_deallocate_arrays( this )
    class(force), intent(inout) :: this

    if( allocated( this%EXT) ) &
      deallocate( this%EXT, this%HF, this%core, this%val, this%total, this%total_save )
  end subroutine 

  impure elemental subroutine force_destructor( this )
    type(force), intent(inout) :: this

    call this%deallocate_arrays
  end subroutine 

  subroutine force_evaluate_total_force( this )
    class(force), intent(inout) :: this

    this%total = this%HF + this%EXT + this%core + this%val
  end subroutine 

  subroutine force_save_total_force( this )
    class(force), intent(inout) :: this

    this%total_save = this%total
  end subroutine 

  subroutine trajectory_deallocate_arrays( this )
    class(trajectory), intent(inout) :: this

    if( allocated(this%positions) ) deallocate( this%positions, this%velocities )
  end subroutine

  impure elemental subroutine trajectory_destructor( this )
    type(trajectory), intent(inout) :: this

    call this%deallocate_arrays
  end subroutine

  subroutine MD_input_keys_parse_input( this )
    class(MD_input_keys), intent(inout) :: this

    this%on = associated( input%MD )
    if( this%on ) then
      this%print_all_force_components = input%MD%printAllForces
      this%update_overlap = input%MD%updateOverlap
      this%update_pmat = input%MD%updatePmat
      this%basis_derivative = input%MD%basisDerivative
      this%valence_corrections = input%MD%coreCorrections
      this%core_corrections = input%MD%valenceCorrections
      this%MD_type = input%MD%type
      this%integration_algorithm = input%MD%integrationAlgorithm
      this%time_step = input%MD%timeStep
    end if
  end subroutine

  !> Obtain force_ext due to an electric field as 
  !> \( \mathbf{F}_{ext} = q\mathbf{E} \)
  subroutine obtain_force_ext( charge, e_field, force_ext )
    !> Particle charge
    real(dp), intent(in)  :: charge
    !> Electric field with its x, y, and z components
    real(dp), intent(in)  :: e_field(3)
    !> Force due to the electric field
    real(dp), intent(out) :: force_ext(3)
    
    force_ext = charge * e_field
  end subroutine

  !> Obtain the Hellmann-Feynman force acting on a specific atom
  !> \[ \mathbf{F}_{HF} = Z\lim_{\mathbf{r}\to 0} \nabla V_C(\mathbf{r}) \]
  !> If we expand \(V_C\) in spherical harmonics, we just need the components
  !> with \(l=1\), due to the limit \(\mathbf{r}\to 0\). This limit also implies 
  !> that only \(l=0\) of \( \nabla V_C(\mathbf{r}) \) is relevant.
  !> The trick here is to use also the component \(l=0\), which is a constant, 
  !> for the gradient. When the gradient is carried out, we take the \(l=0\) 
  !> component, and evaluate the limit with the point with smallest \(r\).
  subroutine obtain_Hellmann_Feynman_force( Z, radial_grid, vc_lm, force_HF )
    !> Atomic number
    real(dp), intent(in)  :: Z
    !> Radial grid of the corresponding muffin-tin sphere
    real(dp), intent(in)  :: radial_grid(:)
    !> lm-components of the \(V_C\) potential, where
    !> \(V_C = V_{Hartree} + V_{nuclear}\)
    real(dp), intent(in)  :: vc_lm(:, :)
    !> Hellmann-Feynman force
    real(dp), intent(out) :: force_HF(3)
    ! Local variables
    integer(i32)           :: nr
    integer(i32),parameter :: l_max = 1
    integer(i32),parameter :: lm_max = (l_max + 1)**2
    real(dp), allocatable  :: grad(:, :, :)
    real (dp), external    :: rfmtinp
    
    nr = size( radial_grid )
    call assert( size(vc_lm, 1) >= lm_max, '1st dim of vc_lm must be >= lm_max')
    call assert( size(vc_lm, 2) == nr, '2nd dim of vc_lm must contain nr elements' )
    
    allocate( grad(lm_max, nr, 3), source = 0._dp )

    call gradrfmt( 1, nr, radial_grid, lm_max, nr, vc_lm(1:lm_max, :), grad )
    force_HF(1:3) = Z * grad(1,1,1:3) * y00
  end subroutine

  !> Obtain the core corrections to the force acting on a given atom
  !> \[ \mathbf{F}_{core} = -\int_{MT_J} \mathrm{d}\mathbf{r} n_c \nabla v_{KS}
  !> = \int_{MT_J} \mathrm{d}\mathbf{r} v_{KS} \nabla n_c \]
  subroutine obtain_core_corrections( radial_grid, rho_core, vKS_MT, force_core )
    !> Radial grid of the corresponding muffin-tin
    real(dp), intent(in)  :: radial_grid(:)
    !> Core density
    real(dp), intent(in)  :: rho_core(:)
    !> Muffin-tin part of the Kohn-Sham potential
    real(dp), intent(in)  :: vKS_MT(:, :)
    !> Core corrections
    real(dp), intent(out) :: force_core(3)
    ! Local variables
    integer(i32)           :: nr, j
    integer(i32),parameter :: l_max = 1
    integer(i32),parameter :: lm_max = (l_max + 1)**2
    real(dp), allocatable  :: rho_mt(:, :)
    real(dp), allocatable  :: grad(:, :, :)
    real (dp), external    :: rfmtinp
    
    nr = size( radial_grid )
    call assert( size(rho_core) == nr, 'Size of rho_core must be nr' )
    call assert( size(vKS_MT, 1) >= lm_max, '1st dim of vKS_MT must be >= lm_max')
    call assert( size(vKS_MT, 2) == nr, '2nd dim of vKS_MT must contain nr elements' )
    
    allocate( rho_mt(lm_max, nr), source = 0._dp )
    allocate( grad(lm_max, nr, 3), source = 0._dp )

    rho_mt(1,:) = rho_core(:)/y00
    call gradrfmt(1, nr, radial_grid, lm_max, nr, rho_mt, grad )
    do j = 1, 3
      force_core(j) = rfmtinp( 1, 1, nr, radial_grid, &
        & lm_max, vKS_MT(1:lm_max, :), grad(:, :, j) )
    end do

  end subroutine

  !> Obtain the 1st part of the valence corrections to the force acting on a specific atom
  !> \[ \mathbf{F}_{val,1}= \int_{MT_J}\mathrm{d}\mathbf{r}(\nabla n_v(\mathbf{r}))
  !> \left( v_{KS}(\mathbf{r})+term \right)\]
  subroutine obtain_valence_corrections_part1( radial_grid, rho_core, rho_MT, &
      vKS_MT, term, force_val1 )
    !> Radial grid of the corresponding muffin-tin sphere
    real(dp), intent(in)  :: radial_grid(:)
    !> Core density
    real(dp), intent(in)  :: rho_core(:)
    !> Muffin-tin part of the total density (core + valence)
    real(dp), intent(in)  :: rho_MT(:, :)
    !> Muffin-tin part of the Kohn-Sham potential
    real(dp), intent(in)  :: vKS_MT(:, :)
    !> Term to be added to the l=0 component of `vKS_MT`
    real(dp), intent(in)  :: term
    !> Valence corrections (1st part)
    real(dp), intent(out) :: force_val1(3)
    ! Local variables
    integer(i32)           :: nr, j
    integer(i32)           :: l_max
    integer(i32)           :: lm_max
    real(dp), allocatable  :: rho_val(:, :), vaux(:, :)
    real(dp), allocatable  :: grad(:, :, :)
    real(dp), external     :: rfmtinp

    nr = size( radial_grid )
    lm_max = size( rho_MT, 1 )
    l_max = int( sqrt( dble(lm_max) ) ) - 1
    call assert( lm_max == (l_max+1)**2, 'lm_max must be a perfect square' )
    call assert( size(rho_core) == nr, 'Size of rho_core must be nr' )
    call assert( size(rho_MT, 1) == lm_max, '1st dim of rho_MT must be == lm_max')
    call assert( size(rho_MT, 2) == nr, '2nd dim of rho_MT must contain nr elements' )
    call assert( size(vKS_MT, 1) == lm_max, '1st dim of vKS_MT must be >= lm_max')
    call assert( size(vKS_MT, 2) == nr, '2nd dim of vKS_MT must contain nr elements' )

    allocate( rho_val, source = rho_MT )
    allocate( vaux, source = vKS_MT)
    allocate( grad(lm_max, nr, 3), source = 0._dp )

    ! Valence charge = total charge - core charge (rho_core has only l=0 component)
    rho_val(1,:) = rho_val(1,:) - rho_core/y00
    call gradrfmt( l_max, nr, radial_grid, lm_max, nr, rho_val, grad )
    vaux(1,:) = vaux(1,:) + term/y00
    do j = 1, 3
      force_val1(j) = rfmtinp( 1, l_max, nr, radial_grid, lm_max, vaux, grad(:,:,j) )
    end do
  end subroutine

  !> Obtain the 2nd part of the valence corrections to the force acting on a specific atom
  !> \[ F = \sum_j f_{j}(\psi_j^\dagger)(\mathcal{H}-\mathcal{S})(\psi_j)
  subroutine obtain_valence_corrections_part2( H, S, psi, occ, F )
    !> \(\mathcal{H}\) matrix
    complex(dp), intent(in)  :: H(:, :, :)
    !> \(\mathcal{S}\) matrix
    complex(dp), intent(in)  :: S(:, :, :)
    !> \(\psi_j\)Wavefunction coefficients
    complex(dp), intent(in)  :: psi(:, :)
    !> \(f_{j}\): Occupation factors
    real(dp), intent(in)     :: occ(:)
    !> \(F\): Contribution to the valence corrections (as given by the formula)
    real(dp), intent(out)    :: F(3)
    ! Local variables
    integer(i32)             :: i, j, m, n
    real(dp), allocatable    :: aux(:)
    complex(dp), allocatable :: diff(:, :), prod(:, :)

    m = size( psi, 1 )
    n = size( psi, 2 )
    allocate( diff(m, m), prod(m, n), aux(n) )
    call assert( size(H,1) == m, 'H must have m elements along 1st dim')
    call assert( size(H,2) == m, 'H must have m elements along 2nd dim')
    call assert( size(H,3) == 3, 'H must have 3 elements along 3rd dim')
    call assert( size(S,1) == m, 'S must have m elements along 1st dim')
    call assert( size(S,2) == m, 'S must have m elements along 2nd dim')
    call assert( size(S,3) == 3, 'S must have 3 elements along 3rd dim')
    call assert( size(occ) == n, 'occ must have n elements')

    ! loop over x, y, z
    do j = 1, 3
      diff = H(:,:,j) - S(:,:,j)
      call hermitian_matrix_multiply( diff, psi, prod )
      forall (i = 1:n) aux(i) = dble( dot_product( psi(:,i), prod(:,i) ) )
      F(j) = dot_product( occ, aux )
    end do
  end subroutine
end module