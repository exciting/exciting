module MD
  use asserts, only: assert
  use constants, only: y00
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  use mod_atoms, only: atposc, idxas, natoms, nspecies
  use modinput, only: input, structure_type
  use rttddft_electric_field, only: Electric_Field
  use precision, only: dp, i32
  use to_char_conversion, only: to_char
  
  implicit none

  private 

  integer(i32), parameter :: n_cartesian = 3
  
  ! Procedures
  public :: obtain_core_corrections, &
            force_ext, &
            obtain_Hellmann_Feynman_force, &
            obtain_valence_corrections_part1, &
            obtain_valence_corrections_part2

  type, public :: force
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

  type, public :: trajectory
    !> Positions of the nuclei in cartesian coordinates
    real(dp), allocatable     :: positions(:, :)
    !> Velocities in cartesian coordinates
    real(dp), allocatable     :: velocities(:, :)
  contains
    private
    procedure, public  :: allocate_arrays => trajectory_allocate_arrays
    procedure, public  :: assert_consistency => trajectory_assert_consistency
    procedure, private :: deallocate_arrays => trajectory_deallocate_arrays
    procedure, public  :: initialize => trajectory_initialize
    procedure, public  :: update_globals => trajectory_update_global_vars
    final              :: trajectory_destructor
  end type

  type, public :: MD_input_keys
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
  type, public :: MD_timing
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
  contains
    procedure :: reset_MD_timing
  end type 

contains
  pure subroutine reset_MD_timing( this )
    class(MD_timing), intent(inout) :: this

    this%t_MD_step = 0._dp
    this%t_MD_1st = 0._dp
    this%t_MD_2nd = 0._dp
    this%t_MD_sumforces = 0._dp
    this%t_MD_moveions = 0._dp
    this%t_MD_updateBasis = 0._dp

  end subroutine reset_MD_timing

  !> Allocate force arrays
  subroutine force_allocate_arrays( this, n_atoms, allocate_total )
    class(force), intent(inout) :: this
    !> Number of atoms
    integer(i32), intent(in) :: n_atoms
    !> If `.true.`, also allocate the `total` and `total_save` components.  
    !> `.false.` is the right option when restarting and old calculation
    logical, intent(in) :: allocate_total

    allocate( this%EXT(3, n_atoms), source = 0._dp )
    allocate( this%HF(3, n_atoms), source = 0._dp )
    allocate( this%core(3, n_atoms), source = 0._dp )
    allocate( this%val(3, n_atoms), source = 0._dp )
    if( allocate_total ) then
      allocate( this%total(3, n_atoms), source = 0._dp )
      allocate( this%total_save(3, n_atoms), source = 0._dp )
    end if
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

  !> Initialize from positions and velocities defined in the input file.
  !> This can be executed even before calling `init0`, which initializes global variables 
  !> such as `natoms(:)` and `idxas(:, :)`
  subroutine trajectory_initialize( this, structure )
    class(trajectory), intent(inout) :: this
    !> `structure` element in the input file
    type(structure_type), intent(in) :: structure

    integer(i32) :: is, ia, ias, n_species, n_atoms
    logical :: cartesian

    cartesian = structure%cartesian
    ias = 0
    n_species = size( structure%speciesarray )
    do is = 1, n_species
      n_atoms = size( structure%speciesarray(is)%species%atomarray )
      do ia = 1, n_atoms
        ias = ias + 1
        this%velocities(:, ias) = structure%speciesarray(is)%species%atomarray(ia)%atom%velocity
        if ( cartesian ) then
          this%positions(:, ias) = structure%speciesarray(is)%species%atomarray(ia)%atom%coord
        else
          this%positions(:, ias) = matmul( structure%crystal%basevect, structure%speciesarray(is)%species%atomarray(ia)%atom%coord )
        end if
      end do
    end do   
  end subroutine

  subroutine trajectory_update_global_vars( this )
    class(trajectory), intent(in) :: this

    integer(i32) :: is, ia, ias

    do is = 1, nspecies
      do ia = 1, natoms (is)
        ias = idxas (ia, is)
        atposc(:, ia, is) = this%positions(: , ias)
      end do
    end do
  end subroutine

  !> Trajectory asserts
  subroutine trajectory_assert_consistency( this )
    class(trajectory), intent(in) :: this

    call assert( size( this%positions, 1 ) == n_cartesian, "positions must have " // to_char(n_cartesian) // " components along 1st dim")
    call assert( size( this%velocities, 1 ) == n_cartesian, "velocities must have " // to_char(n_cartesian) // " components along 1st dim")
    call assert( size( this%velocities, 2 ) == size( this%positions, 2 ), "velocities and positions must have same size along 2nd dim")
  end subroutine

  subroutine trajectory_allocate_arrays( this, n_atoms )
    class(trajectory), intent(inout) :: this
    !> Number of atoms
    integer(i32), intent(in) :: n_atoms
    
    call this%deallocate_arrays( )
    allocate( this%positions(n_cartesian, n_atoms), this%velocities(n_cartesian, n_atoms) )
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
  pure function force_ext( charge, e_field )
    !> Particle charge
    real(dp), intent(in)  :: charge
    !> Electric field with its x, y, and z components
    type(Electric_Field), intent(in)  :: e_field
    real(dp) :: force_ext(3)
    
    force_ext = charge * e_field%components
  end function

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
    real(dp), intent(out) :: force_HF(n_cartesian)
    ! Local variables
    integer(i32)           :: nr
    integer(i32),parameter :: l_max = 1
    integer(i32),parameter :: lm_max = (l_max + 1)**2
    real(dp), allocatable  :: grad(:, :, :)
    real (dp), external    :: rfmtinp
    
    nr = size( radial_grid )
    call assert( size(vc_lm, 1) >= lm_max, '1st dim of vc_lm must be >= lm_max')
    call assert( size(vc_lm, 2) == nr, '2nd dim of vc_lm must contain nr elements' )
    
    allocate( grad(lm_max, nr, n_cartesian), source = 0._dp )
    call gradrfmt( 1, nr, radial_grid, lm_max, nr, vc_lm(1:lm_max, :), grad )
    force_HF = Z * grad(1, 1, :) * y00
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
    real(dp), intent(out) :: force_core(n_cartesian)
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
    allocate( grad(lm_max, nr, n_cartesian), source = 0._dp )

    rho_mt(1,:) = rho_core(:)/y00
    call gradrfmt(1, nr, radial_grid, lm_max, nr, rho_mt, grad )
    do j = 1, n_cartesian
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
    allocate( grad(lm_max, nr, n_cartesian), source = 0._dp )

    ! Valence charge = total charge - core charge (rho_core has only l=0 component)
    rho_val(1,:) = rho_val(1,:) - rho_core/y00
    call gradrfmt( l_max, nr, radial_grid, lm_max, nr, rho_val, grad )
    vaux(1,:) = vaux(1,:) + term/y00
    do j = 1, n_cartesian
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
    real(dp), intent(out)    :: F(n_cartesian)
    ! Local variables
    integer(i32)             :: i, j, m, n
    real(dp), allocatable    :: aux(:)
    complex(dp), allocatable :: diff(:, :), prod(:, :)

    m = size( psi, 1 )
    n = size( psi, 2 )
    allocate( diff(m, m), prod(m, n), aux(n) )
    call assert( size(H,1) == m, 'H must have m elements along 1st dim')
    call assert( size(H,2) == m, 'H must have m elements along 2nd dim')
    call assert( size(H,3) == n_cartesian, 'H must have n_cartesian elements along 3rd dim')
    call assert( size(S,1) == m, 'S must have m elements along 1st dim')
    call assert( size(S,2) == m, 'S must have m elements along 2nd dim')
    call assert( size(S,3) == n_cartesian, 'S must have n_cartesian elements along 3rd dim')
    call assert( size(occ) == n, 'occ must have n elements')

    ! loop over x, y, z
    do j = 1, n_cartesian
      diff = H(:,:,j) - S(:,:,j)
      call hermitian_matrix_multiply( diff, psi, prod )
      forall (i = 1:n) aux(i) = dble( dot_product( psi(:,i), prod(:,i) ) )
      F(j) = dot_product( occ, aux )
    end do
  end subroutine
end module
