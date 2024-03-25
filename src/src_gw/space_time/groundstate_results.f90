!> Hold all ground state results to be used for GW calculation.
module groundstate_results
    use precision, only: dp
    use mod_kpointset, only: k_set
    implicit none
    private
    
    type, public :: groundstate_results_type
      real(dp), allocatable :: eigenvalues(:, :)               !! First-variation KS eigenvalues
      real(dp), allocatable :: occupations(:, :)               !! First-variation Occupations
      real(dp) :: efermi                                       !! Fermi level
      integer  :: n_spin                                       !! Spin occupancy per state 
      integer  :: n_atoms                                      !! Total number of atoms
      real(dp) :: lattice_vectors(3, 3)                        !! Real-space lattice vectors
      integer, allocatable :: atomic_species(:)                !! Atomic species array
      real(dp), allocatable :: atomic_positions(:, :)          !! Atomic position in Bohr unit
      ! real(dp), allocatable :: atomic_position_lat(:, :)       !! Atomic position in lattice coord.
      real(dp), allocatable :: kgrid_cart(:,:)                 !! k-grid in cartesian coord.
      ! real(dp), allocatable :: kgrid_lat(:,:)                  !! k-grid in lattice coord.
      real(dp), allocatable :: kgrid_weight(:)                 !! weight of k-grid points

    contains
      procedure :: init_from_globals   !! Initialise the ground state results from globals
      procedure :: n_ks_states         !! Get the number of KS states == number of bands
      procedure :: n_kpoints           !! Get the number of k-points (reducible or irreducible?)
      procedure :: n_occupied          !! Number of occupied states
    end type
    
    contains
    
    ! TODO(Alex/Man).
    ! See how existing GW sets the fermi level and shifts the KS energies
    ! - we should be consistent.
    
    ! TODO(Alex/Man).
    ! Add additional, relevant data.

    !> Encapsulation all global quantities from an exciting ground state calculation.
    !> 
    !>  Note, on a second pass, where data is defined from input file
    !>  we should assign from that rather than the global variables.
    !>
    !> Detail how each quantity is obtained
    !> --------------------------------------
    !>  evalfv. call init_dft_eigenvalues()
    !>    Reads first variation KS eigenvalues from EVALFV.OUT, which is a binary file.
    !>    Computes efermi and shifts KS eigenvalues such that efermi = 0.
    !>  efermi.
    !>     ADD ME
    !>  occfv.
    !>     ADD ME
    !>  nspnfv. 
    !>    Set to 1 or 2 in init0, according to `isspinspiral()`.
    !>    I actually think associated(input%groundstate%spin) should be used.
    !>  
    subroutine init_from_globals(this)
        use mod_bands, only: evalfv, occfv
        use mod_eigenvalue_occupancy, only: efermi
        use mod_spin, only: nspnfv
        use mod_lattice, only: avec
        use mod_atoms, only: nspecies, natoms, natmtot, atposc 
        use mod_misc_gw, only: atposl
        type(k_set) :: kset
        class(groundstate_results_type), intent(inout) :: this

        !> Local variables 
        integer :: is, ia, i_atoms
        
        ! Energies
        this%eigenvalues = evalfv
        this%efermi = efermi
        ! Shifting Fermi energy to zero 
        this%eigenvalues = this%eigenvalues - efermi

        ! Occupations
        this%n_spin = nspnfv
        this%occupations = occfv

        ! Structure
        this%lattice_vectors = avec

        ! Radial mesh info could be required

        !> Index to atoms and species "1st index->(atom), 2nd index->(species)"
        ! this%atom_species(:,:) = idxas(:,:)

        !> Atomic position : indexing->(3, n_atoms)
        this%n_atoms = natmtot
        allocate(this%atomic_species(this%n_atoms))

        i_atoms = 0
        do is = 1, nspecies
            do ia = 1, natoms(is)
                i_atoms = i_atoms + 1
                this%atomic_species(i_atoms) = is
                ! needs to check if atposc(:,:,:) in Bohr or in other units
                ! if not then convert to Bohr and proceed 
                this%atomic_positions(:, i_atoms) = atposc(:, ia, is)
            enddo
        enddo
        print*, 'space-time testing MH'

        !  APW info
        this%kgrid_cart(:,:) = kset%vkc(:,:)
        ! this%kgrid_lat(:,:) = kset%vkl(:,:)
        this%kgrid_weight(:) = kset%wkpt(:)
        !
        ! TO ADD from init2
        !  Possibly k+q set. See what existing GW initialisation
        !  retains.

    end subroutine

    !> Class destructor to free memory.
    subroutine finalize(this)
        class(groundstate_results_type), intent(inout) :: this
        deallocate(this%eigenvalues)
        deallocate(this%occupations)
        deallocate(this%atomic_species)
        deallocate(this%atomic_positions)
        deallocate(this%kgrid_cart)
        deallocate(this%kgrid_weight)
    end subroutine

    !> Get the number of Kohn-Sham states (bands).
    !> Also see `nstfv`
    integer function n_ks_states(this)
        class(groundstate_results_type), intent(inout) :: this
        n_ks_states = size(this%eigenvalues, 1) * this%n_spin
    end function

    !> Get the number of occupied states at zero temperature.
    !> Assumes NO partial occupations.
    !>
    ! TODO(Alex) See mod_bands, and how these variables are used in `init_dft_eigenvalues.f90`
    ! It could replace what is implemented below.
    integer function n_occupied(this)
        class(groundstate_results_type), intent(inout) :: this
        !TODO(Alex/Man) Implement!
        n_occupied = -1
    end function
    
    !> Get the number of k-points (reducible or irreducible?)
    ! TODO(Alex) Usually named kset%nkpt. Check init_kqpoint_set and see what quantity it is (assume reduced grid)
    integer function n_kpoints(this)
       class(groundstate_results_type), intent(inout) :: this
       n_kpoints = size(this%eigenvalues, 2)
    end function

end module
