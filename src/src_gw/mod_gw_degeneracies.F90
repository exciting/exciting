!> For the self-energy any degenerate complexes of states must be summed over
!> in order to prevent symmetry breakdown. There are two reasons for that:
!> 1) For degenerate states what is invariant to symmetry operations is the
!>    whole subspace not just a part of it. Therefore, all the whole degenerate subspace must be summed over [1,2]. 
!> 2) The limited precision in the wavefunctions from the ground state calculation can cause symmetry breakdown for the self-energy [2].
!>
!> References:
!> [1] https://journals.aps.org/prb/abstract/10.1103/PhysRevB.35.5585
!> [2] https://doi.org/10.1016/j.cpc.2011.12.006
!>
!> Overall, this module contains procedures and variables to properly tackle degeneracies
!> Notice that all public variables are protected (i.e. read-only for the outer-module world)
module mod_gw_degeneracies

    use precision, only: i32, dp
    use modmpi, only: terminate
#include "offload.fpp"

    implicit none 

    private
    public :: get_degenerate_limits_qp_interval_ikp, &
              enforce_degeneracy, &
              initialize_degeneracy_module, &
              delete_degeneracy_module, &
              absolute_tolerance_gw_degeneracy, &
              relative_tolerance_gw_degeneracy, &
              degenerate_subspaces, &
              band_degeneracy, &
              ibgw_including_degeneracy, &
              nbgw_including_degeneracy

    !> Absolute tolerance for the degenerate states
    real(dp), protected :: absolute_tolerance_gw_degeneracy
    !> Relative tolerance for the degenerate states
    real(dp), protected :: relative_tolerance_gw_degeneracy
    !> Lower and upper indexes, as well size of the degenerated subspaces, for each irreducible point
    !> This list of the degenerate subspaces is as follows: [each subspace info is
    !> giving an individual rows providing the starting [row 1] and ending [row 2]
    !> indices for each of them, as well their size [row 3].
    integer(i32), protected, allocatable :: degenerate_subspaces(:,:,:)
    !> Lower state id to which the QP corrections are computed taking into account lower degenerate states with ibgw 
    !> If not taken into account the result will be dependent on the window because of the averaging
    integer(i32), protected :: ibgw_including_degeneracy
    !> Upper state id to which the QP corrections are computed taking into account upper degenerate states with nbgw 
    !> If not taken into account the result will be dependent on the window because of the averaging
    integer(i32), protected :: nbgw_including_degeneracy
    !> List providing the degeneracy of a band
    integer(i32), protected, allocatable :: band_degeneracy(:,:)

interface enforce_degeneracy
    module procedure :: enforce_degeneracy_real_dp, &
                        enforce_degeneracy_complex_dp
end interface

contains

    !> Release module state so a subsequent GW initialization can start cleanly.
    subroutine delete_degeneracy_module()
        if (allocated(degenerate_subspaces)) deallocate(degenerate_subspaces)
        if (allocated(band_degeneracy)) then
            OMP_OFFLOAD target exit data map(delete: band_degeneracy)
            deallocate(band_degeneracy)
        end if
    end subroutine delete_degeneracy_module
    
    ! Let the compiler to automatically generate enforce_degeneracy(ikp, input_vector) for different types
    ! See enforce_degeneracy_template.inc: The function name is concatenated by _TYPE_PRECISION.
    !
    ! Template inputs:
    ! ikp: irreducible k-point index.
    ! input_vector: data at ikp where degeneracy must be strictly enforced.
    !
    ! Generating enforce_degeneracy_real_dp(ikp, input_vector)
#define TYPE1 real
#define PRECISION1 dp
#include "enforce_degeneracy_template.inc"
    ! Generating enforce_degeneracy_complex_dp(ikp, input_vector)
#define TYPE1 complex
#define PRECISION1 dp
#include "enforce_degeneracy_template.inc"


    !> Obtains the limits for degenerate subspaces for the giving 
    !> irreducible point. Notice that degenerate_subspaces = -1 is defined
    !> in this array as not defined values. See init_dft_eigenvalues.f90 for more information
    !> Also, note that any subspace including the limit bands is explicitly treated
    !> with all their members.  
    pure subroutine get_degenerate_limits_qp_interval_ikp(ikp, init_space, final_space)
            
        ! Using the intervals for which the QP corrections are computed, i.e. [ibgw, nbgw].
        use modgw, only: ibgw, nbgw

        implicit none

        !> Id of the irreducible k-point
        integer(i32), intent(in) :: ikp 
        !> Index of the degenerate subspace containing ibgw
        integer(i32), intent(out) :: init_space
        !> Index of the degenerate subspace containing nbgw
        integer(i32), intent(out) :: final_space 

        init_space  = minloc(degenerate_subspaces(2,:,ikp), dim=1, &
                            mask = degenerate_subspaces(2,:,ikp) >= ibgw)
        final_space = maxloc(degenerate_subspaces(2,:,ikp), dim=1, &
                        mask = degenerate_subspaces(1,:,ikp) <= nbgw)

    end subroutine get_degenerate_limits_qp_interval_ikp

    !> Initialize the module protected variables.
    subroutine initialize_degeneracy_module(kset, nbands, evalfv, occfv, enforce_degeneracy_flag)
        use modinput, only: input
        use math_utils, only: get_degeneracies
        use mod_kpointset, only: k_set
        use modgw, only: ibgw, nbgw

        implicit none

        !> k reduced mesh
        type(k_set), intent(in) :: kset
        !> number of bands in the NSCF calculation
        integer(i32), intent(in) :: nbands
        !> Eigenenergies
        real(dp), intent(inout) :: evalfv(:,:)
        !> Occupations
        real(dp), intent(inout) :: occfv(:,:)
        !> Flag indicating if we enforce the degeneracy or not
        logical, intent(in) :: enforce_degeneracy_flag

        ! Id of the irreducible k-point
        integer(i32) :: ikp, ib
        ! Indexes for the degenerate space to consider when computing the QP correction 
        integer(i32) :: init_space, final_space
        ! Index to iterate over degenerate spaces
        integer(i32) :: ispace
        ! Degenerate subspaces per irreducible k-point
        integer(i32), allocatable :: degenerate_subspaces_ikp(:,:)
        ! Big negative number
        real(dp), parameter :: big_negative_number = -1.0e+42_dp

        ! Obtain the toleraces from the input
        absolute_tolerance_gw_degeneracy = input%gw%degeneracyAbsoluteTolerance
        relative_tolerance_gw_degeneracy = input%gw%degeneracyRelativeTolerance

        ! In case degeneracy is not enforced set tolerance to negative value; so that 
        ! no degeneracies are found then
        if (.not. enforce_degeneracy_flag) then
            absolute_tolerance_gw_degeneracy = big_negative_number
            relative_tolerance_gw_degeneracy = 0.0_dp 
        end if

        ! Init values for searching smaller and bigger values in the irreducible wedge
        ibgw_including_degeneracy =  huge(1_i32)
        nbgw_including_degeneracy = -huge(1_i32)

        allocate(band_degeneracy(nbands,kset%nkpt), source=0_i32)
        ! Notice that this is initialized to the maximum possible size for all the k-points,
        ! that is the number of bands in the NSCF calculation. The default value is -1; that 
        ! means that if in a given k-point there are degenerate states, there will be entries 
        ! with -1 in the first dimension. In other words we will have uninitalized values represented
        ! with -1 because the space contains degeneracies and it is thus smaller that the number of bands in the NSCF calculation.
        allocate(degenerate_subspaces(3, nbands, kset%nkpt), source=-1_i32)
        
        do ikp = 1, kset%nkpt
            degenerate_subspaces_ikp = get_degeneracies(evalfv(:nbands, ikp), absolute_tolerance_gw_degeneracy, relative_tolerance_gw_degeneracy)
            degenerate_subspaces(:,:size(degenerate_subspaces_ikp,2),ikp) =  degenerate_subspaces_ikp(:,:)

            call get_degenerate_limits_qp_interval_ikp(ikp, init_space, final_space)
            
            ! Use those to obtain ibgw and nbgw that consider degeneracies
            ! This is important as otherwise the correction can cut degenerate states 
            ! making result dependent on the choice of ibgw and nbgw when band averaging for
            ! degenerate subspaces (Note that is global among all the k-points)
            ibgw_including_degeneracy = min(degenerate_subspaces(1,init_space,ikp),ibgw_including_degeneracy)
            nbgw_including_degeneracy = max(degenerate_subspaces(2,final_space,ikp),nbgw_including_degeneracy)

            ! Here we are correcting the eigenenergies, so we set two degenerate bands to
            ! exactly the same energy value so the bands are clean of numerical noise
            ! and the results can be safely used for self-consistent approaches which are
            ! rather sensitive to small diferences.
            call enforce_degeneracy(ikp, evalfv(:nbands,ikp))
            call enforce_degeneracy(ikp, occfv(:nbands,ikp))

            ! Check the degeneracy of each band
            do ispace = 1, size(degenerate_subspaces_ikp,2)
                band_degeneracy(degenerate_subspaces_ikp(1,ispace):degenerate_subspaces_ikp(2,ispace),ikp) = degenerate_subspaces_ikp(3,ispace)
            end do

        end do

        OMP_OFFLOAD target enter data map(always, to: band_degeneracy)

    end subroutine initialize_degeneracy_module

end module mod_gw_degeneracies
