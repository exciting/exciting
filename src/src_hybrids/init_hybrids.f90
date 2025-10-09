
subroutine init_hybrids()
    use modinput
    use modgw
    use mod_bands, only: evalfv
    use modmpi, only : rank, mpiglobal
    use vx_enums, only: HYB_PBE0, HYB_HSE
    use mod_lattice, only : omega !unit cell volume
    use errors_warnings, only: terminate_if_false
    use precision, only: dp, str_16 
    use hse_singularity , only : hse_singularity_exact_solution, hse_singularity_Taylor_expansion    
#include "offload.fpp"
    implicit none
    integer :: lmax, ik
    integer :: number_k_points
    real(dp) :: omega_hse
    real(dp) :: integral_singularity_hse
    real(dp) :: volume_unit_cell
    character(str_16) :: hse_singularity_method
    ! Time for the self-energy calculations
    call init_timing()
    !---------------------------------------
    ! MB parameters are taken from GW
    !---------------------------------------
    if (.not.associated(input%gw)) &
    &  input%gw => getstructgw(emptynode)

    !---------------------------------------
    ! Options for the mixed basis functions
    !---------------------------------------
    if (.not.associated(input%gw%MixBasis)) &
    &  input%gw%MixBasis => getstructmixbasis(emptynode)
    input%gw%mixbasis%lmaxmb = input%groundstate%Hybrid%lmaxmb
    input%gw%mixbasis%epsmb  = input%groundstate%Hybrid%epsmb
    input%gw%mixbasis%gmb    = input%groundstate%Hybrid%gmb

    !---------------------------------------
    ! Parameters for the bare coulomb potential
    !---------------------------------------
    if (.not.associated(input%gw%BareCoul)) &
    &  input%gw%BareCoul => getstructbarecoul(emptynode)
    input%gw%barecoul%basis = input%groundstate%Hybrid%BasisBareCoulomb

    !---------------------------------------------------------
    ! Initialize auxiliary arrays used further for convenience
    !---------------------------------------------------------
    call init_misc_gw()

    !---------------------------------------
    ! Initialize k/q grids
    !---------------------------------------
    input%gw%ngridq  = input%groundstate%ngridk
    input%gw%vqloff  = input%groundstate%vkloff
    input%gw%reduceq = input%groundstate%reducek

    call init_kqpoint_set()

    !--------------------------------------------------------------
    ! Calculate the integrals to treat the singularities at G+q->0
    !--------------------------------------------------------------
    call terminate_if_false(mpiglobal, (xctype(1) == HYB_HSE) .or. (xctype(1) == HYB_PBE0), "The functional selected is not an hybrid functional.")
    if (xctype(1) == HYB_HSE) then
        omega_hse = input%groundstate%Hybrid%omega
        number_k_points = kqset%nkpt
        hse_singularity_method = input%groundstate%Hybrid%HSEsingularity
        if (hse_singularity_method=="Exact") then
          volume_unit_cell = omega
          integral_singularity_hse = hse_singularity_exact_solution(omega_hse, volume_unit_cell, number_k_points)  
        else if (hse_singularity_method=="Taylor") then 
          integral_singularity_hse = hse_singularity_Taylor_expansion(omega_hse)   
        endif
        singc2= integral_singularity_hse / (4.d0*pi*dble(number_k_points))
    else if (xctype(1) == HYB_PBE0) then
        call setsingc()
    end if

    ! Gaunt coefficients
    lmax = max(input%groundstate%lmaxapw+1, 2*(input%gw%mixbasis%lmaxmb+1))
    call calcgauntcoef(lmax)

    !-------------------------------------------------------------------------------
    ! Matrix block size
    !-------------------------------------------------------------------------------
    mblksiz = input%groundstate%Hybrid%mblksiz
    if (mblksiz <= 0) then
        mblksiz = 1000000 ! just a big number to account for all available states
    end if

    !-------------------------------------------------------------------------------
    ! Batching size for the computation of the expansion coefficients in the 
    ! interstitial region
    !-------------------------------------------------------------------------------
    input%gw%GBatchCount = input%groundstate%Hybrid%GBatchCount

    !----------------------------------
    ! Upload GS globals to the devices
    !----------------------------------
    OMP_OFFLOAD target enter data map(always, to: idxas, idxlo, idxlm, lorbl, apword, nlorb, corind, evalcr, evalfv)

    return
end subroutine
