
subroutine init_hybrids()
    use modinput
    use modgw
    use modmpi, only : rank
    implicit none
    integer :: lmax, ik

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
    if (xctype(1) == 408) then
        ! singular term in HSE
        singc2 = pi/input%groundstate%Hybrid%omega**2 / (4.d0*pi*dble(kqset%nkpt))
    else
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

    return
end subroutine
