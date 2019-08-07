
subroutine init_hybrids()
    use modinput
    use modgw
    use modmpi, only : rank
    ! use mod_mpi_gw, only: myrank
    implicit none
    integer :: lmax, ik

    ! initialize GW MPI environment
    ! myrank = rank
    ! call init_mpi_gw

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

!---------------------------------------
! Parameters for the bare coulomb potential
!---------------------------------------
    if (.not.associated(input%gw%BareCoul)) &
    &  input%gw%BareCoul => getstructbarecoul(emptynode)

!---------------------------------------------------------
! Intialize auxiliary arrays used further for convenience
!---------------------------------------------------------
    call init_misc_gw()

!---------------------------------------
!   Initialize k/q grids
!---------------------------------------
    input%gw%ngridq  = input%groundstate%ngridk
    input%gw%vqloff  = input%groundstate%vkloff
    input%gw%reduceq = input%groundstate%reducek

    call init_kqpoint_set()

!--------------------------------------------------------------
! Calculate the integrals to treat the singularities at G+q->0
!--------------------------------------------------------------
    call setsingc()

! Gaunt coefficients
    lmax = max(input%groundstate%lmaxapw+1, 2*(input%gw%mixbasis%lmaxmb+1))
    call calcgauntcoef(lmax)

    return
end subroutine
