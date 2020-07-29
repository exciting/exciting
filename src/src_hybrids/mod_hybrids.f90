!
!  Main module defining the variables for Hartree-Fock hybrids.
!
module mod_hybrids

    use modmain
    use modgw
    use mod_coulomb_potential, only: delete_coulomb_potential
    use mod_misc_gw, only : gammapoint
    use modmpi, only: rank

    implicit none

    ! set true if HF-Hybrids are used as starting point
    Logical :: hybridhf
    data hybridhf / .false. /

    ! non-local exchange energy
    real(8) :: exnl

    ! non-local exchange potential
    complex(8), allocatable :: vxnl(:,:,:)

    ! APW matrix elements of the non-local potential
    complex(8), allocatable :: vnlmat(:,:,:)

    ! File names
    character(80) :: fname_vxnl
    data fname_vxnl / 'VXNL.OUT' /
    character(80) :: fname_vxnlmat
    data fname_vxnlmat / 'VXNLMAT.OUT' /

!*******************************************************************************
contains

    ! deallocate hybrids related data
    subroutine exit_hybrids()

        ! deallocate global
        if (allocated(vxnl)) deallocate(vxnl)
        if (allocated(vnlmat)) deallocate(vnlmat)

        ! deallocate mixed-basis stuff
        call delete_product_basis
        call delete_core_states

        nullify(input%gw%MixBasis)
        nullify(input%gw%BareCoul)
        nullify(input%gw)
        call rereadinput()

        return
    end subroutine

end module
