!
!  Main module defining the variables for Hartree-Fock hybrids.
!
module mod_hybrids

    use modmain
    use modgw
    use mod_coulomb_potential, only: delete_coulomb_potential
    use mod_misc_gw, only : gammapoint
    use mod_kpointset, only: delete_Gk_vectors, delete_k_vectors, delete_kq_vectors, delete_G_vectors
    use modmpi, only: rank
#include "offload.fpp"

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

        ! Deallocate all the reciprocal space meshes
#if defined(FLANG_OPENMP_DERIVED_TYPE_MAP_BUG_WORKAROUND)
        OMP_OFFLOAD target exit data map(delete: Gset%ivg, Gkqset%igkig, Gset%intgv, Gset%ivgig, Gqbarc%igigk)
#endif
        OMP_OFFLOAD target exit data map(delete: kset, Gset, Gkset, Gkqset, Gqset, Gqbarc, kqset)
       
        call delete_k_vectors(kset)
        call delete_G_vectors(Gset)
        call delete_Gk_vectors(Gkset)
        call delete_Gk_vectors(Gkqset)
        call delete_Gk_vectors(Gqset)
        call delete_Gk_vectors(Gqbarc)
        call delete_kq_vectors(kqset)

        nullify(input%gw%MixBasis)
        nullify(input%gw%BareCoul)
        nullify(input%gw)
        call rereadinput()

        return
    end subroutine

end module
