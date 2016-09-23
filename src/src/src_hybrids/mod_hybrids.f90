!
!  Main module defining the variables for Hartree-Fock hybrids.
!
module mod_hybrids

    use modinput
    use modmain
    use modgw
    use modmpi, only: rank

    implicit none

    ! number of HF cycles
    integer :: ihyb

    ! non-local exchange energy
    real(8) :: exnl

    ! non-local exchange potential
    complex(8), allocatable :: vxnl(:,:,:)
    complex(8), allocatable :: vxnlcc(:,:)
    complex(8), allocatable :: bxnl(:,:,:)

    ! APW matrix elements of the non-local potential
    complex(8), allocatable :: vnlmat(:,:,:)

    ! complex(8), allocatable :: minmmat(:,:,:)
    ! complex(8), allocatable :: mincmat(:,:,:)
    ! complex(8), allocatable :: micmmat(:,:,:)
    complex(8), allocatable :: miccmat(:,:,:)

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

        if (allocated(minmmat)) deallocate(minmmat)
        if (allocated(mincmat)) deallocate(mincmat)
        if (allocated(micmmat)) deallocate(micmmat)
        if (allocated(miccmat)) deallocate(miccmat)

        return
    end subroutine

end module
