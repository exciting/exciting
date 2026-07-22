module mgga_init 
#include "asserts.fpp"
   use precision, only: dp

    implicit none 

    private 
    
    !> mgga exchange correlation type 
    integer, public :: xctype_mgga(3)
    !> number of GGA iterations 
    logical, public :: mgga_read_in
    !> mgga exchange-correlation functional description
    character (len=512), public :: xcdescr_mgga

    public :: init_mgga, set_mgga_potential
    
    contains

        !> Initializes all routines required for a meta-GGA calculation. 
        subroutine init_mgga()
            use kinetic_energy_density_vars, only: ked_var_init
            use kinetic_energy_density, only: allocate_ked
            use mgga_poteff, only: init_poteff_gga
            use mgga_potxc, only:init_potxc_non_mult_mgga

            call ked_var_init
            call allocate_ked
            call init_poteff_gga
            call init_potxc_non_mult_mgga
        end subroutine 

        !> Configures the `xctype_mgga` array based on the provided input file,
        !> ensuring proper handling of exchange, correlation, or combined XC functionals.
        subroutine set_mgga_potential()
            use modinput

            CALL_ASSERT(associated(input%groundstate%mgga),  message="Meta-GGA potential setup requires meta-GGA parameters in the input file.")

            xctype_mgga(1) = 100 
            input%groundstate%xctype = "LibXC"  

            xctype_mgga(2) = input%groundstate%mgga%exchangenumber

            ! Assign LDA-correlation functional in case of TASK 
            if (xctype_mgga(2) == 707) then
                xctype_mgga(3) = 12 
            else
                xctype_mgga(3) = input%groundstate%mgga%correlationnumber
            end if

            if (input%groundstate%mgga%xcnumber .ne. 0) then
                xctype_mgga(2) = input%groundstate%mgga%xcnumber
                xctype_mgga(3) = 0  ! Combined XC functional; correlation is already included
            end if
        end subroutine set_mgga_potential
 
end module 
