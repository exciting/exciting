!> Routines to calculate the mgga or gga potential. 
module mgga_poteff
    use exciting_mpi, only: xmpi_bcast
    use precision, only: dp, i32, long_int
    use modinput, only: input
#include "asserts.fpp"
    use mgga_potxc, only: potxc_ir_spinunpolarised, potxc_mt_spinunpolarised
    use mod_potential_and_density, only: vclir, vclmt, vhalfir, vhalfmt
    use mod_muffin_tin, only: nrmtmax, nrmt, lmmaxvr, nrmtinr, lmmaxinr
    use mod_atoms, only: nspecies, natoms, idxas, natmtot
    use mod_Gvector, only: ngrtot
    use modmpi, only: mpiglobal

    implicit none 

    private 

    !> muffin-tin effective GGA potential used to construct the basis functions for meta-GGA calculations
    real(dp), public, allocatable :: veffmt_gga(:, :, :)
    !> muffin-tin GGA potential used to construct the basis functions for meta-GGA calculations
    real(dp), public, allocatable :: vxcmt_gga(:, :, :)
    !> muffin-tin GGA exchange energy density 
    real(dp), public, allocatable :: exmt_gga(:, :, :)
    !> muffin-tin GGA correlation energy density 
    real(dp), public, allocatable :: ecmt_gga(:, :, :)
    
    !> interstitial effective GGA-potential used to construct the basis functions for meta-GGA calculations
    real(dp), public, allocatable :: veffir_gga(:)
    !> interstitial GGA-potential 
    real(dp), public,  allocatable :: vxcir_gga(:)  
    !> interstitial GGA exchange energy density 
    real(dp), public, allocatable :: exir_gga(:)
    !> interstitial GGA correlation energy density 
    real(dp), public, allocatable :: ecir_gga(:)

    public :: calc_poteff_mgga, calc_poteff_gga
    public :: init_poteff_gga
    public :: mgga_mixer

    contains 

        !> Allocate all arrays needed to calculate the GGA potential
        subroutine init_poteff_gga()
            if ( allocated( veffmt_gga) ) deallocate(veffmt_gga )
            allocate(veffmt_gga(lmmaxvr, nrmtmax, natmtot))
            if ( allocated( vxcmt_gga) ) deallocate( vxcmt_gga)
            allocate(vxcmt_gga(lmmaxvr, nrmtmax, natmtot))
            if ( allocated( exmt_gga) ) deallocate(exmt_gga )
            allocate(exmt_gga(lmmaxvr, nrmtmax, natmtot))
            if ( allocated( ecmt_gga) ) deallocate(ecmt_gga )
            allocate(ecmt_gga(lmmaxvr, nrmtmax, natmtot))
            if ( allocated( veffir_gga) ) deallocate(veffir_gga )
            allocate(veffir_gga(ngrtot))
            if ( allocated( vxcir_gga) ) deallocate(vxcir_gga )
            allocate(vxcir_gga(ngrtot))
            if ( allocated( exir_gga) ) deallocate(exir_gga )
            allocate(exir_gga(ngrtot))
            if ( allocated( ecir_gga) ) deallocate(ecir_gga)
            allocate(ecir_gga(ngrtot))
        end subroutine 

        !> Calculate the effective meta-GGA potential. Returns the multiplicative potential `vxcmt_mgga` and `vxcir_mgga` and the mon-multiplicative potential 
        !> `vxcmt_mgga_nonmult` and `vxcir_mgga_nonmult`. 
        !> In the spin-polarised case will return the multiplicative magnetic potential `bxcmt_mgga` and `bxcir_mgga` abd the non-multiplicative 
        !> magnetic potential `bxcmt_mgga_nonmult` and `bxcir_mgga_nonmult`. 
        subroutine calc_poteff_mgga(veffmt_mgga, veffir_mgga, xcgrad, xctype, rhomt, rhoir, exmt_mgga, ecmt_mgga, ecir_mgga, exir_mgga, vxcmt_mgga, & 
                                    vxcmt_mgga_nonmult, vxcir_mgga, vxcir_mgga_nonmult, ked_ir, ked_mt)
            !> Effective (multiplicative) potential in the muffin-tin region 
            real(dp), intent(inout) :: veffmt_mgga(:, :, :)
            !> Effective (multiplicative) potential in the interstitial region 
            real(dp), intent(inout) :: veffir_mgga(:)
            !> Degree of exchange-correlation potential (can be 2 - for GGA or 3 - for meta-GGA)
            integer(i32), intent(in) :: xcgrad
            !> Exchange-correlation type
            integer(i32), intent(in) :: xctype(:)
            !> Density in the muffin-tin region
            real(dp), intent(in) :: rhomt(:, :, :)
            !> Density in the interstitial region
            real(dp), intent(in) :: rhoir(:)
            !> multiplicative exchange-correlation potential in the muffin-tin region
            real(dp), intent(inout) :: vxcmt_mgga(:, :, :)
            !> multiplicative exchange-correlation potential in the interstitial region
            real(dp), intent(inout) :: vxcir_mgga(:)
            !> exchange energy density in the muffin-tin region
            real(dp), intent(inout) :: exmt_mgga(:, :, :)
            !> correlation energy density in the muffin-tin region
            real(dp), intent(inout) :: ecmt_mgga(:, :, :)
            !> exchange energy density in the interstitial region
            real(dp), intent(inout) :: exir_mgga(:)
            !> correlation energy density in the interstitial region
            real(dp), intent(inout) :: ecir_mgga(:)
            !> kinetic energy density in the muffin-tin region
            real(dp), intent(in) :: ked_mt(:, :, :)
            !> kinetic energy density in the interstitial region
            real(dp), intent(in) :: ked_ir(:)
            !> mon-multiplicative exchange-correlation potential in the muffin-tin region
            real(dp), intent(inout) :: vxcmt_mgga_nonmult(:, :, :)
            !> mon-multiplicative exchange-correlation potential in the interstitial region
            real(dp), intent(inout) :: vxcir_mgga_nonmult(:)

            CALL_ASSERT((xcgrad == 3), message='subroutine for mGGAs only.')
            
            call potxc_mt_spinunpolarised(xcgrad, xctype, rhomt, vxcmt_mgga, exmt_mgga, ecmt_mgga, ked_mt, vxcmt_mgga_nonmult)
            call potxc_ir_spinunpolarised(xcgrad, xctype, rhoir, vxcir_mgga, exir_mgga, ecir_mgga, ked_ir, vxcir_mgga_nonmult)

            call symrf(1, vxcmt_mgga, vxcir_mgga)
            call symrf(1, vxcmt_mgga_nonmult, vxcir_mgga_nonmult)

            call potcoul
            call add_coulomb_and_xc_potentials(veffmt_mgga, veffir_mgga, vclmt, vclir, vxcmt_mgga, vxcir_mgga)
        end subroutine 

        !> Calculate the effective GGA potential. Returns the multiplicative potential `vxcmt_gga' and `vxcir_gga` and the mon-multiplicative potential.
        subroutine calc_poteff_gga(veffmt_gga, veffir_gga, xcgrad, xctype, rhomt, rhoir, vxcmt_gga, vxcir_gga, exmt_gga, ecmt_gga, ecir_gga, exir_gga)
            !> Effective (multiplicative) potential in the muffin-tin region (MTR)
            real(dp), intent(inout) :: veffmt_gga(:, :, :)
            !> Effective (multiplicative) potential in the interstitial region (IR)
            real(dp), intent(inout) :: veffir_gga(:)
            !> degree of exchange-correlation potential (can be 2 - for GGA or 3 - for meta-GGA)
            integer(i32), intent(in) :: xcgrad
            !> exchange-correlation type
            integer(i32), intent(in) :: xctype(:)
            !> density in the muffin-tin region 
            real(dp), intent(in) :: rhomt(:, :, :)
            !> density in the interstitial region 
            real(dp), intent(in) :: rhoir(:)
            !> exchange-correlation potential in the muffin-tin region
            real(dp), intent(inout) :: vxcmt_gga(:, :, :)
            !> exchange-correlation potential in the interstitial region
            real(dp), intent(inout) :: vxcir_gga(:)
            !> exchange energy density in the muffin-tin region
            real(dp), intent(inout) :: exmt_gga(:, :, :)
            !> correlation energy density in the muffin-tin region
            real(dp), intent(inout) :: ecmt_gga(:, :, :)
            !> exchange energy density in the interstitial region
            real(dp), intent(inout) :: exir_gga(:)
            !> correlation energy density in the interstitial region
            real(dp), intent(inout) :: ecir_gga(:)

            CALL_ASSERT((xcgrad == 2), message='Only with libxc GGAs implemented.')
            
            call potxc_mt_spinunpolarised(xcgrad, xctype, rhomt, vxcmt_gga, exmt_gga, ecmt_gga)
            call potxc_ir_spinunpolarised(xcgrad, xctype, rhoir, vxcir_gga, exir_gga, ecir_gga)

            call symrf(1, vxcmt_gga, vxcir_gga)

            call potcoul
            call add_coulomb_and_xc_potentials(veffmt_gga, veffir_gga, vclmt, vclir, vxcmt_gga, vxcir_gga)

        end subroutine calc_poteff_gga

        !> Adds the exchange correlation potential `vxcmt` and the Coulomb potential `vclmt`
        !> in the muffin-tin region to obtain the effective potential `veffmt`. 
        !> The same is done in the interstital region to obtain the effective potential `veffir`. 
        subroutine add_coulomb_and_xc_potentials(veffmt, veffir, vclmt, vclir, vxcmt, vxcir)
            use constants, only: y00
            !> Effective potential in the muffin-tin region
            real(dp), intent(inout) :: veffmt(:, :, :)
            !> Effective potential in the interstitial region 
            real(dp), intent(inout) :: veffir(:)
            !> Coulomb potential in the muffin-tin region 
            real(dp), intent(inout) :: vclmt(:, :, :)
            !> Coulomb potential in the interstitial region  
            real(dp), intent(inout) :: vclir(:)
            !> Exchange-correlation potential in the muffin-tin region 
            real(dp), intent(inout) :: vxcmt(:, :, :)
            !> Exchange-correlation potential in the interstitial region 
            real(dp), intent(inout) :: vxcir(:)
    
            integer(i32) :: is, ias, ia, lmmax, ir, lm 
            real(dp) :: shift 
            logical :: dfthalf_on
            
            shift=input%groundstate%energyref
            dfthalf_on = associated( input%groundstate%dfthalf )
            if( dfthalf_on ) dfthalf_on = .not. input%groundstate%dfthalf%NSCF

            vclmt(1,:,:) = vclmt(1,:,:)+shift/y00
            Do is = 1, nspecies
                Do ia = 1, natoms (is)
                    ias = idxas (ia, is)
                    lmmax = lmmaxinr
                    Do ir = 1, nrmt (is)
                    If (ir .Gt. nrmtinr(is)) lmmax = lmmaxvr
                    Do lm = 1, lmmax
                        if ( dfthalf_on ) then
                            veffmt(lm,ir,ias) = vclmt(lm,ir,ias) + vxcmt(lm,ir,ias) + vhalfmt (lm, ir, ias)
                        else
                            veffmt(lm,ir,ias) = vclmt(lm,ir,ias) + vxcmt(lm,ir,ias)
                        endif
                    End Do
                    Do lm = lmmax + 1, lmmaxvr
                        veffmt(lm,ir,ias) = 0.d0
                    End Do
                    End Do
                End Do
            End Do
        
            ! interstitial part
            vclir(:) = vclir(:) + shift
            
            if ( dfthalf_on ) then
                veffir(:) = vclir(:) + vxcir(:) + vhalfir(:)
            else
                veffir(:) = vclir(:) + vxcir(:)
            endif
        end subroutine 

        !> Main routine to mix the meta-GGA potential. In the first scf-iteration `iscl = 1` the mixing routines 
        !> are initialized. If `iscl > 1`, the multiplicative and non-multiplicative potentials for both muffin-tin 
        !> and interstitial region are mixed. 
        subroutine mgga_mixer(iscl, nu, mode, currentconvergence, vcurrentconvergence)
            use modmpi,         only: rank
            use modinput
            use mod_spin,       only: ndmag
            use mod_muffin_tin, only: lmmaxvr, nrmtmax
            use mod_Gvector,    only: ngrtot
            use mod_atoms,      only: natmtot

            !> scf-iteration 
            integer(i32), intent(in)                 :: iscl
            !> (un)packed function
            real(dp), allocatable, intent(inout):: nu(:)
            !> mode for mixing: -1: call initialisation routines, -2: call destructor, else ignore
            integer(i32), intent(inout)             :: mode
            !> Convergence of current scf-iteration
            real(dp), intent(inout)            :: currentconvergence
            !> Array to save consecutive convergence values 
            real(dp), intent(inout)            :: vcurrentconvergence(:)

            integer(long_int) :: n
            integer(i32) :: id

            if (input%groundstate%mixernumber == 3) then
                write(*,*) "Error(mgga mixer): meta-GGA not tested with Pulay mixer"
                stop
            end if

            if (iscl == 1) then
                n = 2 * (lmmaxvr*nrmtmax*natmtot + ngrtot)
                if (associated(input%groundstate%spin)) n = n * (1 + ndmag)

                allocate(nu(n))
                mode = -1

                call mgga_pot_mixpack(.true., n, nu)
                if (rank == 0) then
                    call mixerifc(input%groundstate%mixernumber, n, nu, currentconvergence, mode, iscl)
                end if
                call mgga_pot_mixpack(.false., n, nu)

            else
                call mgga_pot_mixpack(.true., n, nu)

                if (rank == 0) then
                    call mixerifc(input%groundstate%mixernumber, n, nu, currentconvergence, mode, iscl)

                    do id = 1, input%groundstate%niterconvcheck - 1
                        vcurrentconvergence(id) = vcurrentconvergence(id + 1)
                    end do
                    vcurrentconvergence(input%groundstate%niterconvcheck) = currentconvergence
                end if

                call xmpi_bcast(mpiglobal, nu)

                call mgga_pot_mixpack(.false., n, nu)
            end if
        end subroutine mgga_mixer

        !> This subroutine (un)packs the mgga potentials (multiplicative and non-mult) for the mixing.
        subroutine mgga_pot_mixpack( tpack, n, nu )
            use modinput  
            use mixer_pack, only: pack_fun, unpack_fun
            use mod_potential_and_density, only: veffmt, veffir
            use mgga_potxc, only: vxcir_mgga_nonmult, vxcmt_mgga_nonmult
            !> `.true.` for packing and `.false.` for unpacking
            logical, intent(in) :: tpack
            !> number of elements per function
            integer(long_int), intent(out) :: n
            !> (un)packed function
            real(dp), intent(inout) :: nu(*)
            
            integer(i32) :: idm

            n = 0
            if (tpack) then 
                call pack_fun( veffmt, input%groundstate%lmaxvr, 1, veffir, n, nu)
                call pack_fun( vxcmt_mgga_nonmult, input%groundstate%lmaxvr, 1, vxcir_mgga_nonmult, n, nu)
            else 
                call unpack_fun( veffmt, input%groundstate%lmaxvr, 1, veffir, n, nu)
                call unpack_fun( vxcmt_mgga_nonmult, input%groundstate%lmaxvr, 1, vxcir_mgga_nonmult, n, nu)
            end if 
        end subroutine 

end module 