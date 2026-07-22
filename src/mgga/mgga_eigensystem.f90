!> Module sets up meta-GGA Hamiltonian and overlap matrices for both interstitial and muffin-tin region. 
module mGGA_eigensystem
    use precision, only: dp
    use kinetic_energy_density_vars
    use matrix_elements
    use modmpi, only: terminate_if_false
#include "asserts.fpp"
    use mod_Gkvector
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    use mod_muffin_tin, only: lmmaxapw, nrmt, nrmtmax
    use constants, only: zzero, zone, y00
    use physical_constants, only: alpha
    use modinput
    use mod_APW_LO, only: apwfr, lofr, nlotot, apwordmax
    use muffin_tin_basis, only: mt_basis_type
    use mgga_poteff, only: veffmt_gga, veffir_gga
    use mgga_potxc, only: vxcir_mgga_nonmult

    implicit none 

    private 

    !> Hamiltonian matrix
    complex(dp), allocatable, public :: mGGA_H(:, :)
    !> Overlap matrix
    complex(dp), allocatable, public :: mGGA_S(:, :)

    !> radial muffin-tin integrals times Gaunt coefficients for overlap matrix \(S\)
    complex(dp), allocatable :: Smat_mt_basis(:,:,:)
    !> radial muffin-tin integrals times Gaunt coefficients for Hamiltonian matrix \(H\)
    complex(dp), allocatable :: Hmat_mt_basis(:,:,:)
    
    !> characteristic function in reciprocal space
    complex(dp), allocatable :: cfun_ig(:)
    !> interstitial effective potential times characteristic function in reciprocal space
    complex(dp), allocatable :: pot_cfun_ig(:)
    !> interstitial (scalar relativistic) kinetic energy times characteristic function in reciprocal space
    complex(dp), allocatable :: kin_cfun_ig(:)
    !> interstitial non-multiplicative mgga potential times characteristic function in reciprocal space
    complex(dp), allocatable :: pot_mgga_cfun_ig(:)
    !> interstitial non-multiplicative mgga potential times gradient of characteristic function in reciprocal space
    complex(dp), allocatable :: pot_mgga_dcfun_ig(:)
    
    public :: mGGA_eig_init, mGGA_eig_free, gen_mGGA_H_and_S

    contains 


        !> This subroutine initializes variables for the calculation of matrix elements
        !> that remain constant during the entire calculation.
        subroutine mGGA_eig_init(pot_mt, pot_ir, pot_mt_mgga, pot_ir_mgga )            
            real(dp), optional, intent(in) :: pot_ir_mgga(:)
            real(dp), optional, intent(in) :: pot_mt_mgga(:, :, :)
            real(dp), allocatable :: kin_ir(:)
            real(dp), intent(in) :: pot_ir(:)
            real(dp), intent(in) :: pot_mt(:, :, :)
            complex(dp) :: grad(3, 3)
            complex(dp), allocatable :: zfft(:), dcfun_ir(:,:,:), dcfun_ig(:,:,:)
            integer :: i, ia, is, ias, ig

            ! Update MT basis 
            ked_mt_basis%apw_rad_fun = apwfr     
            ked_mt_basis%lo_rad_fun = lofr

            call me_init( ked_mt_basis, input%groundstate%lmaxvr, ked_Gset )

            ! compute MT radial integrals times Gaunt coefficients
            if ( present(pot_mt_mgga) ) then 
                call gen_overlap_hamiltonian_mt_basis( Smat_mt_basis, Hmat_mt_basis, pot_mt, pot_mt_mgga=pot_mt_mgga, &
                                                    lmax_apw=input%groundstate%lmaxapw,  lmax_pot=input%groundstate%lmaxvr, kinetic=.true., &
                                                    potential=.true., potential_mgga=.true.)
            else 
               call gen_overlap_hamiltonian_mt_basis( Smat_mt_basis, Hmat_mt_basis, pot_mt, lmax_apw=input%groundstate%lmaxapw, &
                                                    lmax_pot=input%groundstate%lmaxvr, kinetic=.true., potential=.true., potential_mgga=.false.)
            end if 
            
            ! compute interstitial representation of characteristic function, 
            ! effective potential and kinetic energy
            call me_ir_alloc( pot_cfun_ig )
            call me_ir_prepare( zone, pot_ir, zzero, pot_cfun_ig )
            
            if ( present(pot_ir_mgga) ) then 
                 call me_ir_alloc( pot_mgga_cfun_ig )
                 call me_ir_prepare( zone, 0.5_dp * pot_ir_mgga, zzero, pot_mgga_cfun_ig )
            end if 

            if (allocated(cfun_ig)) deallocate(cfun_ig) 
            allocate( cfun_ig(ked_Gset%ngvec) ) 
            allocate( kin_ir(ked_Gset%ngrtot) )
            call gencfunig( ked_Gset%ngvec, ked_Gset%gc, ked_Gset%vgc, cfun_ig )
            select case( input%groundstate%ValenceRelativity )
                case( 'iora*' )
                call terminate_if_false( .false., '(mgga_eig_init) &
                   & metaGGAs in combination with `iora*` valence relativity not implemented.' )
                case( 'none' )
                kin_ir = 0.5_dp
                case default
                kin_ir = 0.5_dp / (1.0_dp - 0.5_dp * alpha**2 * veffir_gga )
            end select
            call me_ir_alloc( kin_cfun_ig )
            call me_ir_prepare( zone, kin_ir, zzero, kin_cfun_ig )
                    
        end subroutine mGGA_eig_init

        !> This subroutine frees memory from the module variables
        !> and cleans up the matrix elements module.
        subroutine mGGA_eig_free
            if( allocated( Smat_mt_basis ) ) deallocate( Smat_mt_basis )
            if( allocated( Hmat_mt_basis ) ) deallocate( Hmat_mt_basis )
            if( allocated( pot_cfun_ig ) ) deallocate( pot_cfun_ig )
            if( allocated( kin_cfun_ig ) ) deallocate( kin_cfun_ig )
            if( allocated( pot_mgga_cfun_ig ) ) deallocate( pot_mgga_cfun_ig )
            if( allocated( pot_mgga_dcfun_ig ) ) deallocate( pot_mgga_dcfun_ig )
        end subroutine mGGA_eig_free


        !> This subroutine computes the radial muffin-tin integrals times Gaunt coefficients
        !> for the overlap and Hamiltonian matrix.
        subroutine gen_overlap_hamiltonian_mt_basis( Smat_mt_basis, Hmat_mt_basis, pot_mt, pot_mt_mgga, lmax_apw, lmax_pot, kinetic, potential, potential_mgga )
            !> overlap radial muffin-tin integrals times Gaunt coefficients
            complex(dp), allocatable, intent(inout) :: Smat_mt_basis(:,:,:)
            !> Hamiltonian radial muffin-tin integrals times Gaunt coefficients
            complex(dp), allocatable, intent(inout) :: Hmat_mt_basis(:,:,:)
            !> muffin-tin effective potential as real spherical harmonics expansion
            real(dp), intent(in) :: pot_mt(:,:,:)
            !> non-multiplicative meta-GGA muffin-tin effective potential as real spherical harmonics expansion
            real(dp), optional, intent(in) :: pot_mt_mgga(:,:,:)
            !> maximum angular momentum \(l\) for APWs and potential expansion (default: from input file)
            integer, optional, intent(in) :: lmax_apw, lmax_pot
            !> include kinetic energy contribution (default: true)
            logical, optional, intent(in) :: kinetic
            !> include effective potential contribution (default: true)
            logical, optional, intent(in) :: potential
            !> include non multiplicative potential 
            logical, optional, intent(in) :: potential_mgga

            integer :: lmaxapw, lmaxpot, is, ia, ias, i, ip
            real(dp), allocatable :: rfun(:,:)
            logical :: kin, pot, pot_mgga

            lmaxapw = ked_lmaxapw
            if( present( lmax_apw ) ) lmaxapw = lmax_apw
            lmaxpot = ked_lmaxvr
            if( present( lmax_pot ) ) lmaxpot = lmax_pot
            kin = .true.
            if( present( kinetic ) ) kin = kinetic
            pot = .true.
            if( present( potential ) ) pot = potential
            pot_mgga = .true.
            if( present( potential_mgga ) ) pot_mgga = potential_mgga

            ! allocate integrals
            call me_mt_alloc( Smat_mt_basis )
            call me_mt_alloc( Hmat_mt_basis )

            allocate( rfun(1, nrmtmax) )

            do is = 1, nspecies
                do ia = 1, natoms(is)
                ias = idxas(ia, is)

                ! overlap
                rfun = 1.0_dp / y00
                call me_mt_prepare( is, ias, 0, zone, rfun, zzero, Smat_mt_basis(:, :, ias) )

                ! Hamiltonian
                ! kinetic energy
                if( kin ) then
                    rfun = 0.5_dp / y00
                    if( input%groundstate%ValenceRelativity /= 'none' ) &
                    rfun(1, 1:nrmt(is)) = rfun(1, 1:nrmt(is)) / (1.0_dp - 0.5_dp * alpha**2 * veffmt_gga(1, 1:nrmt(is), ias) * y00) 
                    call me_mt_prepare( is, ias, 0, zone, rfun, zone, Hmat_mt_basis(:, :, ias), gradient_product=.true. )
                end if
                ! potential
                if( pot ) then
                    call me_mt_prepare( is, ias, lmaxpot, zone, pot_mt(:, :, ias), zone, Hmat_mt_basis(:, :, ias) )
                end if
                !added non-multiplicative part. 
                if (pot_mgga) then 
                    call me_mt_prepare( is, ias, lmaxpot, zone, 0.5_dp * pot_mt_mgga(:, :, ias), zone, Hmat_mt_basis(:, :, ias), gradient_product=.true. )
                end if 

                end do
            end do
            deallocate( rfun )
        end subroutine gen_overlap_hamiltonian_mt_basis


        !> This subroutine sets up the Kohn-Sham Hamiltonian and the overlap matrix 
        !> for the given wavevector \({\bf k}\).
        subroutine gen_mGGA_H_and_S( ik, mGGA_H, mGGA_S )
            use constants, only: zzero
            !> index of the wavevector \({\bf k}\) in the set
            integer, intent(in) :: ik
            !> Hamilton matrix 
            complex(8), intent(inout) :: mGGA_H(:, :)
            !> Overlap matrix
            complex(8), intent(inout) :: mGGA_S(:, :)

            integer :: i, ip, ig
            integer :: nmatk, nmatmax, n
            integer, target :: ked_ngkmax
            integer :: is, ia, ias

            complex(dp), allocatable :: apwalm_k(:,:,:,:)
            integer, pointer :: ngkmax_ptr      
            complex(dp), allocatable :: zfft(:), dcfun_ir(:,:), dcfun_ig(:, :)
            complex(dp) :: grad(3,3)

            ! set matrix sizes
            nmatk = ked_Gkset%ngk(1, ik) + nlotot
            ked_ngkmax = ked_Gkset%ngkmax
            ngkmax_ptr => ked_ngkmax

            ! allocate local variables
            allocate( apwalm_k(ngkmax_ptr, apwordmax, ked_lmmaxapw, natmtot), source=zzero )
            CALL_ASSERT(size(mGGA_H,1) == nmatk .and. size(mGGA_H,2) == nmatk, message='Incorrect matrix size for `mGGA_H`.')
            CALL_ASSERT(size(mGGA_S,1) == nmatk .and. size(mGGA_S,2) == nmatk, message='Incorrect matrix size for `mGGA_S`.')
            mGGA_H = zzero
            mGGA_S = zzero

            ! * set up overlap and Hamiltonian matrix
            ! muffin-tin contribution
            call match(ked_Gkset%ngk(1, ik), ked_Gkset%gkc(:, 1, ik), ked_Gkset%tpgkc(:, :, 1, ik), ked_Gkset%sfacgk(:, :, 1, ik), apwalm_k )
            do is = 1, nspecies
                do ia = 1, natoms(is)
                ias = idxas(ia, is)
                call me_mt_mat( is, ias, ked_Gkset%ngk(1, ik), ked_Gkset%ngk(1, ik), apwalm_k(:, :, :, ias), apwalm_k(:, :, :, ias), zone, Smat_mt_basis(:, :, ias), zone, mGGA_S )
                call me_mt_mat( is, ias, ked_Gkset%ngk(1, ik), ked_Gkset%ngk(1, ik), apwalm_k(:, :, :, ias), apwalm_k(:, :, :, ias), zone, Hmat_mt_basis(:, :, ias), zone, mGGA_H )
                end do
            end do

            ! interstitial contribution
            call me_ir_mat( ked_Gkset, ik, ked_Gkset, ik, zone, cfun_ig, zone, mGGA_S, Gset_op=ked_Gset )
            call me_ir_mat( ked_Gkset, ik, ked_Gkset, ik, zone, pot_cfun_ig, zone, mGGA_H, Gset_op=ked_Gset )
            call me_ir_mat( ked_Gkset, ik, ked_Gkset, ik, zone, kin_cfun_ig, zone, mGGA_H, Gset_op=ked_Gset, gradient_product=.true. )

            ! added non-multiplicative part
            if ( allocated(pot_mgga_cfun_ig) ) call me_ir_mat( ked_Gkset, ik, ked_Gkset, ik, zone, pot_mgga_cfun_ig, zone, &
                                                              & mGGA_H, Gset_op=ked_Gset, gradient_product=.true. )

            ! deallocate local variables
            deallocate( apwalm_k )
        end subroutine gen_mGGA_H_and_S

        end module 
