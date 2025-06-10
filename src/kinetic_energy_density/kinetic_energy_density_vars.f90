module kinetic_energy_density_vars
    use precision, only: dp
    use modinput
    use mod_kpointset, only: k_set, G_set, Gk_set
    use muffin_tin_basis, only: mt_basis_type

    implicit none
    private

    !> maximum \(l\) for basis functions
    integer, public :: ked_lmaxapw
    !> maximum number of \(l,m\) pairs in basis functions
    integer, public :: ked_lmmaxapw
    !> maximum \(l\) for kinetic energy density vector
    integer, public :: ked_lmaxvr
    !> maximum number of \((l,m)\) pairs in kinetic energy expansion
    integer, public :: ked_lmmaxvr
    !> muffin-tin basis object
    type(mt_basis_type), public :: ked_mt_basis
    !> number of non-zero Clebsch-Gordan coefficients for a given \((l,m)\) pair
    integer, allocatable, public :: ked_cg_num(:, :, :)
    !> \((l',m')\) pair of non-zero Clebsch-Gordan coefficients for a given \((l,m)\) pair
    integer, allocatable, public :: ked_cg_lm(:, :, :, :)
    !> value non-zero Clebsch-Gordan coefficients for a given \((l,m)\) pair
    complex(dp), allocatable, public :: ked_cg_val(:, :, :, :)
    !> Timing for total kinetic energy density
    real(dp), public :: timeked

    !> set of electronic \({\bf k}\) vectors
    type(k_set), public :: ked_kset

    !> set of \({\bf G}\) vectors for the expansion of Fourier series
    !> with cutoff `gmaxvr`
    type(G_set), public :: ked_Gset
    !!> set of \({\bf G}\) vectors with twice the cutoff `2*gmaxvr`
    type(G_set), public :: ked_2Gset
    !> set of \({\bf G+k}\) vectors for the (L)APW expansion
    type(Gk_set), public :: ked_Gkset

    public :: ked_var_init, ked_var_free, & 
            gen_kpoints_bandstructure

    contains

    !> Initializes all global variables that remain constant during the calculation
    !> needed to calculate the kinetic energy density.
    subroutine ked_var_init
        use muffin_tin_basis, only: mt_basis_type, generate_non_zero_clebsch_gordan
        use mod_atoms, only: spr, nspecies
        use mod_muffin_tin, only: nrmt
        use modinput, only: input
        use mod_APW_LO, only: apword, nlorb, lorbl, apwfr, lofr
        use mod_lattice, only: bvec
        use mod_Gkvector, only: gkmax
        use mod_kpointset, only: generate_k_vectors, generate_G_vectors, generate_Gk_vectors
        use gaunt
        use mod_misc, only: task

        call initialise_mgga_timings()

        ked_lmaxapw = input%groundstate%lmaxapw
        ked_lmmaxapw = (ked_lmaxapw+1)**2
        ked_lmaxvr = input%groundstate%lmaxvr
        ked_lmmaxvr = (ked_lmaxvr + 1)**2

        ked_mt_basis = mt_basis_type(spr(:, 1:nspecies), nrmt(1:nspecies), apwfr, lofr, ked_lmaxapw, &
                                     apword(:, 1:nspecies), nlorb(1:nspecies), lorbl(:, 1:nspecies))
                                     
        ! generate Clebsch-Gordan coefficients
        call generate_non_zero_clebsch_gordan( ked_lmaxapw, 1e-64_dp, ked_cg_num, ked_cg_lm, ked_cg_val )

        ! check if Gaunt coefficients are available and create them if not
        if(.not. gaunt_coeff_yry%check_bounds(ked_lmaxapw, ked_lmaxvr, ked_lmaxapw) ) then 
            gaunt_coeff_yry = non_zero_gaunt_yry( ked_lmaxapw+1, ked_lmaxvr+1, ked_lmaxapw+1 ) 
        end if

        ! generate G-vectors for plane wave expansion
        call generate_G_vectors( ked_Gset,  bvec, [[0,0,0],[0,0,0]],   input%groundstate%gmaxvr, auto_intgv=.true. )
        !call generate_G_vectors( ked_2Gset, bvec, [[0,0,0],[0,0,0]], 2 * input%groundstate%gmaxvr, auto_intgv=.true. )


        if (task == 20) then 
            ! for band structure plots generate k-points along a line
            call gen_kpoints_bandstructure
        else
            ! generate k-point set
            call generate_k_vectors( ked_kset, bvec, &
                    input%groundstate%ngridk, &
                    input%groundstate%vkloff, &
                    input%groundstate%reducek, &
                    uselibzint=.false. )
        end if 

        ! generate G+k vectors
        call generate_Gk_vectors( ked_Gkset, ked_kset, ked_Gset, gkmax )
    end subroutine 

    subroutine gen_kpoints_bandstructure
        use modinput, only: input
        use mod_lattice, only: bvec
        use mod_kpointset, only: generate_k_vectors

        integer :: iq
        integer :: ked_nvp1d, ked_npp1d
        real(dp), allocatable :: ked_dvp1d(:), ked_vplp1d(:, :), ked_dpp1d(:)
        
        ked_nvp1d = size(input%properties%bandstructure%plot1d%path%pointarray)
        ked_npp1d = input%properties%bandstructure%plot1d%path%steps
        
        allocate(ked_dvp1d(ked_nvp1d))
        allocate(ked_vplp1d(3, ked_npp1d))
        allocate(ked_dpp1d(ked_npp1d))

        call connect(bvec, ked_nvp1d, ked_npp1d, input%properties%bandstructure%plot1d%path%pointarray, ked_vplp1d, ked_dvp1d, ked_dpp1d)
        call generate_k_vectors(ked_kset, bvec, (/1, 1, ked_npp1d/), (/0.d0, 0.d0, 0.d0/), .false.)
        
        do iq = 1, ked_kset%nkpt
            ked_kset%vkl( :, iq) = ked_vplp1d( :, iq)
            call r3mv(bvec, ked_kset%vkl( :, iq), ked_kset%vkc( :, iq))
        end do
    end subroutine 

    !> Initialises the timings for mgga calculations.
    subroutine initialise_mgga_timings()
        use mod_timing
        timeked = 0.0_dp
    end subroutine

    subroutine ked_var_free
        use mod_kpointset, only: delete_k_vectors, delete_G_vectors, delete_Gk_vectors
        call delete_k_vectors( ked_kset )
        call delete_G_vectors( ked_Gset )
        call delete_Gk_vectors( ked_Gkset )
        call ked_mt_basis%destroy
    end subroutine 

end module kinetic_energy_density_vars
