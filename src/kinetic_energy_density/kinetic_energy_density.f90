!> Module provides routines to calculate the kinetic energy density (KED) \(\tau\) for both the interstitial region (`ked_ir`)
!> and the muffin-tin region (`ked_mt`). The kinetic energy density in the muffin-tin contains the KED from the core electrons (`ked_cr`) and the valence electrons 
!> from the muffin-tin region.
!> For the spin-polarised case, the spin kinetic energy density vectors (`ked_magmt` or `ked_magir`) are calculated as well. 
!> These vectors are defined analogously to the spin density (magnetisation) vectors (`magmt` or `magir`, see [[mod_potential_and_density(module)]]).
!>
!> Generally, the kinetic energy density for valence electrons is defined as:
!> \[
!>  \tau(\mathbf{r}) = \frac{1}{2} \sum_{n \mathbf{k}} w_{n \mathbf{k}} \nabla \Psi_{n \mathbf{k}}^{\dagger} \nabla \Psi_{n \mathbf{k}}
!> \]
!> and for the core electrons as 
!> \[
!>  \tau(\mathbf{r}) = \frac{1}{2} \sum_{k M} \nabla \Psi_{k M}^{\dagger} \nabla \Psi_{k M}
!> \]
!> A more detailed derivation of the calculation of the kinetic energy density can be found in the respective modules [[kinetic_energy_density_ir(module)]] for the interstitial region and 
!> [[kinetic_energy_density_mt(module)]] for the valence electrons in the muffin-tin region.
!> For the core elctrons a more detailed derivation can be found in: [[kinetic_energy_density_cr(module)]].
module kinetic_energy_density 
    use kinetic_energy_density_vars
    use kinetic_energy_density_mt, only: gen_ked_mt, gen_denmat_k
    use kinetic_energy_density_cr, only: gen_ked_cr
    use kinetic_energy_density_ir, only: gen_ked_ir
    use precision, only: dp

    use constants, only: zzero, zone
    use mod_atoms, only: nspecies, natmtot, spnrmax, spnr, natoms, idxas
    use modinput, only: input
    use mod_APW_LO, only: apwordmax, apwfr, lofr
    use modmpi, only : mpi_env_k, distribute_loop
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use mod_getoccsv, only: getoccsv
    use mod_spin, only: ncmag, nspnfv, ndmag
    use mod_muffin_tin, only : nrcmtmax, nrmt
    use muffin_tin_basis, only: mt_basis_type

    implicit none     
    private 

        !> Muffin-tin kinetic energy density
        real(dp), public, allocatable :: ked_mt(:, :, :)
        !> Muffin-tin spin kinetic energy density vector
        real(dp), public, allocatable :: ked_magmt(:, :, :, :)
        !> Core electrons kinetic energy density 
        real(dp), public, allocatable :: ked_cr(:, :, :)

        ! Interstitial basis
        !> Interstitial kinetic energy density
        real(dp), public,  allocatable :: ked_ir(:)
        !> Interstitial spin kinetic energy density vectors
        real(dp), public, allocatable :: ked_magir(:, :)


    public :: gen_ked, allocate_ked 
    
    contains 
    
    !>  Interface to calculate total kinetic energy density
    subroutine gen_ked
        if (associated(input%groundstate%spin)) then
            call gen_kinetic_energy_density_spin_polarised()
        else 
            call gen_kinetic_energy_density_spin_unpolarised()
        end if 
    end subroutine 

    subroutine allocate_ked
        integer :: ked_nrmtmax
        ked_nrmtmax = maxval(ked_mt_basis%n_rad_grid)

        ! Allocation of arrays
        if( allocated( ked_mt ) ) deallocate( ked_mt )
        allocate(ked_mt(ked_lmmaxvr, ked_nrmtmax, natmtot), source= 0.0_dp)
        
        if( allocated( ked_cr ) ) deallocate( ked_cr )
        allocate(ked_cr(ked_lmmaxvr, spnrmax, natmtot), source= 0.0_dp)
        
        if( allocated( ked_ir ) ) deallocate( ked_ir )
        allocate(ked_ir(ked_Gset%ngrtot), source= 0.0_dp)

        if (associated(input%groundstate%spin)) then
            if( allocated( ked_magmt ) ) deallocate( ked_magmt )
            allocate(ked_magmt(ked_lmmaxvr, ked_nrmtmax, natmtot, 3), source = 0.0_dp)

            if( allocated( ked_magir ) ) deallocate( ked_magir )
            allocate(ked_magir(ked_Gset%ngrtot, ndmag), source = 0.0_dp)
        end if 
    end subroutine 
    
    !> Interface to calculate total spin unpolarised kinetic energy density
    subroutine gen_kinetic_energy_density_spin_unpolarised()
        use mod_APW_LO, only: apword, nlorb, lorbl, apwfr, lofr

        ! First and second variational eigenvectors
        complex(dp), allocatable :: evecfv(:, :, :), evecsv(:, :)
        ! Matching coefficients
        complex(dp), allocatable :: apwalmk(:, :, :, :)
        ! Density matrics for the kinetic energy density and umt
        complex(dp), allocatable :: ked_mat_alpha(:, :, :)
        
        real(dp), allocatable :: ked_ir_k(:)
        real(dp), allocatable :: occsvk(:)
        integer, target :: ngkmax
        integer :: ik, firstk, lastk, ked_nrmtmax 
        integer, pointer :: ngkmax_ptr      
        real(dp) :: ts0, ts1, ts1_mat, ts0_mat, ts0_kedmt, ts1_kedmt, ts1_kedcr, ts0_kedcr
        real(dp) :: ts0_kedir, ts1_kedir
        integer:: is, ia, ias, ir , i, lm

        ! Update MT basis 
        ked_mt_basis%apw_rad_fun = apwfr     
        ked_mt_basis%lo_rad_fun = lofr

        ked_nrmtmax = maxval(ked_mt_basis%n_rad_grid)

        ! Allocation of arrays
        call allocate_ked()

        if (allocated( ked_mat_alpha )) deallocate(ked_mat_alpha)
        allocate(ked_mat_alpha(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max, natmtot), source=zzero)

        allocate(evecfv(nmatmax, nstfv, nspnfv))
        allocate(evecsv(nstsv, nstsv))
        
        if (allocated(occsvk)) deallocate(occsvk)
        allocate(occsvk(nstsv), source= 0.0_dp)

        allocate(ked_ir_k(ked_Gset%ngrtot), source = 0.0_dp)
        
        ngkmax = ked_Gkset%ngkmax
        ngkmax_ptr => ngkmax            

        if (allocated(apwalmk)) deallocate(apwalmk)
        allocate(apwalmk(ngkmax_ptr, apwordmax, ked_lmmaxapw, natmtot), source=zzero)

        call distribute_loop( mpi_env_k, ked_kset%nkpt, firstk, lastk )
        do ik = firstk, lastk
            ! get the eigenvectors from file
            call Getevecfv(ked_kset%vkl(:, ik), ked_Gkset%vgkl(:, :, :, ik), evecfv)
            call Getevecsv(ked_kset%vkl(:, ik), evecsv)

            call match(ked_Gkset%ngk(1, ik), ked_Gkset%gkc(:, 1, ik), ked_Gkset%tpgkc(:, :, 1, ik), ked_Gkset%sfacgk(:, :, 1, ik), apwalmk)

            call getoccsv( ked_kset%vkl(:,ik), occsvk(:))
            call gen_denmat_k(ik, occsvk(:), evecfv(:, :, 1), evecsv, apwalmk, ked_mat_alpha)    
            
            ! call genked_ir
            call gen_ked_ir(ik, evecfv(:, :, 1), ked_ir_k)
            ked_ir(:) = ked_ir(:) + ked_ir_k(:)
        end do
    
        call gen_ked_mt(ked_mt, ked_mat_alpha) 
        deallocate(ked_ir_k)

        call mpisum_ked_and_ked_mag(mpi_env_k)
        
        ! generate ked_cr
        call gen_ked_cr(ked_cr)

        ! add rhocr to rhomt
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia, is)
                do ir = 1, nrmt(is)
                    do lm = 1, ked_lmmaxvr
                        ked_mt(lm, ir, ias) = ked_mt(lm, ir, ias) + ked_cr(lm,ir, ias)
                    end do 
                end do 
            end do 
        end do

    end subroutine 

    !> Interface to calculate total spin-unpolarised kinetic energy density
    subroutine gen_kinetic_energy_density_spin_polarised()
        ! First and second variational eigenvectors
        complex(dp), allocatable :: evecfv(:, :, :), evecsv(:, :)
        ! Matching coefficients
        complex(dp), allocatable :: apwalmk(:, :, :, :)
        ! Density matrics for the kinetic energy density and umt
        complex(dp), allocatable :: ked_mat_alpha(:, :, :), ked_mat_beta(:, :, :), ked_mat_ab(:, :, :), & 
                                    ked_mat_alpha_plus_beta(:, :, :), ked_mat_alpha_min_beta(:, :, :)
        
        real(dp), allocatable :: ked_ir_k(:), ked_magir_k(:, :)
        real(dp), allocatable :: occsvk(:)
        integer :: is, ia, ias, ir , lm
        integer, target :: ngkmax
        integer :: ik, firstk, lastk, ked_nrmtmax 
        integer, pointer :: ngkmax_ptr      

        ! Update MT basis 
        ked_mt_basis%apw_rad_fun = apwfr     
        ked_mt_basis%lo_rad_fun = lofr

        ked_nrmtmax = maxval(ked_mt_basis%n_rad_grid)
        
        ! Allocation of arrays
        call allocate_ked()
        
        if (allocated( ked_mat_alpha )) deallocate(ked_mat_alpha)
        allocate(ked_mat_alpha(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max, natmtot), source=zzero)

        if( allocated( ked_mat_beta ) ) deallocate( ked_mat_beta )
        allocate(ked_mat_beta(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max, natmtot), source=zzero)

        if (ncmag) then 
            if( allocated( ked_mat_ab ) ) deallocate( ked_mat_ab )
            allocate(ked_mat_ab(ked_mt_basis%n_basis_fun_max, ked_mt_basis%n_basis_fun_max, natmtot), source=zzero)
        end if 

        allocate(evecfv(nmatmax, nstfv, nspnfv))
        allocate(evecsv(nstsv, nstsv))

        allocate(ked_ir_k(ked_Gset%ngrtot), source = 0.0_dp)
        allocate(ked_magir_k(ked_Gset%ngrtot, ndmag), source = 0.0_dp)
        
        if (allocated(occsvk)) deallocate(occsvk)
        allocate(occsvk(nstsv), source= 0.0_dp)

        ngkmax = ked_Gkset%ngkmax
        ngkmax_ptr => ngkmax
        
        if (allocated(apwalmk)) deallocate(apwalmk)
        allocate(apwalmk(ngkmax_ptr, apwordmax, ked_lmmaxapw, natmtot), source=zzero)

        call distribute_loop( mpi_env_k, ked_kset%nkpt, firstk, lastk )
        do ik = firstk, lastk
            ! get the eigenvectors from file
            call Getevecfv(ked_kset%vkl(:, ik), ked_Gkset%vgkl(:, :, :, ik), evecfv)
            call Getevecsv(ked_kset%vkl(:, ik), evecsv)

            ! generate kinetic energy density matrix
            call match(ked_Gkset%ngk(1, ik), ked_Gkset%gkc(:, 1, ik), ked_Gkset%tpgkc(:, :, 1, ik), ked_Gkset%sfacgk(:, :, 1, ik), apwalmk)

            call getoccsv( ked_kset%vkl(:,ik), occsvk(:))

            if (ncmag) then
                call gen_denmat_k(ik, occsvk(:), evecfv(:, :, 1), evecsv, apwalmk, ked_mat_alpha, ked_mat_beta, ked_mat_ab )    
            else 
                call gen_denmat_k(ik, occsvk(:), evecfv(:, :, 1), evecsv, apwalmk, ked_mat_alpha, ked_mat_beta )    
            end if 

            ! call genked_ir
            call gen_ked_ir(ik, evecfv(:, :, 1), evecsv, ked_ir_k, ked_magir_k)
            ked_ir(:) = ked_ir(:) + ked_ir_k(:)
            ked_magir(:, :) = ked_magir(:, :) + ked_magir_k(:, :)
        end do

        if (ncmag) then
            call gen_ked_mt(ked_mt, ked_magmt, ked_mat_alpha, ked_mat_beta, ked_mat_ab ) 
        else
            call gen_ked_mt(ked_mt, ked_magmt, ked_mat_alpha, ked_mat_beta ) 
        end if 

        call mpisum_ked_and_ked_mag(mpi_env_k)

        deallocate(ked_ir_k, ked_magir_k)

        ! generate ked_cr
        call gen_ked_cr(ked_cr)
        
        !add rhocr to rhomt
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia, is)
                ked_mt(1:ked_lmmaxvr, 1:nrmt(is), ias) = ked_mt(1:ked_lmmaxvr, 1:nrmt(is), ias) + ked_cr(1:ked_lmmaxvr, 1:nrmt(is), ias)
            end do 
        end do 
    end subroutine 

    subroutine mpisum_ked_and_ked_mag(mpi_env)
        use modinput, only: input
        use modmpi, only: mpiinfo, ierr
        use mod_spin, only: ndmag
        use mod_atoms, only: natmtot
#ifdef MPI
        use modmpi, only: MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
#endif
        type(mpiinfo), intent(in) :: mpi_env

        integer :: ked_nrmtmax        
        ked_nrmtmax = maxval(ked_mt_basis%n_rad_grid)

        ! quick exit
        if ( mpi_env%procs == 1 ) return
#ifdef MPI
        ! Muffin-tin kinetic energy density
        call MPI_barrier(mpi_env%comm, ierr)
        call MPI_allreduce(mpi_in_place, ked_mt, ked_lmmaxvr*ked_nrmtmax*natmtot, &
                            & MPI_DOUBLE_PRECISION, MPI_SUM, mpi_env%comm, ierr)

        ! Interstitial kinetic energy density
        call MPI_barrier(mpi_env%comm, ierr)
        call MPI_allreduce(mpi_in_place, ked_ir, ked_Gset%ngrtot, MPI_DOUBLE_PRECISION, &
                            & MPI_SUM, mpi_env%comm, ierr)
        
        if (associated(input%groundstate%spin)) Then
            ! Muffin-tin spin kinetic energy density
            call MPI_barrier(mpi_env%comm, ierr)
            call MPI_allreduce(mpi_in_place, ked_magmt, &
                                & ked_lmmaxvr*ked_nrmtmax*natmtot*ndmag, MPI_DOUBLE_PRECISION, MPI_SUM, &
                                & mpi_env%comm, ierr)

            ! Insterstitial spin kinetic energy density
            call MPI_barrier(mpi_env%comm, ierr)
            call MPI_allreduce(mpi_in_place, ked_magir, ked_Gset%ngrtot*ndmag, &
                                & MPI_DOUBLE_PRECISION, MPI_SUM, mpi_env%comm, ierr)
        end if
#endif
    end subroutine 

end module
