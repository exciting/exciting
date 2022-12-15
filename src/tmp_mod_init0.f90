!> Short-lived module to house refactored initialisation routines for init0.f90
module tmp_mod_init0
    use precision, only: dp
    use modinput, only: input_type, groundstate_type, isspinorb
    use modmpi, only: terminate_if_false, terminate

    implicit none
    private
    public :: initialise_xc_mixing_coefficients, allocate_coulomb_potentials, &
              initialise_groundstate_timings, map_atoms_per_species_to_atomic_index

contains

    !> Map atoms per species to an atomic index over all atoms in the system
   subroutine map_atoms_per_species_to_atomic_index (nspecies, natoms, idxas, natmmax, natmtot)
      use constants, only: maxatoms, maxspecies

      !> Number of species
      integer, intent(in) :: nspecies
      !> Number of atoms for each species
      integer, intent(in) :: natoms (:)
      !> Maximum number of atoms over all the species
      integer, intent(out) :: natmmax
      !> Total number of atoms
      integer, intent(out):: natmtot
      !> Map atoms per species to an atomic index over all atoms in the system
      integer, intent(out) :: idxas (maxatoms, maxspecies)

      !> Local variables
      integer :: ias, is, ia

      natmmax = 0
      ias = 0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = ias + 1
            idxas (ia, is) = ias
         End Do
         natmmax = Max (natmmax, natoms(is))
      End Do
      natmtot = ias

    end subroutine


    !> Initialise ground state timings.
    !>
    !> This deliberately uses globals because
    !> it does not make sense to pass all as args OR even wrap
    !> into a type. Would be better-implemented as a dict.
    subroutine initialise_groundstate_timings()
        use mod_timing
        timeinit = 0.d0
        timemat = 0.d0
        timefv = 0.d0
        timesv = 0.d0
        timerho = 0.d0
        timepot = 0.d0
        timefor = 0.d0
        timeio=0d0
        timemt=0d0
        timemixer=0d0
        timematch=0d0
        time_hmlaan=0d0
        time_hmlalon=0d0
        time_hmllolon=0d0
        time_olpaan=0d0
        time_olpalon=0d0
        time_olplolon=0d0
        time_hmlistln=0d0
        time_olpistln=0d0
        time_rdirac=0d0
        time_rschrod=0d0
        time_oepvnl=0.0d0
        time_oep_iter=0.0d0
    end subroutine


    !> Initialise XC mixing coefficients.
    !>
    !> Check mod_potential_and_density for declarations
    subroutine initialise_xc_mixing_coefficients(gs_input, xctype, xcdescr, xcspin, &
                                                           xcgrad, ex_coef, ec_coef)
        use modxcifc, only: getxcdata
        use vx_enums, only: HYB_PBE0, HYB_HSE
        !> Ground state input XML object
        type(groundstate_type), intent(in) :: gs_input
        !> exchange-correlation functional type
        integer, intent(out) :: xctype(3)
        !> exchange-correlation functional description
        character (len=512) , intent(out) :: xcdescr
        !> exchange-correlation functional spin treatment
        integer, intent(out) :: xcspin
        !> exchange-correlation functional density gradient treatment
        integer, intent(out) :: xcgrad
        !> exchange mixing parameter for hybrid functionals
        real(dp), intent(out) :: ex_coef
        !> correlation mixing parameter for hybrid functionals
        real(dp), intent(out) :: ec_coef

        if  (associated(gs_input%HartreeFock) .and. associated(gs_input%OEP)) then
           write (*,*)
           write (*, '("Error(init0): illegal choice for exact exchange")')
           write (*, '("You cannot use HF and OEP simultaneously")')
           write (*,*)
           call terminate()
        endif

        call getxcdata(xctype, xcdescr, xcspin, xcgrad, ex_coef)

        if (isspinorb()) then
            if (xctype(1)==HYB_PBE0 .or. xctype(1)==HYB_HSE) then
                call terminate_if_false(gs_input%spin%realspace,'("Error(init0): &
                  & `input%groundstate%spin%realspace` needs to be set to true, for &
                  & calculations with hybrid functionals which account for spin-orbit coupling effects.")')
            end if
        endif

        ! Default - No mixing
        ex_coef = 0.0_dp
        ec_coef = 1.0_dp

        if (gs_input%xctypenumber < 0) ex_coef = 1.0_dp

        ! Mixing parameters when using hybrids
        if  (xctype(1) == HYB_PBE0 .or. xctype(1) == HYB_HSE) then
            ex_coef = gs_input%Hybrid%excoeff
            ec_coef = gs_input%Hybrid%eccoeff
        endif

    end subroutine


    !> Initialise Coulomb potentials: INT, MT and Madeulung.
    !>
    !> Integer sizes are declared in mod_muffin_tin.F90 and .mod_Gvector.F90
    subroutine allocate_coulomb_potentials(lmmaxvr, nrmtmax, natmtot, ngrtot, &
        vclmt, vclir, vmad)

        !> Maximum angular momentum for potentials and densities, squared.
        integer, intent(in) :: lmmaxvr
        !> Max number of muffin-tin radial points, for all species.
        integer, intent(in) :: nrmtmax
        !> Total number of atoms in system.
        integer, intent(in) :: natmtot
        !> Total number of G vectors.
        integer, intent(in) :: ngrtot

        !> Muffin-tin Coulomb potential
        real(dp), allocatable, intent(inout) :: vclmt(:, :, :)
        !> Interstitial real-space Coulomb potential
        real(dp), allocatable, intent(inout) :: vclir(:)
        !> Madeulung potential (excluding the on-site nuclear contribution)
        !> in the innermost radial point
        real(dp), allocatable, intent(inout) :: vmad(:)

        if (allocated(vclmt)) deallocate(vclmt)
        allocate (vclmt(lmmaxvr, nrmtmax, natmtot))

        if (allocated(vclir)) deallocate(vclir)
        allocate (vclir(ngrtot))

        if (allocated(vmad)) deallocate(vmad)
        allocate (vmad(natmtot))
    end subroutine

end module 
