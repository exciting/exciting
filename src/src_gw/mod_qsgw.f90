!> Module for handling and manipulating quasiparticle self-consistent
!> GW (QSGW) quantities, such as optimized XC potentials, and file
!> preparation between iterations.
module mod_qsgw

    use precision,                  only: i32, dp, str_512
    use modmpi,                     only: terminate_if_false, mpiglobal, barrier, terminate
    use gw_io,                      only: read_from_file, write_to_file, open_file
    use mod_band_to_lapw_transform, only: overlap_rootname
    use mod_offdiagonal_selfenergy, only: read_offdiagonal_selfenergy_correlation, read_offdiagonal_selfenergy_exchange, &
                                          delete_offdiagonal_selfenergy
#include "offload.fpp"

    implicit none

    private
    public  :: prepare_next_iteration_qsgw, prepare_current_iteration_qsgw, compute_optimized_vxc, &
               write_optimized_vxc_to_a_file, read_optimized_vxc_to_a_file, check_convergence_qsgw


    !> Convergence limit for the QSGW cycle
    real(dp), public, protected :: qsgw_eps

    !> Whether to compute the off-diagonal self-energy
    logical, public, protected :: offdiagonal_selfenergy = .false.

    !> Whether to add the optimized potential
    logical, public, protected :: add_optimized_vxc      = .false.

    !> Rootname for the optimized potential files
    character(len=*), parameter :: optimized_vxc_rootname = "OPTIMIZED_VXC_K"

    !> File suffix for the initial ground-state run
    character(len=*), public, parameter :: initial_gs_suffix = '_INITIAL_GS.OUT'

    !> File suffix for the previous ground-state run
    character(len=*), private, parameter :: previous_gs_suffix = '_PREVIOUS_GS.OUT'

    !> File suffix for the standard ground-state run
    character(len=*), private, parameter :: standard_gs_suffix = '.OUT'

    !> Reference densities of the previous run in the MT
    real(dp), allocatable :: rho_mt_ref(:,:,:)
    !> Reference densities of the previous run in the MT
    real(dp), allocatable :: rho_ir_ref(:)
    

contains

    !> Save/move/link QSGW quantities for the **next iteration**.
    !> Prepares files and optimized potentials needed in the
    !> `next_iteration` directory.
    subroutine prepare_next_iteration_qsgw(iteration, file_format)

        use mod_atoms, only: nspecies
        use modinput,  only: input

        !> Current iteration step
        integer(i32), intent(in) :: iteration
        !> File format for large I/O files
        character(len=*), intent(in) :: file_format

        integer(i32) :: ierror, ispecies, ipos
        character(str_512) :: buffer
        character(str_512) :: speciesfile, scffile
        logical :: file_exists

        character(len=*), parameter :: optimized_command = &
            'for f in OPTIMIZED_VXC_K*.OUT; do ' // &
            'cp "$f" "next_iteration/CURRENT_${f}"; done'

        if (mpiglobal%rank == 0) then

            ! Copy the optimized exchange-correlation potential
            call execute_command_line(trim(optimized_command), exitstat=ierror)
            call terminate_if_false(ierror == 0, message="Error(prepare_next_iteration_qsgw): copying optimized potential failed")

            ! Copy input.xml (or argument if provided)
            buffer = "input.xml"
            if (command_argument_count() == 1) call getarg(1, buffer)
            call copy_files(".", buffer, "", "next_iteration", buffer, "")

            ! Copy species files
            do ispecies = 1, nspecies
                call copy_files(input%structure%speciespath, &
                                 input%structure%speciesarray(ispecies)%species%speciesfile, &
                                 "", "next_iteration", "", "")
            end do

            ! Copy STATE.OUT from this run for the next one
            call copy_files(".", "STATE", standard_gs_suffix, "next_iteration", "STATE", previous_gs_suffix)

            ! Copy initial state information
            ! we need this one to recompute the radial
            ! function in each step
            if (iteration == 0) then
                if (associated(input%groundstate%hybrid)) then
                    ! For hybrid runs copy the PBE information
                    ! as it is the one needed to regenerate the basis set
                    call copy_files(".", "STATE_PBE", standard_gs_suffix, "next_iteration", "STATE", initial_gs_suffix)
                else
                    call copy_files(".", "STATE", standard_gs_suffix, "next_iteration", "STATE", initial_gs_suffix)
                end if
            else
                call copy_files(".", "STATE", initial_gs_suffix, "next_iteration", "STATE", initial_gs_suffix)
            end if

        end if

        call barrier(mpiglobal)

    end subroutine prepare_next_iteration_qsgw


    !> Generates all the GS information needed for 
    !> the QSGW
    subroutine prepare_current_iteration_qsgw(iteration)

        use modinput,  only: input
        use mod_misc,  only: task, filext, versionname
        use mod_potential_and_density, only: xctype
        use modxs, only: isreadstate0
        use mod_potential_and_density, only: rhoir, rhomt
        use mod_gen_lo, only: genlofr
        use mod_timing, only: timeinit, timemat, timefv, timesv, timerho, timepot, timematch, &
                              timeio, timemt, timemixer, time_rdirac, time_rschrod

        !> Current iteration step
        integer(i32), intent(in) :: iteration

        integer(i32) :: ierror
        character(str_512) :: buffer
        logical :: restart_gs
        integer(i32), parameter :: infodotout_file_id = 60
        integer(i32), parameter :: totenergy_file_id = 61
        integer(i32), parameter :: rmsdveff_file_id = 65

        ! Consistency checks
        if (associated(input%groundstate%spin)) &
            call terminate("QSGW is not compatible with spin-polarized calculations")
        if (input%groundstate%tforce) then
            call warning("QSGW is not compatible with the computation of forces. " // &
                         "Setting tforce to .false.")
            input%groundstate%tforce = .false.
        end if
        if (associated(input%groundstate%oep)) &
            call terminate("QSGW is not compatible with OEP calculations")

        ! Set some values
        offdiagonal_selfenergy = .true.
        add_optimized_vxc      = .true.

        ! Importantly, current implementation uses exclusively
        ! "Extended linear tetrahedron method for the calculation of q-dependent
        !  dynamical response functions", to be published in Comp. Phys. Commun. (2010)
        input%groundstate%stypenumber = -1

        ! Set the empty states to the proper value
        input%groundstate%nempty = input%gw%nempty  

        ! Beyond the first iteration:
        ! prepare/modify globals so that the XC is replaced
        ! by the optimized potential.
        if (iteration /= 0) then

            ! Set the species path to the current directory
            input%structure%speciespath = "."           

            ! Allow for reading files with different extensions
            isreadstate0 = .false.

            ! Init global variables
            call init0
            call init1

            ! Prepare initial KS GS information
            ! Same logic than in hybrids when starting from the
            ! PBE functional.
            ! Note that this is always initialized
            filext = initial_gs_suffix
            call readstate()
            call gencore()          ! generate the core wavefunctions and densities
            call linengy()          ! find the new linearization energies
            call genapwfr()         ! generate the APW radial functions
            call genlofr()          ! generate the local-orbital radial functions
            call olprad()           ! compute the overlap radial integrals
            call energykncr()       ! core kinetic energy
            filext = standard_gs_suffix

            ! Reset references to hybrids or OEP
            ! That is from step 1 onwads all will 
            ! be a normal SCF but with an external
            ! fixed xc correlation potential
            if (associated(input%groundstate%hybrid)) nullify(input%groundstate%hybrid)
            if (associated(input%groundstate%oep))    nullify(input%groundstate%oep)

            ! Check if we are restarting the GS calculation
            restart_gs = input%groundstate%do == "fromfile"

            ! For QSGW beyond the first iteration `fromscratch`
            ! requires actually reading the previous iteration 
            ! information so that if not set to `skip`
            ! the code internally treats `fromscratch` as a
            ! `fromfile`.
            if (input%groundstate%do /= "skip") then
                ! For actual restarts we already have a STATE file so we do not to
                ! obtain it from the last run.
                if (.not. restart_gs .and. mpiglobal%rank == 0) then
                    call copy_files(".", "STATE", previous_gs_suffix, ".", "STATE", standard_gs_suffix)
                end if
                input%groundstate%do = "fromfile"

                ! Mix the potential. From Nora's thesis this seems more stable than
                ! density mixing, specially for oxides in which the latter generates
                ! too many oscillations
                input%groundstate%mixerswitch = 1

                ! Energy is no longer meaningfull quatity for the self-consistency
                ! we use the density to check the convergence.
                input%groundstate%scfconv = 'charge' 

                ! XC from input file is ignored. To detail the Vxc of normal DFT
                ! is replaced by an optimized potential
                xctype(1) = 1
                input%groundstate%xctypenumber = 1
                xctype(2) = 0
                xctype(3) = 0

                call barrier(mpiglobal)

                ! Read the previous step density
                call readstate()

                ! Special handling for the second iteration
                ! Note that for that iteration the previous STATE.OUT
                ! contains XC and Hartree, since we only want the latter
                ! we recompute the potential, which sets the vxc potential
                ! to 0.
                if (iteration == 1 .and. .not. restart_gs) then
                    call poteff(.true.)
                end if

                ! Restore the default
                isreadstate0 = .true.

                ! Get the previous densities 
                allocate(rho_mt_ref, source=rhomt)
                allocate(rho_ir_ref, source=rhoir)

                if (mpiglobal%rank == 0) then
                    open(infodotout_file_id, File='INFO'//trim(filext), Action='WRITE', Form='FORMATTED')
                    open(totenergy_file_id, File='TOTENERGY'//trim(filext), Action='WRITE', Form='FORMATTED')
                    open(rmsdveff_file_id, File='RMSDVEFF'//trim(filext), Action='WRITE', Form='FORMATTED')
                    call writeinfo(infodotout_file_id)
                endif 

                ! Executing the self-consistent cycle
                call scf_cycle(input%groundstate%outputlevelnumber)
                
                ! Printing out information
                if  (mpiglobal%rank == 0) then
                    if (input%groundstate%outputlevelnumber > 1) then
                        buffer = ''
                        write(buffer,'("Timings (seconds)")') 
                        call printbox(infodotout_file_id,"-",buffer)
                        Write (infodotout_file_id, '(" Initialisation", T45, ": ", F12.2)') timeinit
                        Write (infodotout_file_id, '(" Hamiltonian and overlap matrix set up", T45, ": ", F12.2)') timemat
                        Write (infodotout_file_id, '(" First-variational secular equation", T45, ": ", F12.2)') timefv
                        If (associated(input%groundstate%spin)) Then
                            Write (infodotout_file_id, '(" Second-variational calculation", T45, ": ", F12.2)') timesv
                        End If
                        Write (infodotout_file_id, '(" Calculation of charge-density", T45, ": ", F12.2)') timerho
                        Write (infodotout_file_id, '(" Calculation of potential", T45, ": ", F12.2)') timepot
                        Write (infodotout_file_id, '(" Muffin-tin manipulations", T45, ": ", F12.2)') timemt
                        Write (infodotout_file_id, '(" APW matching", T45, ": ", F12.2)') timematch
                        Write (infodotout_file_id, '(" Disk reads/writes", T45, ": ", F12.2)') timeio
                        Write (infodotout_file_id, '(" Mixing efforts", T45, ": ", F12.2)') timemixer
                        Write (infodotout_file_id, '(" Solver of Dirac eqn.", T45, ": ", F12.2)') time_rdirac
                        Write (infodotout_file_id, '(" Solver of rel. Schroedinger eqn.", T45, ": ", F12.2)') time_rschrod
                        Write (infodotout_file_id, '(" Total time spent in radial solvers", T45, ": ", F12.2)') time_rdirac+time_rschrod
                    end if
                    call printline(infodotout_file_id, "=")
                    write(buffer, '("EXCITING ", a, " stopped")') trim(versionname)
                    call printtext(infodotout_file_id, "=", buffer)
                    call printline(infodotout_file_id, "=")

                    close(infodotout_file_id)
                    close(totenergy_file_id)
                    close(rmsdveff_file_id)
                end if

                ! Set the task to skip. This prevents
                ! further elements to contamine the 
                ! now correct globals.
                input%groundstate%do = "skip"

                ! Now remove some tasks if present
                ! as beyond the first iteration they make 
                ! nonsense
                if (associated(input%gw%taskGroup%QPEigenvalues)) nullify(input%gw%taskGroup%QPEigenvalues)
                if (associated(input%gw%taskGroup%vxc)) nullify(input%gw%taskGroup%vxc)

            end if

        end if

        call barrier(mpiglobal)

    end subroutine prepare_current_iteration_qsgw

    !> Copy a file using the system `cp` command (Unix/Linux only).
    !> Constructs a file name from root + suffix and copies to destination.
    subroutine copy_files(source_from, fileroot_from, suffix_from, source_to, fileroot_to, suffix_to)
        !> Source directory, root name, and suffix (e.g. ".OUT")
        character(len=*), intent(in) :: source_from, fileroot_from, suffix_from
        !> Destination directory, root name, and suffix (optional)
        character(len=*), optional, intent(in) :: source_to, fileroot_to, suffix_to

        integer(i32) :: ierror
        character(len=:), allocatable :: cmd

        if (present(source_to) .and. present(fileroot_to) .and. present(suffix_to)) then
            cmd = "cp " // trim(source_from) // "/" // trim(fileroot_from) // trim(suffix_from) &
                        // " " // trim(source_to)   // "/" // trim(fileroot_to)   // trim(suffix_to)
        else
            cmd = "cp " // trim(source_from) // "/" // trim(fileroot_from) // trim(suffix_from) // " ."
        end if

        call execute_command_line(cmd, exitstat=ierror)
        call terminate_if_false(ierror == 0, "Error(copy_files): Unable to copy file: " // trim(cmd))

    end subroutine copy_files

    !> Computes the optimized potential for QSGW from the selfenergy
    subroutine compute_optimized_vxc(ik, offdiagonal, vxc_opt, file_format)

        use modinput,                   only: input
        use mod_eigenvalue_occupancy,   only: nstfv
        use mod_bands,                  only: evalfv
        use mod_frequency,              only: frequency, generate_freqgrid, delete_freqgrid
        use mod_selfenergy,             only: freq_selfc, selfec, selfex
        use mod_offdiagonal_selfenergy, only: add_offdiag_selfenergy_at_energies_to_optimized_xc_potential

        !> The irreducible k-point index
        integer(i32), intent(in) :: ik
        !> The optimized potential
        complex(dp),  allocatable, intent(inout) :: vxc_opt(:,:,:)
        !> Use or not the offdiagonal terms of the self-energy to build the optimized potential
        logical, intent(in) :: offdiagonal
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format
      
        integer(i32) :: ie1, i, n
        complex(dp)  :: dummy, sc
        type(frequency) :: real_freq_selfc

        call generate_freqgrid(real_freq_selfc, &
                               input%gw%selfenergy%wgrid%type, &
                               'refreq', &
                               input%gw%selfenergy%wgrid%size, &
                               input%gw%selfenergy%wgrid%wmin, &
                               input%gw%selfenergy%wgrid%wmax)

        !!$omp parallel do default(none) &
        !!$omp shared(ik, nstfv, real_freq_selfc, selfec, selfex, evalfv, vxc_opt) &
        !!$omp private(ie1, sc, dummy)
        do ie1 = 1, nstfv

            call get_selfc(real_freq_selfc%nomeg, real_freq_selfc%freqs, &
              selfec(ie1,:,1), evalfv(ie1, ik), &
              sc, dummy)
      
            vxc_opt(ie1,ie1,ik) = real(sc + selfex(ie1, 1), kind=dp)

        end do
        !!$omp end parallel do

        call delete_freqgrid(real_freq_selfc)
        
        if (offdiagonal) then
           call read_offdiagonal_selfenergy_correlation([ik], file_format)
           call read_offdiagonal_selfenergy_exchange([ik], file_format)
           call add_offdiag_selfenergy_at_energies_to_optimized_xc_potential( &
                 [ik],  &
                 vxc_opt)
           call delete_offdiagonal_selfenergy()
        end if
      
    end subroutine compute_optimized_vxc
      
    subroutine write_optimized_vxc_to_a_file(vxc_opt, ik, file_format)
        use gw_io, only: build_file_name, write_to_file
        !> Optimized potential
        complex(dp),  allocatable, intent(in) :: vxc_opt(:,:)
        !> Index of the current k-point
        integer(i32), intent(in) :: ik
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format
        
        character(len=str_512) :: file_name
        integer(i32) :: lbounds(2)

        lbounds = lbound( vxc_opt )

        call build_file_name( optimized_vxc_rootname, ik, file_name )
        call write_to_file( file_name, vxc_opt(:,:), lbounds(1:2), file_format)

    end subroutine write_optimized_vxc_to_a_file

    subroutine read_optimized_vxc_to_a_file(vxc_opt, ik, file_format, add_current_prefix)
        use gw_io, only: build_file_name, read_from_file
        !> Optimized potential
        complex(dp),  allocatable, intent(out) :: vxc_opt(:,:)
        !> Index of the current k-point
        integer(i32), intent(in) :: ik
        !> Format of the file where to print. It can be e.g. 'text' or 'binary'
        character(len=*), intent(in) :: file_format
        !> Does it has the prefix current
        logical, intent(in) :: add_current_prefix
        
        character(len=*), parameter :: current_prefix = "CURRENT_"
        character(len=str_512) :: file_name
        character(len=str_512) :: rootname

        ! In case we want to add the current prefix to the root
        ! of the filename
        if (add_current_prefix) then
            rootname = current_prefix // optimized_vxc_rootname
        else 
            rootname = optimized_vxc_rootname
        end if
        
        call build_file_name( rootname, ik, file_name )
        call read_from_file( file_name, vxc_opt, file_format)

    end subroutine read_optimized_vxc_to_a_file

    subroutine check_convergence_qsgw(iteration, convergence_tolerance)
        use mod_charge_and_moment, only : chgdst
        use mod_misc, only: filext
        use modinput, only: input
        use modmpi, only: mpiglobal

        integer(i32), intent(in) :: iteration
        real(dp),     intent(in) :: convergence_tolerance
        
        logical      :: check_convergence
        integer(i32) :: iunit

        ! Guard against procedure calls when
        ! the GS reference densities have not
        ! been initialized (e.g., GS run skipped).
        if (.not. allocated(rho_mt_ref) .or. &
            .not. allocated(rho_ir_ref)) return

        if (iteration /= 0 ) then 
            call chgdist(rho_mt_ref , rho_ir_ref)
        else 
            chgdst = huge(1.0_dp)
        end if

        check_convergence = merge(.false., chgdst < convergence_tolerance, iteration == 0)

        if (mpiglobal%rank == 0) then
            open(newunit=iunit, file='QSGW_CONVERGENCE.OUT', &
                 action='WRITE', status='REPLACE')
            write(iunit, *) "GS run within a QSGW calculation : "

            if (check_convergence) then
                write(iunit, *) " Charge density is converged with respect to the last GS"
            else
                write(iunit, *) " Charge density is not converged with respect to the last GS"
            end if

            if (iteration /= 0) then
                write(iunit, *) " Charge distance : ", chgdst, &
                                " ( target : ", convergence_tolerance, ")"
            else 
                write(iunit, *) " First iteration "
            end if
            
            close(iunit)
        end if
    
        if (check_convergence) then
            if (associated(input%gw)) nullify(input%gw)
        end if

    end subroutine check_convergence_qsgw

end module mod_qsgw
