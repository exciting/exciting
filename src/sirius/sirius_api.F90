! Copyright (C) 2013-2022 exciting team and Anton Kozhevnikov 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!> Interfaces to sirius library routines.
!>
!> Development Notes
!> ---------------------
!> Due to the large number of globals in exciting, much of exciting's input data is passed via
!> use statements to this module's routines, however scope is restricted.
!>
!> Alex notes that this API is equivalent to an OO approach, but implemented
!> via a module. The handlers and pointers are private in this scope, and never need to be 
!> exposed to exciting. All interaction with the SIRIUS lib is done in vai routines defined
!> in this module.
!> 
!> As there are no dummy overloads (for when sirius is not present), every routine must include 
!> preprocessor statements around the routine body/implementation. 
!>
!> To improve upon this design:
!>  * One could migrate to passing vaiables to and from these API routines, rather than
!>    passing as globals. This was started in sirius_init.F90.
!>  * This would also involve breaking up routines into smaller units.
!>  * Context, ground state and k-point handlers could be encapulated in an object
!>    instance, rather than this module.
!>  * One could split this module up into a collection of smaller modules.
module sirius_api 
  use iso_c_binding, only: c_ptr, c_null_ptr 
#ifdef SIRIUS
  use sirius, only: sirius_context_handler,    &
                    sirius_ground_state_handler, &
                    sirius_kpoint_set_handler, &
                    sirius_initialize_context, &
                    sirius_create_context,     &
                    sirius_dump_runtime_setup, &
                    sirius_import_parameters,  &
                    sirius_set_parameters,     &
                    sirius_set_lattice_vectors, &
                    sirius_add_atom_type,      &
                    sirius_set_atom_type_radial_grid, &
                    sirius_set_atom_type_configuration,&
                    sirius_add_atom_type_aw_descriptor,&
                    sirius_add_atom_type_lo_descriptor,&
                    sirius_add_atom,&
                    sirius_set_equivalent_atoms,&
                    sirius_add_xc_functional,&
                    sirius_get_step_function,&
                    sirius_get_fv_eigen_values,&
                    sirius_get_fv_eigen_vectors,&
                    sirius_get_sv_eigen_vectors,&
                    sirius_get_band_energies,&
                    sirius_get_forces, &
                    sirius_find_eigen_states,&
                    sirius_set_periodic_function,&
                    sirius_set_periodic_function_ptr,&
                    sirius_get_periodic_function,&
                    sirius_generate_coulomb_potential,&
                    sirius_generate_xc_potential,&
                    sirius_set_radial_function,&
                    sirius_create_ground_state,&
                    sirius_initialize_kset,&
                    sirius_create_kset,&
                    sirius_get_kpoint_inner_comm,&
                    sirius_get_num_gvec,&
                    sirius_get_gvec_arrays,&
                    sirius_get_fft_index,&
                    sirius_print_timers,&
                    sirius_free_handler,&
                    sirius_dump_runtime_setup,&
                    sirius_set_mpi_grid_dims,&
                    sirius_get_kpoint_inter_comm,&
                    sirius_set_band_occupancies,&
                    sirius_generate_density,&
                    sirius_get_max_num_gkvec,&
                    sirius_get_gkvec_arrays,&
                    sirius_finalize
#endif
  use precision, only: dp
  use modinput, only: input_type
  use mod_atoms, only: natmtot, spnst, natoms, idxas, nspecies, spzn, spn, spsymb, spr, spl, spmass, spk, spocc, spcore, spname 
  use mod_spin, only: ndmag, nspnfv
  use mod_APW_LO, only: apwve, apword, lorbve, lorbord, nlorb, apwe0, apwdm, lorbl, lorbe0, lorbdm, maxlorb, maxlorbord
  use mod_symmetry, only: eqatoms
  use mod_charge_and_moment, only: chgval
  use mod_eigenvalue_occupancy, only: nstfv, nstsv, evalsv, occsv, get_nstfv
  use mod_eigensystem, only: nmatmax
  use mod_muffin_tin, only: nrmt
  use mod_kpoint, only: nkpt, vkl, wkpt
  use mod_potential_and_density, only: rhomt, magmt, rhoir, magir, veffmt, vxcmt, vxcir, bxcmt, bxcir
  use mod_timing, only: timefv, stopwatch
  use modmpi, only: mpiinfo, distribute_loop, mpiglobal, mpi_env_k, mpi_env_band

  implicit none

  private

  !> exciting API to sirius
  public :: setup_sirius, get_mpi_comm_sirius, gengvec_sirius, &
      & setup_sirius_gs_handler, set_radial_functions_sirius, &
      & solve_seceqn_sirius, get_eval_sirius, get_evec_sirius, &
      & finalize_sirius, put_occ_sirius, get_step_function, &
      & generate_density_sirius, get_periodic_function_sirius, &
      & get_max_num_gkvec_sirius, get_gkvec_arrays_sirius, &
      & get_forces_sirius, generate_coulomb_potential_sirius, &
      & set_exchange_correlation_sirius, set_xc_magnetic_sirius, &
      & warn_array_sizes_sirius

  ! Module-scoped variables, required for communicating with the sirius lib
  ! but never exposed to exciting routines. One interacts with them purely
  ! through the API in this module.    
#ifdef SIRIUS
  !> Context handler for SIRIUS input initalisation
  type(sirius_context_handler) :: sctx
  !> Ground state handler of SIRIUS
  type(sirius_ground_state_handler), public  :: gs_handler
  !> K-point set handler of SIRIUS
  type(sirius_kpoint_set_handler), public  :: ks_handler
#endif

contains

  !> Initialise Sirius options with exciting settings.
  !>
  !> Data is set in a module instance of the context handler, sctx
  !> and finally, dumped to a JSON file "setup.json"
  !>
  !> Specifically: 
  !>  - Atoms
  !>  - Basis
  !>  - XC type
  !>  - MPI grid 
  subroutine setup_sirius(input, communicator, mpi_grid)
    use mod_Gvector, only: ngrid
    use constants, only: maxspecies

    !> exciting input options
    type(input_type), intent(in) :: input
    !> MPI communicator for exciting
    integer, intent(in) ::  communicator
    !> 2D Cartesian grid used to split band communicator for the distributed linear algebra oprations
    integer, intent(in) :: mpi_grid(2)

#ifdef SIRIUS
    !> local orbital principal quantum number
    integer :: lorbpqn (maxlorbord, maxlorb, maxspecies)
    integer :: ilo, io, ncls, ia1, ia2
    integer :: is, js, ia
    integer :: ist, l, m, lm, iv
    character(100) :: atom_species_label
    integer, allocatable :: icls(:)
    real(8) :: gkmax

    integer, parameter :: sirius_determine_principal_n = -1
    logical, parameter :: call_mpi_init = .false.

    call sirius_create_context( communicator, sctx )
    call sirius_import_parameters( sctx, &
        '{"parameters" : {"electronic_structure_method" : "full_potential_lapwlo"}, &
          "control"    : {"verification" : 0, &
                          "std_evp_solver_name" : "auto", &
                          "gen_evp_solver_name" : "auto" &
                         }&
         }' )

    call Get_gkmax(input, gkmax)

    call sirius_set_parameters( sctx,&
         lmax_apw=input%groundstate%lmaxapw,&
         lmax_rho=input%groundstate%lmaxvr,&
         lmax_pot=input%groundstate%lmaxvr,&
         auto_rmt=0,&
         num_mag_dims=ndmag,& 
         core_rel=trim(adjustl(input%groundstate%CoreRelativity)),&
         pw_cutoff=input%groundstate%gmaxvr,&
         gk_cutoff=gkmax,&
         fft_grid_size=ngrid,&
         sht_coverage=1,&
         verbosity=3 )

    if (input%groundstate%ValenceRelativity.eq.'iora*') then
      call sirius_set_parameters( sctx, valence_rel='iora' )
    else
      call sirius_set_parameters( sctx, valence_rel=trim(adjustl(input%groundstate%ValenceRelativity)) )
    end if

    call sirius_set_lattice_vectors( sctx,&
         input%structure%crystal%basevect(:,1),&
         input%structure%crystal%basevect(:,2),&
         input%structure%crystal%basevect(:,3) )

    ! Set properties of all atoms
    ! Not clear if the order of the calls in the species loop is important
    ! or one can split this up.         
    do is = 1, nspecies

      atom_species_label = trim(adjustl(spname(is)))

      call sirius_add_atom_type(sctx, atom_species_label, zn=nint(-spzn(is)),&
           &symbol=trim(spsymb(is)), mass=spmass(is))

      call sirius_set_atom_type_radial_grid(sctx, atom_species_label, nrmt(is), spr(1, is))

      ! set elecronic configuration for each atom type (species)
      do ist = 1, spnst(is)
        call sirius_set_atom_type_configuration( sctx, atom_species_label, spn(ist, is), spl(ist, is),&
             &spk(ist, is), spocc(ist, is), spcore(ist, is) )
      end do

      ! set apw basis
      do l = 0, input%groundstate%lmaxapw
        do io = 1, apword(l, is)
          call sirius_add_atom_type_aw_descriptor( sctx, atom_species_label, sirius_determine_principal_n,&
               &l, apwe0(io, l, is), apwdm(io, l, is), apwve(io, l, is) )
        end do
      end do

      ! set local orbital basis
      ! Note, is lorbpqn required? It is not used in exciting
      do ilo = 1, nlorb(is)
        do io = 1, lorbord(ilo, is)
          call sirius_add_atom_type_lo_descriptor( sctx, atom_species_label, ilo, lorbpqn(io, ilo, is),&
               &lorbl(ilo, is), lorbe0(io, ilo, is), lorbdm(io, ilo, is), lorbve(io, ilo, is) )
        end do
      end do

      do ia = 1, natoms(is)
        call sirius_add_atom( sctx, atom_species_label,&
             &input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord,&
             &input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt )
      end do
    end do ! is

    ! set equivalent atoms
    ncls = 0
    allocate(icls(natmtot))
    icls = 0
    do is = 1, nspecies
      do ia1 = 1 ,natoms(is)
        if (icls(idxas(ia1, is)).eq.0) then
          ncls = ncls + 1
          do ia2 = 1, natoms(is)
            if (eqatoms(ia1, ia2, is)) then
              icls(idxas(ia2, is)) = ncls
            end if
          end do
        end if
      end do
    end do
    call sirius_set_equivalent_atoms(sctx, icls)

    ! XC functional
    if (associated(input%groundstate%libxc)) then
      if (input%groundstate%libxc%exchange .ne. 'none') then
        call sirius_add_xc_functional(sctx, trim(input%groundstate%libxc%exchange))
      end if
      if (input%groundstate%libxc%correlation .ne. 'none') then
        call sirius_add_xc_functional(sctx, trim(input%groundstate%libxc%correlation))
      end if
    end if

    ! This is for DEBUG purpose
    call sirius_add_xc_functional(sctx, "XC_GGA_X_PBE")
    call sirius_add_xc_functional(sctx, "XC_GGA_C_PBE")

    ! MPI cart grid
    call sirius_set_mpi_grid_dims(sctx, 2, mpi_grid)

    ! Number of (first variation) eigenstates  
    nstfv = get_nstfv(chgval, input%groundstate%nempty)
    call sirius_set_parameters( sctx, num_fv_states=nstfv )

    ! Iterative ground state solver
    if (input%groundstate%solver%type == 'Davidson') then
      call sirius_set_parameters (sctx, iter_solver_type='davidson')
    end if

    ! Initialize global variables
    call sirius_initialize_context ( sctx )

    if ( mpiglobal%rank == 0 ) then
      call sirius_dump_runtime_setup( sctx, "setup.json" )
    end if

#endif
  end subroutine setup_sirius

  !> Free sirius pointers.
  subroutine finalize_sirius()
#ifdef SIRIUS
    call sirius_free_handler(gs_handler)
    call sirius_free_handler(ks_handler)
    call sirius_free_handler(sctx)
#endif
  end subroutine finalize_sirius

  !> Generate G-vectors using SIRIUS.
  !> 
  !> For SIRIUS only G-vectors within a plane-wave cutoff are generated, in an order
  !> consistent with the distributed FFT. This differs from exciting.
  !> See gencfunc.f90, for example.
  subroutine gengvec_sirius(ivg, ivgig, igfft, vgc, gc, ngvec, intgv)
     !> G-vector integer coordinates
     integer, allocatable, intent(inout) :: ivg(:, :)
     !> Map from integer grid to G-vector array
     integer, allocatable, intent(inout) :: ivgig(:, :, :)
     !> Map from G-vector array to FFT array
     integer, allocatable, intent(inout) :: igfft(:)
     !> G-vectors in Cartesian coordinates
     real(dp), allocatable, intent(inout) :: vgc(:, :)
     !> Length of G-vectors
     real(dp), allocatable, intent(inout) :: gc(:)
     !> Number of G-vectors with G < gmaxvr
     integer, intent(out) :: ngvec
     !> integer grid intervals for each direction
     integer, intent(out) :: intgv (3, 2)
#ifdef SIRIUS
    ! allocate global G-vector arrays
    call sirius_get_num_gvec(sctx, ngvec)

    if (allocated(ivg)) deallocate (ivg)
    allocate (ivg(3, ngvec))

    if (allocated(ivgig)) deallocate (ivgig)
    allocate (ivgig(intgv(1, 1):intgv(1, 2), intgv(2, 1):intgv(2, 2), intgv(3, 1):intgv(3, 2)))

    if (allocated(igfft)) deallocate (igfft)
    allocate (igfft(ngvec))

    if (allocated(vgc)) deallocate (vgc)
    allocate (vgc(3, ngvec))

    if (allocated(gc)) deallocate (gc)
    allocate (gc(ngvec))

    call sirius_get_gvec_arrays(sctx, ivg, vgc, gc, ivgig)
    call sirius_get_fft_index(sctx, igfft)
#endif
  end subroutine gengvec_sirius

  !> Get MPI communicators from SIRIUS
  subroutine get_mpi_comm_sirius(mpi_comm_k, mpi_comm_band)
    !> k-point communicator
    integer, intent(out) :: mpi_comm_k
    !> Bands communicator
    integer, intent(out) :: mpi_comm_band
#ifdef SIRIUS
    call sirius_get_kpoint_inner_comm(sctx, mpi_comm_band)
    call sirius_get_kpoint_inter_comm(sctx, mpi_comm_k)
#endif
  end subroutine get_mpi_comm_sirius

  !> Return the number of k-points per process.
  subroutine k_points_per_process(mpi_env_kpoints, nkpt, k_counts)
#ifdef MPI    
    use modmpi, only: MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
#endif
    !> MPI communicator for k-points
    type(mpiinfo), intent(in) :: mpi_env_kpoints
    !> Total number of k-points
    integer, intent(in) :: nkpt
    !> Number of k-points per process of mpi_env_kpoints
    integer, intent(out) :: k_counts(:)

    integer :: firstk, lastk, ierr
#ifdef MPI
    call distribute_loop(mpi_env_kpoints, nkpt, firstk, lastk)
    k_counts = 0
    k_counts(mpi_env_kpoints%rank + 1) = lastk - firstk + 1
    call mpi_allreduce( MPI_IN_PLACE, k_counts, mpi_env_kpoints%procs, &
                      & MPI_INTEGER, MPI_SUM, mpi_env_kpoints%comm, ierr)
#else
    k_counts(1) = nkpt
#endif
  end subroutine

  !> Set a ground state executor, with a given k-set.
  subroutine setup_sirius_gs_handler( input )
    !> exciting input options
    type(input_type), intent(in) :: input
    !> Number of k-points per process
    integer, allocatable :: counts(:)
#ifdef SIRIUS
    !> Indicates whether the kset will be initialized by sirius
    logical, parameter :: initialised_kset = .false.

    allocate(counts(mpi_env_k%procs))
    call k_points_per_process(mpi_env_k, nkpt, counts)

    ! Create a new k-set for density generation
    call sirius_create_kset( sctx, nkpt, vkl, wkpt, initialised_kset, ks_handler)
    call sirius_initialize_kset( ks_handler, count=counts)
    ! Create a ground state executor with the given k-point set
    call sirius_create_ground_state( ks_handler, gs_handler)

    ! set MT pointer aliases to avoid allocating memory twice
    call sirius_set_periodic_function_ptr( gs_handler, "rho", f_mt=rhomt )
    call sirius_set_periodic_function_ptr( gs_handler, "veff", f_mt=veffmt )
    if (ndmag.eq.1) then
      call sirius_set_periodic_function_ptr( gs_handler, "magz", f_mt=magmt(:, :, :, 1) )
    end if
    if (ndmag.eq.3) then
      call sirius_set_periodic_function_ptr( gs_handler, "magx", f_mt=magmt(:, :, :, 1) )
      call sirius_set_periodic_function_ptr( gs_handler, "magy", f_mt=magmt(:, :, :, 2) )
      call sirius_set_periodic_function_ptr( gs_handler, "magz", f_mt=magmt(:, :, :, 3) )
    end if
    deallocate(counts)
#endif

  end subroutine setup_sirius_gs_handler

  !> Set Sirius radial functions
  subroutine set_radial_functions_sirius(input)
    use mod_APW_LO, only: apwfr, lofr
    !> exciting input options
    type(input_type), intent(in) :: input
#ifdef SIRIUS
    integer is, ia, ias, l, i, enu_deriv

    do is = 1, nspecies
      do ia = 1, natoms (is)
        ias = idxas (ia, is)
        do l = 0, input%groundstate%lmaxapw
          do i = 1, apword (l, is)
            do enu_deriv = 0, 1
              call sirius_set_radial_function( sctx, ias, enu_deriv, apwfr(:, enu_deriv + 1, i, l, ias), l=l, o=i )
            end do
          end do
        end do
        do i = 1, nlorb (is)
          do enu_deriv = 0, 1
            call sirius_set_radial_function( sctx, ias, enu_deriv, lofr(:, enu_deriv + 1, i, ias), ilo=i )
          end do
        end do
      end do !ia
    end do !is
#endif
  end subroutine set_radial_functions_sirius

  !> Call SIRIUS to solve band diagonalization problem.
  !>
  !> The result is returned on the sirius side, and must be accessed with
  !> the `sirius_get_band_energies` routine.
  subroutine solve_seceqn_sirius()
    use mod_potential_and_density, only: veffir
#ifdef SIRIUS
    call sirius_set_periodic_function(gs_handler, "veff", f_rg=veffir, f_rg_global=.true.)
    call sirius_find_eigen_states(gs_handler, ks_handler, precompute_pw=.true., &
                                  precompute_rf=.false., precompute_ri=.true. )
#endif
  end subroutine solve_seceqn_sirius

  !> Get eigenvalues from SIRIUS and write to file.
  !>
  !> As in all of the ground state, the (first variation) eigenvalues and 
  !> eigenvectors are never stored in memory w.r.t k-points, and instead are
  !> written to file.
  !>
  !> Second-variation eigenvalues are used from their module declaration
  !> as they are likely already allocated in init1.
  !> Unlike fv eigenvalues, this array does have dims of (nstfv, nkpt)
  !> and so must be treated differently.
  subroutine get_eval_sirius()
    !> First variation eigenvalues
    real(dp), allocatable :: evalfv(:, :)
    !> k-indices
    integer :: ik, firstk, lastk

    call stopwatch("exciting:get_eval_sirius", 1)
#ifdef SIRIUS

    ! Get second-variational eigenvalues from sirius 
    if (ndmag == 0 .or. ndmag == 3) then
       do ik = 1, nkpt
          call sirius_get_band_energies( ks_handler, ik, 1, evalsv(:, ik) )
       end do
    else
      do ik = 1, nkpt
          call sirius_get_band_energies( ks_handler, ik, 1, evalsv(:, ik) )
          call sirius_get_band_energies( ks_handler, ik, 2, evalsv(nstfv + 1:nstsv, ik) )
      enddo
    endif

    ! Write to first+second variational eigenvalues to file, per k-point 
    allocate(evalfv(nstfv, nspnfv))
    call distribute_loop( mpi_env_k, nkpt, firstk, lastk )
    
    if (mpi_env_band%is_root) then
       do ik = firstk, lastk
            ! Get first variation eigenvalues from sirius
           call sirius_get_fv_eigen_values( ks_handler, ik, evalfv(:, 1), nstfv )
           call putevalfv( ik, evalfv(:, :) )
           call putevalsv( ik, evalsv(:, ik) )
       end do
    end if

#endif
    call stopwatch("exciting:get_eval_sirius", 0)
  end subroutine get_eval_sirius

  !> Generate the density with Sirius.
  !>
  !> The gs_handler is passed via this module as its type
  !> is defined on the library side, therefore one would need a
  !> a mock type (for when sirius is not present) if it's to be 
  !> exposed elsewhere in exciting.
  subroutine generate_density_sirius()
#ifdef SIRIUS
    call sirius_generate_density(gs_handler, .false., .true.)
#endif
  end subroutine 

  subroutine get_periodic_function_sirius(rhoir, ngrtot)
    !> Interstitial real-space charge density
    real(dp), intent(inout) :: rhoir(:)     
    !> Total number of real space grid points
    !> One wonders if this is just size(rhoir)
    integer, intent(inout) :: ngrtot
#ifdef SIRIUS
    call sirius_get_periodic_function(gs_handler, "rho", f_rg=rhoir, num_rg_points=ngrtot)
#endif
  end subroutine

  !> Set band/state occupations in sirius.
  !>
  !> Note, the sv (second-variation) naming convention 
  !> of the occupations is misleading.
  subroutine put_occ_sirius() 
    !> Second-variation state occupations
    !> Compilation complained when explicitly passed - revisit now the API works
    ! real(dp), intent(inout) :: occsv(:, :)
    integer :: ik
#ifdef SIRIUS
    do ik = 1, size(occsv, 2)
      if (ndmag.eq.0.or.ndmag.eq.3) then
        call sirius_set_band_occupancies(ks_handler, ik, 1, occsv(:, ik))
      else
        call sirius_set_band_occupancies(ks_handler, ik, 1, occsv(:, ik))
        call sirius_set_band_occupancies(ks_handler, ik, 2, occsv(nstfv+1:nstsv, ik))
      end if
    end do
#endif
  end subroutine put_occ_sirius

  !> Get eigenvectors from SIRIUS and write to file.
  subroutine get_evec_sirius()
    complex(dp), allocatable :: evecfv(:, :, :)
    complex(dp), allocatable :: evecsv(:, :) 
    integer :: ik, firstk, lastk
    call stopwatch("exciting:get_evec_sirius", 1)
#ifdef SIRIUS
    allocate(evecfv(nmatmax, nstfv, nspnfv))
    allocate(evecsv(nstsv, nstsv))

    call distribute_loop( mpi_env_k, nkpt, firstk, lastk )

    do ik = firstk, lastk
      call sirius_get_fv_eigen_vectors( ks_handler, ik, evecfv(:, :, 1), nmatmax, nstfv )
      call sirius_get_sv_eigen_vectors( ks_handler, ik, evecsv, nstsv )
      if ( mpi_env_band%is_root ) then
        call putevecfv( ik, evecfv )
        call putevecsv( ik, evecsv )
      end if
    end do
#endif
    call stopwatch("exciting:get_evec_sirius", 0)
  end subroutine get_evec_sirius

  !> Get step function from SIRIUS
  subroutine get_step_function( cfunig, cfunir, ngrtot )
    !> G-space characteristic function: 0 inside the muffin-tins and 1 outside
    complex(dp), intent(out) :: cfunig(:)
    !> Real-space characteristic function: 0 inside the muffin-tins and 1 outside
    real(dp), intent(out) :: cfunir(:)
    !> Total number of G-vectors
    integer, intent(in) :: ngrtot
#ifdef SIRIUS
    call sirius_get_step_function( sctx, cfunig, cfunir, ngrtot )
#endif
  end subroutine get_step_function

  !> Get maximum number of Gk vectors
  subroutine get_max_num_gkvec_sirius(ngkmax)
    integer, intent(out) :: ngkmax        
#ifdef SIRIUS    
    call sirius_get_max_num_gkvec(ks_handler, ngkmax)
#else 
    ! Dummy value
    ngkmax = 0
#endif     
  end subroutine 

  ! TODO(Alex) Review intent of arguments.
  !> Get Gk vector arrays for a given k-point and spin.
  subroutine get_gkvec_arrays_sirius(ik, ngk, igkig, vgkl, vgkc, gkc, tpgkc)
    !> k-index
    integer, intent(in) :: ik
    !> Number of G+k-vectors for augmented plane waves
    integer, intent(inout)  :: ngk
    !> Index from G+k-vectors to G-vectors
    integer, intent(inout), contiguous :: igkig(:)
    !> G+k-vectors in lattice coordinates
    real(dp), intent(inout) :: vgkl(:, :)
    !> G+k-vectors in Cartesian coordinates
    real(dp), intent(inout) :: vgkc(:, :)
    !> Length of G+k-vectors
    real(dp), intent(inout), contiguous :: gkc(:)
    !> (theta, phi) coordinates of G+k-vectors
    real(dp), intent(inout) :: tpgkc(:, :)
#ifdef SIRIUS  
    call sirius_get_gkvec_arrays(ks_handler, ik, ngk, igkig, vgkl, vgkc, gkc, tpgkc)
#endif 
  end subroutine

  !> Get forces from sirius
  subroutine get_forces_sirius(force_type, forces)
    !> Type of force requested
    character(len=*), intent(in) :: force_type
    !> Force on each atom
    real(dp), intent(inout) :: forces (:, :)
#ifdef SIRIUS    
    !NOTE, was passing as forces(1,1)
    call sirius_get_forces(gs_handler, trim(adjustl(force_type)), forces)
#endif    
  end subroutine

  ! TODO(Alex) Review intent of arguments.
  !> Generate Coulomb Potential for MT and Interstitial Regions.
  !>
  !> Total number of atoms `natmtot` accessed as a global (already loaded in this module).
  subroutine generate_coulomb_potential_sirius(lmmaxvr, nrmtmax, natmmax, ngrtot, rhoir, vmad, vclmt, vclir)
    !> Interstitial real-space charge density
    real(dp), intent(in) :: rhoir(:)  
    !> Madulung potential (excluding the on-site nuclear contribution) 
    !> in the innermost radial point
    real(dp), intent(inout) :: vmad(:)
    !> Muffin-tin Coulomb potential
    real(dp), intent(inout) :: vclmt(:, :, :)
    !> Interstitial real-space Coulomb potential
    real(dp), intent(inout) :: vclir(:)

    !> Maximum angular momentum for potentials and densities.
    integer, intent(in) :: lmmaxvr
    !> Maximum nrmt (radial MT sampling points) over all the species.
    integer, intent(in) :: nrmtmax
    !> Maximum number of atoms over all the species.
    integer, intent(in)  :: natmmax
    !> Total number of G-vectors.
    integer, intent(in)  :: ngrtot

#ifdef SIRIUS
    ! Note, was passing as rhoir(1)
    call sirius_set_periodic_function(gs_handler, "rho", f_rg=rhoir, f_rg_global=.true.)
    ! Note, was passing as vmad(1)
    call sirius_generate_coulomb_potential(gs_handler, vh_el=vmad)
    ! Note, was passing as vclir(1)
    call sirius_get_periodic_function(gs_handler, "vha", f_mt=vclmt, lmmax=lmmaxvr, &
                                      max_num_mt_points=nrmtmax, num_atoms=natmtot, &
                                      f_rg=vclir, num_rg_points=ngrtot)
#endif
  end subroutine

  !> Set exchange and correlation potential for MT and interstitial regions.
  !>
  !> Total number of atoms `natmtot` accessed as a global (already loaded in this module).
  subroutine set_exchange_correlation_sirius(lmmaxvr, nrmtmax, ngrtot, rhoir, vxcmt, vxcir, exmt, exir, ecmt, ecir)

    !> Maximum angular momentum for potentials and densities
    integer, intent(in) :: lmmaxvr
    !> Maximum nrmt (radial MT sampling points) over all the species
    integer, intent(in) :: nrmtmax
    !> Total number of G-vectors
    integer, intent(in) :: ngrtot
    
    !> Interstitial real-space charge density
    real(dp), intent(in) :: rhoir(:)  
    !> Muffin-tin exchange-correlation potential
    real(dp), intent(out) :: vxcmt(:, :, :)
    !> Interstitial real-space exchange-correlation potential
    real(dp), intent(out) :: vxcir(:)

    !> Muffin-tin exchange energy density
    real(dp), intent(out):: exmt (:, :, :)
    !> Interstitial real-space exchange energy density
    real(dp), intent(out) :: exir(:)
    !> Muffin-tin correlation energy density
    real(dp), intent(out) :: ecmt(:, :, :)
    !> Interstitial real-space correlation energy density
    real(dp), intent(out)  :: ecir(:)

#ifdef SIRIUS
    call sirius_set_periodic_function(gs_handler, "rho", f_rg=rhoir, f_rg_global=.true.)
    call sirius_generate_xc_potential(gs_handler)
    call sirius_get_periodic_function(gs_handler, "vxc", f_mt=vxcmt, lmmax=lmmaxvr, &
                                      max_num_mt_points=nrmtmax, num_atoms=natmtot, &
                                      f_rg=vxcir, num_rg_points=ngrtot)
    call sirius_get_periodic_function(gs_handler, "exc", f_mt=exmt, lmmax=lmmaxvr, &
                                      max_num_mt_points=nrmtmax, num_atoms=natmtot,&
                                      f_rg=exir, num_rg_points=ngrtot)
    ecmt = 0.d0
    ecir = 0.d0
#endif
  end subroutine


  !> Set exchange and correlation potential for MT and interstitial regions, under
  !> the influence of a magnetic field.
  !>
  !> Total number of atoms `natmtot` accessed as a global (already loaded in this module).
  !> 
  !> TODO(Alex) Refactor the implementation. One could do this more cleanly looping:
  !>
  !> if (ndmag == 1) field = ['bz']
  !> if (ndmag == 3) field = ['bx', 'by', 'bz']
  !> do i = 1 , ndmag
  !>   call sirius_get_periodic_function(gs_handler, field(i), f_mt=bxcmt(:, :, :, i),&
  !>                                     lmmax=lmmaxvr, max_num_mt_points=nrmtmax, num_atoms=natmtot,&
  !>                                     f_rg=bxcir(:, i), num_rg_points=ngrtot)
  !> endo
  subroutine set_xc_magnetic_sirius(ndmag, lmmaxvr, nrmtmax, ngrtot, bxcmt, bxcir)
    !> Dimension of magnetisation and magnetic vector fields (1 or 3)
    integer, intent(in) :: ndmag
    !> Maximum angular momentum for potentials and densities
    integer, intent(in) :: lmmaxvr
    !> Maximum nrmt (radial MT sampling points) over all the species
    integer, intent(in) :: nrmtmax
    !> Total number of G-vectors
    integer, intent(in) :: ngrtot
    !> Muffin-tin exchange-correlation magnetic field
    real(dp), intent(out) :: bxcmt (:, :, :, :)
    !> Interstitial exchange-correlation magnetic field
    real(dp), intent(out) :: bxcir (:, :)

#ifdef SIRIUS
  if (ndmag == 1) then
    call sirius_get_periodic_function(gs_handler, "bz", f_mt=bxcmt(:, :, :, 1),&
                                      lmmax=lmmaxvr, max_num_mt_points=nrmtmax, num_atoms=natmtot,&
                                      f_rg=bxcir(:, 1), num_rg_points=ngrtot)
  end if

  if (ndmag == 3) then
    call sirius_get_periodic_function(gs_handler, "bx", f_mt=bxcmt(:, :, :, 1),&
                                      lmmax=lmmaxvr, max_num_mt_points=nrmtmax, num_atoms=natmtot,&
                                      f_rg=bxcir(:, 1), num_rg_points=ngrtot)
    call sirius_get_periodic_function(gs_handler, "by", f_mt=bxcmt(:, :, :, 2),&
                                      lmmax=lmmaxvr, max_num_mt_points=nrmtmax, num_atoms=natmtot,&
                                      f_rg=bxcir(:, 2), num_rg_points=ngrtot)
    call sirius_get_periodic_function(gs_handler, "bz", f_mt=bxcmt(:, :, :, 3),&
                                      lmmax=lmmaxvr, max_num_mt_points=nrmtmax, num_atoms=natmtot,&
                                      f_rg=bxcir(:, 3), num_rg_points=ngrtot)
  end if
#endif
  end subroutine


  !> Check array sizes of spherical harmonics and structure factor for the G-vectors.
  !>
  !> ylmg is a huge, non distributed array on Exciting side. 
  !> It is required for the generation of initial density and during the SCF, 
  !> for generating the Hartree potential. If not distributed (requiring a refactor)
  !> this and sfacg become memory bottlenecks. 
  ! sfacg is even worse than ylmg, because it scales as natomtot^2
  !> This routine indicates to the user when their memory requirements become large.
  !>
  !> Total number of atoms `natmtot` accessed as a global (already loaded in this module).
  subroutine warn_array_sizes_sirius(lmmaxvr, ngvec)
    !> Maximum angular momentum for potentials and densities
    integer, intent(in) :: lmmaxvr
    !> Number of G-vectors with G < gmaxvr
    integer, intent(in) :: ngvec

    real(dp) :: mb, ylmg_mb, sfacg_mb

    mb = 1._dp / 1024._dp / 1024._dp
    ylmg_mb = mb * 16 * lmmaxvr * ngvec
    sfacg_mb = mb * 16 * natmtot * ngvec

    if (ylmg_mb > 1000._dp) then
        write(*, *)'Warning: large ylmg array of ',ylmg_mb, 'Mb'
    end if
    if (sfacg_mb > 1000._dp) then
        write(*,*)'Warning: large sfacg array of ',sfacg_mb, 'Mb'
    end if

  end subroutine

end module sirius_api
