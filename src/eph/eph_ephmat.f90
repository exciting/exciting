!> Module handling electron-phonon matrix elements.
module eph_ephmat
  use eph_variables

  use precision, only: dp
  use asserts, only: assert
  use modmpi
  use block_data_file, only: block_data_file_type

  implicit none
  private

  ! standard (Fan-Migdal) matrix elements
  !> name for binary file to save standard (Fan-Migdal) EPH matrix in reciprocal space phonon Hamiltonian gauge for later access
  character(*), parameter :: eph_ephmat_gFMkq_filename = "EPH_GKQ.OUT"
  !> name for binary file to save standard (Fan-Migdal) EPH matrix in real space atomic Wannier gauge for later access
  character(*), parameter :: eph_ephmat_gFMRR_filename = "EPH_GRR.OUT"
  
  public :: eph_ephmat_gen_coarse, eph_ephmat_setup_interpolation, eph_ephmat_interpolate, eph_ephmat_average

contains

  !================================================================================ 
  ! GENERATE MATRIX ELEMENTS ON COARSE GRIDS
  !
  !> Generate the electron-phonon matrix elements in atomic Hamiltonian gauge
  !> on the coarse electron grid \({\bf k}\) and phonon grid \({\bf q}\)
  !>
  !> \[ g_{mn,\kappa \alpha}({\bf k},{\bf q}) = 
  !>    \langle \psi_{m{\bf k}+{\bf q}} | \delta^{\bf q}_{\kappa \alpha} V({\bf r}) | \psi_{n{\bf k}} \rangle \;. \]
  subroutine eph_ephmat_gen_coarse( standard, mpicomm )
    use constants, only: zzero
    use dfpt_variables, only: dfpt_lmaxvr, dfpt_lmmaxvr, dfpt_kset, dfpt_Gset, dfpt_Gkset, mt_basis, fevalk0, feveck0
    use dfpt_eigensystem, only: dfpt_eig_geteval, dfpt_eig_getevec
    use matrix_elements, only: me_init, me_finit, me_mt_alloc
    use mod_kpointset
    use mod_atoms, only: natmtot
    use mod_muffin_tin, only: nrmtmax
    use mod_APW_LO, only: nlotot
    !> compute standard e-ph matrix elements (default: `.true.`)
    logical, optional, intent(in) :: standard
    !> MPI communicator (default: global MPI communicator)
    type(mpiinfo), optional, intent(inout) :: mpicomm

    integer :: iq, ik, irec, ias, ip, imode, iq_range(2), nmat
    logical :: do_standard
    type(mpiinfo) :: mpi
    type(block_data_file_type) :: gFMkq_file
    type(k_set) :: kqset
    type(G_set) :: Gqset
    type(Gk_set) :: Gkqset

    integer, allocatable :: ik_range(:,:)
    real(dp), allocatable :: evalk(:), evalkq(:)
    complex(dp), allocatable :: gpot_mt_basis(:,:,:,:), &
                                dpot_mt(:,:,:,:,:), dpot_ir(:,:,:), dpot_mt_basis(:,:,:,:), dpot_cfun_ig(:,:), &
                                gFM_atomic(:,:,:), &
                                eveck(:,:), eveckq(:,:)

    ! set defaults
    do_standard = .true.
    if (present(standard)) do_standard = standard
    mpi = mpiglobal
    if (present(mpicomm)) mpi = mpicomm

    ! initialize matrix elements module
    call me_init( mt_basis, dfpt_lmaxvr, dfpt_Gset )

    ! set up interstitial basis
    if (.not. allocated(eph_Gkset_el%ngk)) &
      call generate_Gk_vectors( eph_Gkset_el, eph_kset_el, dfpt_Gset, dfpt_Gkset%gkmax )

    ! allocate arrays
    allocate( dpot_mt(dfpt_lmmaxvr, nrmtmax, natmtot, 3, natmtot) )
    allocate( dpot_ir(dfpt_Gset%ngrtot, 3, natmtot) )

    !******************************************************************************** 
    ! standard electron-phonon matrix elements in atomic Hamiltonian gauge
    !******************************************************************************** 
    ! Note: Called `gFM` due their use in the Fan-Migdal self-energy and to distinguish 
    !       them from the Debye-Waller matrix elements.
    if (do_standard) then
      ! set up and open binary files
      gFMkq_file = block_data_file_type( eph_ephmat_gFMkq_filename, [eph_nst, eph_nst, 3*natmtot], cmplx( 0, 0, dp ) )
      call gFMkq_file%open( mpi, delete_existing=.true. )

      ! set loop limits
      call distribute_double_loop( mpi%rank, mpi%procs, [1, eph_qset_ph%nkpt], [1, eph_kset_el%nkpt], iq_range, ik_range )
      
      ! allocate arrays
      allocate( gFM_atomic(eph_fst:eph_lst, eph_fst:eph_lst, 3*natmtot) )
      call me_mt_alloc( dpot_mt_basis, 3*natmtot )

      ! generate radial integrals times Gaunt coefficients for potential gradient
      ! using Gauss' theorem and regularized effective potential
      call eph_ephmat_gen_grad_pot_mt_basis( dfpt_lmaxvr, gpot_mt_basis )
      
      ! compute matrix elements
      do iq = iq_range(1), iq_range(2)
        ! generate k+q, G+q, and G+k+q vectors
        kqset = eph_kset_el
        do ik = 1, eph_kset_el%nkpt
          kqset%vkl(:, ik) = kqset%vkl(:, ik) + eph_qset_ph%vkl(:, iq)
          kqset%vkc(:, ik) = kqset%vkc(:, ik) + eph_qset_ph%vkc(:, iq)
        end do
        call generate_G_vectors( Gqset, dfpt_Gset%bvec, dfpt_Gset%intgv, dfpt_Gset%gmaxvr, vpl=eph_qset_ph%vkl(:, iq) )
        call generate_Gk_vectors( Gkqset, kqset, dfpt_Gset, dfpt_Gkset%gkmax )
        ! allocate arrays
        allocate( dpot_cfun_ig(Gqset%ngvec, 3*natmtot) )
        ! read potential response from file
        call eph_ephmat_read_dpot( eph_qset_ph, eph_qset_ph%vkl(:, iq), dfpt_lmaxvr, dpot_mt, Gqset, dpot_ir )
        ! prepare k-independent part of matrix elements
        do ias = 1, natmtot
          do ip = 1, 3
            imode = (ias -1) * 3 + ip
            call eph_ephmat_prepare_gmat( dfpt_lmaxvr, dpot_mt(:, :, :, ip, ias), Gqset, dpot_ir(:, ip, ias), &
              dpot_mt_basis(:, :, :, imode), dpot_cfun_ig(:, imode) )
            ! add gradient contribution
            dpot_mt_basis(:, :, ias, imode) = dpot_mt_basis(:, :, ias, imode) + gpot_mt_basis(:, :, ias, ip)
          end do
        end do
        ! compute matrix elements for all k-points
        do ik = ik_range(1, iq), ik_range(2, iq)
          nmat = Gkqset%ngk(1, ik) + nlotot
          ! read eigenvalues and eigenvectors at k
          call dfpt_eig_geteval( eph_kset_el%vkl(:, ik), fevalk0, dfpt_kset, [eph_fst, eph_lst], evalk )
          call dfpt_eig_getevec( eph_kset_el%vkl(:, ik), eph_Gkset_el%vgkl(:, :, 1, ik), feveck0, dfpt_kset, dfpt_Gkset, [eph_fst, eph_lst], eveck )
          ! read eigenvalues and eigenvectors at k+q
          call dfpt_eig_geteval( kqset%vkl(:, ik), fevalk0, dfpt_kset, [eph_fst, eph_lst], evalkq )
          call dfpt_eig_getevec( kqset%vkl(:, ik), Gkqset%vgkl(:, :, 1, ik), feveck0, dfpt_kset, dfpt_Gkset, [eph_fst, eph_lst], eveckq )
      
          ! STANDARD MATRIX ELEMENTS
          do ias = 1, natmtot
            do ip = 1, 3
              imode = (ias - 1) * 3 + ip
              call eph_ephmat_gen_gmat( ik, eph_Gkset_el, Gkqset, Gqset, eph_fst, eph_lst, eph_fst, eph_lst, &
                dpot_mt_basis(:, :, :, imode), dpot_cfun_ig(:, imode), eveckq, eveck, gFM_atomic(:, :, imode) )
            end do
          end do
      
          ! write eph matrix elements in atomic coordinates
          irec = (iq - 1) * eph_kset_el%nkpt + ik
          call gFMkq_file%write( irec, gFM_atomic )
        end do
        ! deallocate arrays
        deallocate( dpot_cfun_ig )
      end do

      deallocate( gFM_atomic )
      
      ! close binary file
      call gFMkq_file%close( mpi )
    end if

    ! free memory
    deallocate( dpot_mt, dpot_ir )
    if (allocated(dpot_mt_basis)) deallocate( dpot_mt_basis )
    call delete_k_vectors( kqset )
    call delete_G_vectors( Gqset )
    call delete_Gk_vectors( Gkqset )
    call me_finit
  end subroutine eph_ephmat_gen_coarse
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! SETUP WANNIER-FOURIER INTERPOLATION OF EPH MATRIX
  !
  !> Set up the Wannier-Fourier interpolation of the electron-phonon matrix.
  !>
  !> This includes
  !> 
  !>   * reading \(g_{mn,\kappa\alpha}({\bf k}, {\bf q})\) in atomic Hamiltonian gauge on coarse grids
  !>   * transformation to atomic Wannier gauge \(\mathcal{g}_{mn,\kappa\alpha}({\bf k}, {\bf q})\)
  !>   * subtraction of long-range Fröhlich coupling for polar materials \(\mathcal{\bf g}^\mathcal{S} = \mathcal{\bf g} - \mathcal{\bf g}^\mathcal{L}\)
  !>   * double Fourier transform to real space \(\mathcal{g}^\mathcal{S}_{mn,\kappa\alpha}({\bf R}_e, {\bf R}_{ph})\)
  subroutine eph_ephmat_setup_interpolation( write_localization )
    use eph_electrons, only: eph_el_mfi, eph_el_evec_k
    use eph_phonons, only: eph_ph_mfi
    use constants, only: zone
    use mod_atoms, only: natmtot
    use modinput
    !> write spatial localization of \(\mathcal{\bf g}^\mathcal{S}({\bf R}_e, 0)\) and \(\mathcal{\bf g}^\mathcal{S}(0, {\bf R}_{ph})\) to file (default: `.false.`)
    logical, optional, intent(in) :: write_localization

    integer :: iq, iq1, iq2, ik, ik0, ikq0, ire, ire1, ire2, irec, irp, isymk, isymkq, ivg(3), un
    logical :: write_loc
    real(dp) :: vkl(3), vkql(3)
    type(block_data_file_type) :: gFMkq_file, gFMRq_file, gFMRR_file

    complex(dp), allocatable :: gFM_aW_RR(:,:,:,:), gFM_aH_kq(:,:,:), gFM_aW_kq(:,:,:,:), gFM_aW_Rq(:,:,:,:), &
                                g_lr_pref_a(:,:)

    write_loc = .false.
    if (present(write_localization)) write_loc = write_localization

    !******************************************************************************** 
    ! standard electron-phonon matrix elements
    !******************************************************************************** 
    ! set up binary files
    gFMRR_file = block_data_file_type( eph_ephmat_gFMRR_filename, [eph_nwf_tot, eph_nwf_tot, 3*natmtot], cmplx( 0, 0, dp ) )
    gFMkq_file = block_data_file_type( eph_ephmat_gFMkq_filename, [eph_nst, eph_nst, 3*natmtot], cmplx( 0, 0, dp ) )

    ! try to read last g_aW(R,R) from file
    if (gFMRR_file%exists()) then
      allocate( gFM_aW_RR(eph_nwf_tot, eph_nwf_tot, 3*natmtot, 1) )
      irec = eph_el_mfi%nr * eph_ph_mfi%nr
      call gFMRR_file%open( mpiglobal )
      call gFMRR_file%read( irec, gFM_aW_RR(:, :, :, 1) )
      call gFMRR_file%close( mpiglobal )
      deallocate( gFM_aW_RR )
    ! compute g_aW(R,R) and write to file
    else if (gFMkq_file%exists()) then
      ! open file for g_aH(k,q)
      call gFMkq_file%open( mpiglobal )

      ! create and open temporary file for g_aW(Re,q)
      gFMRq_file = block_data_file_type( 'EPH_GRQ.TMP', [eph_nwf_tot, eph_nwf_tot, 3*natmtot], cmplx( 0, 0, dp ) )
      call gFMRq_file%open( mpiglobal, delete_existing=.true. )

      !................................................................................ 
      ! Fourier transform k -> Re
      !
      allocate( gFM_aH_kq(eph_nst, eph_nst, 3*natmtot) )
      allocate( gFM_aW_kq(eph_nwf_tot, eph_nwf_tot, 3*natmtot, eph_el_mfi%np) ) 
      allocate( gFM_aW_Rq(eph_nwf_tot, eph_nwf_tot, 3*natmtot, eph_el_mfi%nr) ) 
      if (eph_polar) allocate( g_lr_pref_a(3, natmtot) )

      ! set loop limits
      iq1 = firstofset( mpiglobal%rank, eph_ph_mfi%np, mpiglobal%procs )
      iq2 = lastofset( mpiglobal%rank, eph_ph_mfi%np, mpiglobal%procs )

      do iq = iq1, iq2
        ! generate prefactors for long-range matrix elements in atomic gauge
        if (eph_polar) then
          if (input%eph%elphbolt) then
            call eph_ephmat_gen_lr_prefactor_atomic_elphbolt( eph_ph_mfi%bvec, eph_ph_mfi%vpl(:, iq), eph_dielten, eph_borncharge, g_lr_pref_a, eph_ph_mfi%ngrid, 14.0_dp )
          else
            call eph_ephmat_gen_lr_prefactor_atomic( eph_ph_mfi%bvec, eph_ph_mfi%vpl(:, iq), eph_dielten, eph_borncharge, g_lr_pref_a )
          end if
        end if
        do ik = 1, eph_el_mfi%np
          ! find k in electronic k-point set
          vkl = eph_el_mfi%vpl(:, ik)
          call r3frac( 1e-6_dp, vkl, ivg )
          call findkptinset( vkl, eph_kset_el, isymk, ik0 )
          ! find k+q in electronic k-point set
          vkql = eph_el_mfi%vpl(:, ik) + eph_ph_mfi%vpl(:, iq)
          call r3frac( 1e-6_dp, vkql, ivg )
          call findkptinset( vkql, eph_kset_el, isymkq, ikq0 )
          ! read g_aH(k,q) from file
          call eph_ephmat_read_coarse( eph_kset_el, eph_qset_ph, eph_el_mfi%vpl(:, ik), eph_ph_mfi%vpl(:, iq), gFMkq_file, gFM_aH_kq )
          ! transform to Wannier gauge
          call eph_ephmat_transform_Hamiltonian_Wannier( eph_el_evec_k(:, :, ik0), eph_el_evec_k(:, :, ikq0), gFM_aH_kq, gFM_aW_kq(:, :, :, ik), 3*natmtot, 1 )
          ! subtract long range matrix elements in atomic Wannier gauge
          if (eph_polar) &
            call eph_ephmat_add_lr( -zone, g_lr_pref_a, gFM_aW_kq(:, :, :, ik) )
        end do
        ! transform from k to Re
        call eph_el_mfi%transform_p2R( [eph_nwf_tot, eph_nwf_tot], 3*natmtot, gFM_aW_kq, eph_nwf_tot**2, 3*natmtot, gFM_aW_Rq, eph_nwf_tot**2, 3*natmtot )
        ! write g_aW(Re,q) to temporary file
        do ire = 1, eph_el_mfi%nr
          irec = (iq - 1) * eph_el_mfi%nr + ire
          call gFMRq_file%write( irec, gFM_aW_Rq(:, :, :, ire) )
        end do
      end do

      deallocate( gFM_aH_kq, gFM_aW_kq, gFM_aW_Rq )
      if (eph_polar) deallocate( g_lr_pref_a )
      !................................................................................ 

      ! close file for g_aH(k,q)
      call gFMkq_file%close( mpiglobal )

      ! open file for g_aW(Re,Rp)
      call gFMRR_file%open( mpiglobal )

      !................................................................................ 
      ! Fourier transform q -> Rp
      !
      allocate( gFM_aW_Rq(eph_nwf_tot, eph_nwf_tot, 3*natmtot, eph_ph_mfi%np) )
      allocate( gFM_aW_RR(eph_nwf_tot, eph_nwf_tot, 3*natmtot, eph_ph_mfi%nr) )

      ! set loop limits
      ire1 = firstofset( mpiglobal%rank, eph_el_mfi%nr, mpiglobal%procs )
      ire2 = lastofset( mpiglobal%rank, eph_el_mfi%nr, mpiglobal%procs )

      do ire = ire1, ire2
        ! read g_aW(Re,q) from temporary file
        do iq = 1, eph_ph_mfi%np
          irec = (iq - 1) * eph_el_mfi%nr + ire
          call gFMRq_file%read( irec, gFM_aW_Rq(:, :, :, iq) )
        end do
        ! transform from q to Rp
        call eph_ph_mfi%transform_p2R( [eph_nwf_tot, eph_nwf_tot], 3*natmtot, gFM_aW_Rq, eph_nwf_tot**2, 3*natmtot, gFM_aW_RR, eph_nwf_tot**2, 3*natmtot )
        ! write g_aW(Re,Rp) to file
        do irp = 1, eph_ph_mfi%nr
          irec = (irp - 1) * eph_el_mfi%nr + ire
          call gFMRR_file%write( irec, gFM_aW_RR(:, :, :, irp) )
        end do
      end do
      !................................................................................ 

      ! close and delete temporary file for g_aW(Re,q)
      call gFMRq_file%close( mpiglobal )
      call gFMRq_file%delete( mpiglobal )

      ! write spatial localization to file
      if (mpiglobal%rank == 0 .and. write_loc) then
        open( newunit=un, file='eph_ephmat_loc_el.dat', action='write', form='formatted' )
        write( un, '("#",a5,a26,a26)' ) 'iR', '|Re|', 'max(|g(Re,0)|)'
        irp = 1
        do ire = 1, eph_el_mfi%nr
          irec = (irp - 1) * eph_el_mfi%nr + ire
          call gFMRR_file%read( irec, gFM_aW_RR(:, :, :, irp) )
          write( un, '(i6,2g26.16)' ) ire, eph_el_mfi%rlen(ire), maxval( abs( gFM_aW_RR(:, :, :, irp) ) )
        end do
        close( un )

        open( newunit=un, file='eph_ephmat_loc_ph.dat', action='write', form='formatted' )
        write( un, '("#",a5,a26,a26)' ) 'iR', '|Rp|', 'max(|g(0,Rp)|)'
        ire = 1
        do irp = 1, eph_ph_mfi%nr
          irec = (irp - 1) * eph_el_mfi%nr + ire
          call gFMRR_file%read( irec, gFM_aW_RR(:, :, :, irp) )
          write( un, '(i6,2g26.16)' ) irp, eph_ph_mfi%rlen(irp), maxval( abs( gFM_aW_RR(:, :, :, irp) ) )
        end do
        close( un )
      end if

      ! close file for g_aW(Re,Rp)
      call gFMRR_file%close( mpiglobal )

      deallocate( gFM_aW_Rq, gFM_aW_RR )
    else
      call terminate_if_false( .false., '(eph_ephmat_setup_interpolation) &
        Electron-phonon matrix elements on coarse grids not found.' )
    end if
  end subroutine eph_ephmat_setup_interpolation
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! INTERPOLATE EPH MATRIX ELEMENTS
  !
  !> Get electron-phonon matrix elements \(g_{mn,\nu}({\bf k}', {\bf q}')\) for a given
  !> set of wave vectors \({\bf k}'\) and \({\bf q}'\) by Wannier-Fourier interpolation.
  !>
  !> This is done by Fourier interpolating the EPH matrix to \({\bf k}'\) and \({\bf q}'\) by
  !> \[ \mathcal{g}({\bf k}', {\bf q}') = \sum_{{\bf R}_e, {\bf R}_{ph}} {\rm e}^{{\rm i} 
  !> ({\bf k}'\cdot{\bf R}_e + {\bf q}'\cdot{\bf R}_{ph})} \mathcal{g}({\bf R}_e, {\bf R}_{ph}) \]
  !> and transforming it from atomic Wannier gauge to phonon Hamiltonian gauge by
  !> \[ g_{mn,\nu}({\bf k}',{\bf q}') = \sum\limits_{m',n',\kappa,\alpha} U_{m'm}({\bf k}'+{\bf q}')\, 
  !>    g_{m'n',\kappa,\alpha}({\bf k'}, {\bf q}')\, U_{n'n}^\ast({\bf k}')\, 
  !>    \frac{e_{\kappa\alpha,\nu}({\bf q}')}{\sqrt{2\,M_\kappa\,\omega_{\nu{\bf q}'})}} \; .\]
  !> See also [[transform_R2p(subroutine)]].
  !>
  !> MPI parallelization is over the product of \({\bf k}'\) and \({\bf R}_{ph}\) points for the Fourier
  !> transform \({\bf R}_e \rightarrow {\bf k}'\)
  !> and over the product of \({\bf k}'\) and \({\bf q}'\) points for the Fourier transform
  !> \({\bf R}_{ph} \rightarrow {\bf q}'\), 
  !> but each process can specify a different band and mode ranges `irange`, `frange`, and `mrange`.
  subroutine eph_ephmat_interpolate( vkl, vql, Umnk, Umnkq, ph_energy_q, ph_evec_q, ephmat, &
      electron_gauge, phonon_gauge, include_polar, irange, frange, mrange, mpicomm )
    use eph_electrons, only: eph_el_mfi
    use eph_phonons, only: eph_ph_mfi
    use constants, only: zone
    use exciting_mpi, only: xmpi_allreduce
    use mod_atoms, only: natmtot
    use modinput
#ifdef MPI
    use mpi_f08, only: MPI_Send, MPI_Recv, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, MPI_Comm
#endif
    !> set of wave vectors \({\bf k}'\) in lattice coordinates
    real(dp), intent(in) :: vkl(:,:)
    !> set of wave vectors \({\bf q}'\) in lattice coordinates
    real(dp), intent(in) :: vql(:,:)
    !> Wannier transformation matrices \(U_{mn}({\bf k}')\)
    complex(dp), intent(in) :: Umnk(:,:,:)
    !> Wannier transformation matrices \(U_{mn}({\bf k}'+{\bf q}')\)
    complex(dp), intent(in) :: Umnkq(:,:,:,:)
    !> phonon energies \(\omega_{\nu{\bf q}'}\)
    real(dp), intent(in) :: ph_energy_q(:,:)
    !> phonon eigenvectors \(e_{\kappa\alpha,\nu}({\bf q}')\)
    complex(dp), intent(in) :: ph_evec_q(:,:,:)
    !> EPH matrix \(g_{mn,\nu}({\bf k}', {\bf q}')\)
    complex(dp), allocatable, intent(out) :: ephmat(:,:,:,:,:)
    !> electron gauge (`'H'` for Hamiltonian (default), `'W'` for Wannier)
    character, optional, intent(in) :: electron_gauge
    !> phonon gauge (`'p'` for phonon (default), `'a'` for atomic)
    character, optional, intent(in) :: phonon_gauge
    !> include long-range part in polar materials (default: `.true.`)
    logical, optional, intent(in) :: include_polar
    !> range of band index for initial state \(|n{\bf k}'\rangle\) (default: all rows of \(U({\bf k}')\))
    integer, optional, intent(in) :: irange(2)
    !> range of band index for final state \(|m{\bf k}'+{\bf q}'\rangle\) (default: all rows of \(U({\bf k}'+{\bf q}')\))
    integer, optional, intent(in) :: frange(2)
    !> range of phonon modes \(\nu\) (default: all columns of \(e({\bf q}')\))
    integer, optional, intent(in) :: mrange(2)
    !> MPI communicator (default: global MPI communicator)
    type(mpiinfo), optional, intent(inout) :: mpicomm

    integer :: nk, nq, nwfk, nwfkq, nmode, nelem, irng(2), frng(2), mrng(2)
    integer :: ik, ik1, ik2, iq, iq1, iq2, ire, irp, ir1, ir2, irec, irank, jrank
    logical :: polar, hgauge, pgauge
    type(mpiinfo) :: mpi
    type(block_data_file_type) :: gFMRR_file

    integer, allocatable :: ik_range(:,:), iq_range(:,:), ir_range(:,:), irng_list(:,:), frng_list(:,:), mrng_list(:,:)
    complex(dp), allocatable :: g_aW_RR(:,:,:,:), g_aW_kR(:,:,:,:,:), g_aW_kq(:,:,:,:,:), &
                                g_lr_pref_a(:,:,:), g_xH(:,:,:), g_pX(:,:,:)

    ! set number of interpolation points
    nk = size( vkl, dim=2 )
    nq = size( vql, dim=2 )
    
    ! check input
    call assert( size( vkl, dim=1 ) == 3, &
      '`vkl` must be a set of vectors of length 3.' )
    call assert( size( vql, dim=1 ) == 3, &
      '`vql` must be a set of vectors of length 3.' )
    call assert( size( Umnk, dim=3 ) == nk, &
      '3rd dimension of `Umnk` must equal number of k-vectors.' )
    call assert( size( Umnkq, dim=3 ) == nk, &
      '3rd dimension of `Umnkq` must equal number of k-vectors.' )
    call assert( size( Umnkq, dim=4 ) == nq, &
      '4th dimension of `Umnkq` must equal number of q-vectors.' )
    call assert( size( ph_energy_q, dim=2 ) == nq, &
      '2nd dimension of `ph_energy_q` must equal number of q-vectors.' )
    call assert( size( ph_evec_q, dim=3 ) == nq, &
      '3rd dimension of `ph_evec_q` must equal number of q-vectors.' )
    call assert( size( Umnk, dim=2 ) == eph_nwf_tot, &
      '2nd dimension of `Umnk` must equal total number of Wannier functions.' )
    call assert( size( Umnk, dim=2 ) == eph_nwf_tot, &
      '2nd dimension of `Umnkq` must equal total number of Wannier functions.' )
    call assert( size( ph_energy_q, dim=1 ) == size( ph_evec_q, dim=2 ), &
      '1st dimension of `ph_energy_q` must equal 2nd dimension of `ph_evec_q`.' )
    call assert( size( ph_evec_q, dim=1 ) == 3*natmtot, &
      '1st dimension of `ph_evec_q` must equal `3*natmtot`.' )
    if (present( electron_gauge )) &
      call assert( any( electron_gauge == ['H', 'h', 'W', 'w'] ), &
        '`electron_gauge` must be either `H` (Hamiltonian gauge) or `W` (Wannier gauge).' )
    if (present( phonon_gauge )) &
      call assert( any( phonon_gauge == ['p', 'P', 'a', 'A'] ), &
        '`phonon_gauge` must be either `p` (phonon gauge) or `a` (atomic gauge).' )

    ! set defaults
    polar = .true.
    if (present( include_polar )) polar = include_polar
    hgauge = .true.
    if (present( electron_gauge )) hgauge = .not. any( electron_gauge == ['W', 'w'] )
    pgauge = .true.
    if (present( phonon_gauge )) pgauge = .not. any( phonon_gauge == ['a', 'A'] )
    mpi = mpiglobal
    if (present(mpicomm)) mpi = mpicomm

    ! set matrix size
    nwfk = size( Umnk, dim=1 )
    if (.not. hgauge) nwfk = eph_nwf_tot
    nwfkq = size( Umnkq, dim=1 )
    if (.not. hgauge) nwfkq = eph_nwf_tot
    nmode = size( ph_energy_q, dim=1 )
    if (.not. pgauge) nmode = 3*natmtot

    ! set band and mode ranges
    irng = [1, nwfk]
    if (present(irange)) irng = irange
    frng = [1, nwfkq]
    if (present(frange)) frng = frange
    mrng = [1, nmode]
    if (present(mrange)) mrng = mrange

    ! set total number of matrix elements
    nelem = eph_nwf_tot**2 * 3*natmtot

    ! communicate information on band and mode distribution
    allocate( irng_list(2, mpi%procs), frng_list(2, mpi%procs), mrng_list(2, mpi%procs) )
    irng_list = 0; irng_list(:, mpi%rank+1) = irng
    frng_list = 0; frng_list(:, mpi%rank+1) = frng
    mrng_list = 0; mrng_list(:, mpi%rank+1) = mrng
    call xmpi_allreduce( irng_list, mpi )
    call xmpi_allreduce( frng_list, mpi )
    call xmpi_allreduce( mrng_list, mpi )

    ! open file with g_aW(Re,Rp)
    gFMRR_file = block_data_file_type( eph_ephmat_gFMRR_filename, [eph_nwf_tot, eph_nwf_tot, 3*natmtot], cmplx( 0, 0, dp ) )
    call gFMRR_file%open( mpi )

    !................................................................................ 
    ! Fourier transform Re -> k'
    !
    ! distribute k and Rp among processes
    call patchwork_distribution( eph_ph_mfi%nr, nk, mpi%procs, ir_range, ik_range )
    ir1 = ir_range(1, mpi%rank+1)
    ir2 = ir_range(2, mpi%rank+1)
    ik1 = ik_range(1, mpi%rank+1)
    ik2 = ik_range(2, mpi%rank+1)
    allocate( g_aW_RR(eph_nwf_tot, eph_nwf_tot, 3*natmtot, eph_el_mfi%nr) )
    if (mpi%rank == 0) then
      allocate( g_aW_kR(eph_nwf_tot, eph_nwf_tot, 3*natmtot, nk, eph_ph_mfi%nr) )
    else
      allocate( g_aW_kR(eph_nwf_tot, eph_nwf_tot, 3*natmtot, ik1:ik2, ir1:ir2) )
    end if

    do irp = ir1, ir2
      ! read g_aW(Re,Rp) from file
      do ire = 1, eph_el_mfi%nr
        irec = (irp - 1) * eph_el_mfi%nr + ire
        call gFMRR_file%read( irec, g_aW_RR(:, :, :, ire) )
      end do
      ! transform to g_aW(k,Rp)
      if (ik_range(3, mpi%rank+1) > 0) then
        call eph_el_mfi%transform_R2p( [eph_nwf_tot, eph_nwf_tot], 3*natmtot, &
          g_aW_RR, eph_nwf_tot**2, 3*natmtot, &
          g_aW_kR(1, 1, 1, ik1, irp), eph_nwf_tot**2, 3*natmtot, vkl(:, ik1:ik2), &
          minimal_distances=.true. )
      end if
    end do
    deallocate( g_aW_RR )

    ! collect g_aW_kR at rank 0
#ifdef MPI
    if (mpi%rank == 0) then
      do irank = 1, mpi%procs-1
        do irp = ir_range(1, irank+1), ir_range(2, irank+1)
          call MPI_Recv( g_aW_kR(1, 1, 1, ik_range(1, irank+1), irp), nelem*ik_range(3, irank+1), MPI_DOUBLE_COMPLEX, irank, irank*eph_ph_mfi%nr+irp, MPI_Comm(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
        end do
      end do
    else
      do irp = ir1, ir2
        call MPI_Send( g_aW_kR(1, 1, 1, ik1, irp), nelem*ik_range(3, mpi%rank+1), MPI_DOUBLE_COMPLEX, 0, mpi%rank*eph_ph_mfi%nr+irp, MPI_Comm(mpi%comm), mpi%ierr )
      end do
      deallocate( g_aW_kR )
    end if
#endif
    !................................................................................ 

    !................................................................................ 
    ! Fourier transform Rp -> q'
    !
    ! distribute k and q among processes
    call patchwork_distribution( nk, nq, mpi%procs, ik_range, iq_range )
    ik1 = ik_range(1, mpi%rank+1)
    ik2 = ik_range(2, mpi%rank+1)
    iq1 = iq_range(1, mpi%rank+1)
    iq2 = iq_range(2, mpi%rank+1)
    if (mpi%rank /= 0) allocate( g_aW_kR(eph_nwf_tot, eph_nwf_tot, 3*natmtot, ik1:ik2, eph_ph_mfi%nr) )
    allocate( g_aW_kq(eph_nwf_tot, eph_nwf_tot, 3*natmtot, ik1:ik2, iq1:iq2) )

    ! distribute g_aW_kR among processes
#ifdef MPI
    if (mpi%rank == 0) then
      do irank = 1, mpi%procs-1
        if (ik_range(3, irank+1) < 1) cycle
        do irp = 1, eph_ph_mfi%nr
          call MPI_Send( g_aW_kR(1, 1, 1, ik_range(1, irank+1), irp), nelem*ik_range(3, irank+1), MPI_DOUBLE_COMPLEX, irank, irank*eph_ph_mfi%nr+irp, MPI_Comm(mpi%comm), mpi%ierr )
        end do
      end do
    else if (ik_range(3, mpi%rank+1) > 0) then
      do irp = 1, eph_ph_mfi%nr
        call MPI_Recv( g_aW_kR(1, 1, 1, ik1, irp), nelem*ik_range(3, mpi%rank+1), MPI_DOUBLE_COMPLEX, 0, mpi%rank*eph_ph_mfi%nr+irp, MPI_Comm(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
      end do
    end if
#endif

    ! transform to g_aW_(k,q)
    if (iq_range(3, mpi%rank+1) > 0) then
      call eph_ph_mfi%transform_R2p( [eph_nwf_tot**2, 3*natmtot], ik_range(3, mpi%rank+1), &
        g_aW_kR(1, 1, 1, ik1, 1), nelem, size( g_aW_kR, dim=4 ), &
        g_aW_kq, nelem, size( g_aW_kq, dim=4 ), vql(:, iq1:iq2), &
        minimal_distances=.false. )
    end if

    deallocate( g_aW_kR )
    !................................................................................ 

    !................................................................................ 
    ! add long-range contribution and transform gauge
    !
    ! generate long range prefactors in atomic gauge
    if (eph_polar .and. polar) then
      allocate( g_lr_pref_a(3, natmtot, iq1:iq2) )
      !$omp parallel default( shared )
      !$omp do
      do iq = iq1, iq2
        if (input%eph%elphbolt) then
          call eph_ephmat_gen_lr_prefactor_atomic_elphbolt( eph_ph_mfi%bvec, vql(:, iq), eph_dielten, eph_borncharge, g_lr_pref_a(:, :, iq), eph_ph_mfi%ngrid, 14.0_dp )
        else
          call eph_ephmat_gen_lr_prefactor_atomic( eph_ph_mfi%bvec, vql(:, iq), eph_dielten, eph_borncharge, g_lr_pref_a(:, :, iq) )
        end if
      end do
      !$omp end do
      !$omp end parallel
    end if

    ! transform gauge
    allocate( g_xH(nwfkq, nwfk, 3*natmtot) )
    allocate( g_pX(eph_nwf_tot, eph_nwf_tot, nmode) )
    do iq = iq1, iq2
      do ik = ik1, ik2
        ! add long range matrix elements in atomic Wannier gauge
        if (eph_polar .and. polar) &
          call eph_ephmat_add_lr( zone, g_lr_pref_a(:, :, iq), g_aW_kq(:, :, :, ik, iq) )
        ! transform to Hamiltonian gauge
        if (hgauge) then
          call eph_ephmat_transform_Hamiltonian_Wannier( Umnk(:, :, ik), Umnkq(:, :, ik, iq), g_xH, g_aW_kq(:, :, :, ik, iq), 3*natmtot, -1 )
          g_aW_kq(:nwfkq, :nwfk, :, ik, iq) = g_xH
        end if
        ! transform to phonon gauge
        if (pgauge) then
          call eph_ephmat_transform_atomic_phonon( ph_energy_q(:, iq), ph_evec_q(:, :, iq), g_aW_kq(:, :, :, ik, iq), g_pX, eph_nwf_tot**2, 1 )
          g_aW_kq(:, :, :nmode, ik, iq) = g_pX
        end if
      end do
    end do
    deallocate( g_xH, g_pX )
    if (eph_polar .and. polar) deallocate( g_lr_pref_a )
    !................................................................................ 

    ! distribute final result
    if (allocated(ephmat)) deallocate( ephmat )
    allocate( ephmat(frng(1):frng(2), irng(1):irng(2), mrng(1):mrng(2), nk, nq) )
    nelem = size( ephmat(:, :, :, 1, 1) )
    ephmat(:, :, :, ik1:ik2, iq1:iq2) = g_aW_kq(frng(1):frng(2), irng(1):irng(2), mrng(1):mrng(2), ik1:ik2, iq1:iq2)
#ifdef MPI
    do jrank = 0, mpi%procs-1
      if (mpi%rank == jrank) then
        do irank = 0, mpi%procs-1
          if (irank == jrank) cycle
          do iq = iq_range(1, irank+1), iq_range(2, irank+1)
            do ik = ik_range(1, irank+1), ik_range(2, irank+1)
              call MPI_Recv( ephmat(:, :, :, ik, iq), nelem, MPI_DOUBLE_COMPLEX, irank, irank*nq*nk+iq*nk+ik, MPI_Comm(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
            end do
          end do
        end do
      else
        do iq = iq1, iq2
          do ik = ik1, ik2
            g_xH = g_aW_kq(frng_list(1, jrank+1):frng_list(2, jrank+1), irng_list(1, jrank+1):irng_list(2, jrank+1), mrng_list(1, jrank+1):mrng_list(2, jrank+1), ik, iq)
            call MPI_Send( g_xH, size( g_xH ), MPI_DOUBLE_COMPLEX, jrank, mpi%rank*nq*nk+iq*nk+ik, MPI_Comm(mpi%comm), mpi%ierr )
          end do
        end do
      end if
    end do
#endif
    deallocate( g_aW_kq )

    ! close file with g_aW(Re,Rp)
    call gFMRR_file%close( mpi )
  end subroutine eph_ephmat_interpolate
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! READ MATRIX ELEMENTS ON COARSE GRIDS FROM FILE
  !
  !> Read electron-phonon matrix elements in atomic Hamiltonian gauge from file.
  !>
  !> The matrix elements at \({\bf k}\) and \({\bf q}\) are obtained from the 
  !> matrix elements at the symmetry equivalent point \({\bf q}_0\), which is related
  !> to \({\bf q}\) via the symmetry \(\mathcal{S}=\lbrace \mathrm{\bf R} | {\bf \tau} \rbrace\) 
  !> such that \({\bf q} = \mathrm{\bf R}{\bf q}_0 + {\bf G}\).
  !> The matrix elements are rotated using the phonon symmetry matrix \(\Gamma(\mathcal{S};{\bf q})\)
  !> (see [[ph_util_symmetry_G(subroutine)]])
  !> \[ g_{mn,\kappa\alpha}({\bf k},{\bf q}) = \sum_{\lambda,\beta} 
  !>    g_{mn,\lambda\beta}(\mathrm{\bf R}^{-1}{\bf k},{\bf q}_0)\, 
  !>    \Gamma_{\kappa\alpha,\lambda\beta}^\ast(\mathcal{S};{\bf q}) \;. \]
  subroutine eph_ephmat_read_coarse( kset, qset, vkl, vql, file, g_aH_kq )
    use constants, only: zzero, zone
    use dfpt_variables, only: dfpt_kset, dfpt_Gset, dfpt_Gkset, feveck0
    use dfpt_eigensystem, only: dfpt_eig_getevec
    use phonons_util, only: ph_util_symmetry_G
    use mod_kpointset, only: k_set, Gk_set, generate_k_vectors, generate_Gk_vectors, delete_k_vectors, delete_Gk_vectors
    use mod_atoms, only: natmtot
    use mod_symmetry, only: lsplsymc, symlat
    use m_linalg, only: zlsp
    !> set of \({\bf k}_0\)-vectors the matrix elements have been calculated for
    type(k_set), intent(in) :: kset
    !> set of \({\bf q}_0\)-vectors the matrix elements have been calculated for
    type(k_set), intent(in) :: qset
    !> \({\bf k}\)-vector (in lattice Coordinates)
    real(dp), intent(in) :: vkl(3)
    !> \({\bf q}\)-vector (in lattice Coordinates)
    real(dp), intent(in) :: vql(3)
    !> file that contains matrix elements at \({\bf k}_0\) and \({\bf q}_0\)
    type(block_data_file_type), intent(inout) :: file
    !> matrix elements \(g_{mn,\kappa\alpha}({\bf k},{\bf q})\)
    complex(dp), intent(out) :: g_aH_kq(:,:,:)

    integer :: iq0, ik0, ik, imode, isymq, isymk, irec, ilspl, f_rank
    real(dp) :: vql0(3), vkl0(3)
    type(k_set) :: kqset
    type(Gk_set) :: Gkqset

    integer, allocatable :: f_shape(:)
    real(dp), allocatable :: vgkql(:,:)
    complex(dp), allocatable :: sym_G(:,:), g_aH_k0q0(:,:,:), Uk(:,:), Ukq(:,:), evec0(:,:), evec(:,:)

    ! check input
    f_rank = file%get_block_rank()
    f_shape = file%get_block_shape()
    call assert( f_rank == 3, &
      'Data blocks in file must have rank 3.' )
    call assert( all( f_shape == shape( g_aH_kq ) ), &
      'Shape of data blocks in file and shape of matrix elements are not compatible.' )
    call assert( f_shape(3) == 3*natmtot, &
      'Size of dimension 3 must be `3*natmtot`.' )

    ! generate G+k vectors if not yet done
    if (.not. allocated(eph_Gkset_el%ngk)) &
      call generate_Gk_vectors( eph_Gkset_el, eph_kset_el, dfpt_Gset, dfpt_Gkset%gkmax )

    ! find equivalent q-point q0 and connecting symmetry
    call findkptinset( vql, qset, isymq, iq0 )
    vql0 = qset%vkl(:, iq0)
    ! find rotated k-point k0=S^-1.k
    ilspl = lsplsymc(isymq)
    call r3mtv( dble( symlat(:, :, ilspl) ), vkl, vkl0 )
    call findkptinset( vkl0, kset, isymk, ik0 )
    call terminate_if_false( isymk == 1, '(eph_ephmat_read_coarse) &
      Rotated k-point not included in k-set.' )

    ! read matrix elements at k0 and q0
    allocate( g_aH_k0q0(f_shape(1), f_shape(2), f_shape(3)) )
    irec = (iq0 - 1) * kset%nkpt + ik0
    call file%read( irec, g_aH_k0q0 )

    !................................................................................ 
    ! phonon rotation
    !
    ! get symmetry matrix Gamma
    sym_G = ph_util_symmetry_G( isymq, qset%vkl(:, iq0) )
    ! apply phonon rotation to matrix elements
    call zgemm( 'n', 'c', f_shape(1)*f_shape(2), 3*natmtot, 3*natmtot, zone, &
      g_aH_k0q0, f_shape(1)*f_shape(2), &
      sym_G, 3*natmtot, zzero, &
      g_aH_kq, f_shape(1)*f_shape(2) )
    !................................................................................ 

    !................................................................................ 
    ! electron rotation
    !
    if (isymq /= 1) then
      ! get unitary rotation between wavefunctions at k and k0
      allocate( Uk(eph_nst, eph_nst) )
      call findkptinset( vkl0, kset, isymk, ik0 )
      call dfpt_eig_getevec( kset%vkl(:, ik0), eph_Gkset_el%vgkl(:, :, 1, ik0), feveck0, dfpt_kset, dfpt_Gkset, [eph_fst, eph_lst], evec0 )

      call findkptinset( vkl, kset, isymk, ik )
      call dfpt_eig_getevec( kset%vkl(:, ik), eph_Gkset_el%vgkl(:, :, 1, ik), feveck0, dfpt_kset, dfpt_Gkset, [eph_fst, eph_lst], evec )

      call rotate_evecfv( isymq, kset%vkl(:, ik0), kset%vkl(:, ik), &
             eph_Gkset_el%ngk(1, ik0), eph_Gkset_el%vgkl(:, :, 1, ik0), eph_Gkset_el%vgkl(:, :, 1, ik), &
             evec0, size( evec0, dim=1 ), size( evec0, dim=2 ) )
      call zlsp( evec0, evec, Uk )
      ! get unitary rotation between wavefunctions at k+q and k0+q0
      allocate( Ukq(eph_nst, eph_nst) )
      call generate_k_vectors( kqset, kset%bvec, [1, 1, 1], [0.0_dp, 0.0_dp, 0.0_dp], .false., .false. )
      kqset%vkl(:, 1) = vkl0 + vql0 
      call r3mv( kqset%bvec, kqset%vkl(:, 1), kqset%vkc(:, 1) )
      call generate_Gk_vectors( Gkqset, kqset, dfpt_Gset, dfpt_Gkset%gkmax )
      call dfpt_eig_getevec( vkl0+vql0, Gkqset%vgkl(:, :, 1, 1), feveck0, dfpt_kset, dfpt_Gkset, [eph_fst, eph_lst], evec0 )
      vgkql = Gkqset%vgkl(:, :, 1, 1)
      call delete_k_vectors( kqset )
      call delete_Gk_vectors( Gkqset )

      call generate_k_vectors( kqset, kset%bvec, [1, 1, 1], [0.0_dp, 0.0_dp, 0.0_dp], .false., .false. )
      kqset%vkl(:, 1) = vkl + vql 
      call r3mv( kqset%bvec, kqset%vkl(:, 1), kqset%vkc(:, 1) )
      call generate_Gk_vectors( Gkqset, kqset, dfpt_Gset, dfpt_Gkset%gkmax )
      call dfpt_eig_getevec( vkl+vql, Gkqset%vgkl(:, :, 1, 1), feveck0, dfpt_kset, dfpt_Gkset, [eph_fst, eph_lst], evec )

      call rotate_evecfv( isymq, vkl0+vql0, vkl+vql, &
             Gkqset%ngk(1, 1), vgkql, Gkqset%vgkl(:, :, 1, 1), &
             evec0, size( evec0, dim=1 ), size( evec0, dim=2 ) )
      call zlsp( evec0, evec, Ukq )
      ! apply rotations
      call zgemm( 'c', 'n', eph_nst, f_shape(2)*3*natmtot, eph_nst, zone, &
        Ukq, eph_nst, &
        g_aH_kq, f_shape(1), zzero, &
        g_aH_k0q0, f_shape(1) )
      do imode = 1, 3*natmtot
        call zgemm( 'n', 'n', eph_nst, eph_nst, eph_nst, zone, &
          g_aH_k0q0(:, :, imode), f_shape(1), &
          Uk, eph_nst, zzero, &
          g_aH_kq(:, :, imode), f_shape(1) )
      end do

      ! free memory
      deallocate( Uk, Ukq, evec, evec0, vgkql )
      call delete_k_vectors( kqset )
      call delete_Gk_vectors( Gkqset )
    end if
    !................................................................................ 

    deallocate( f_shape, g_aH_k0q0, sym_G )
  end subroutine eph_ephmat_read_coarse
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! READ POTENTIAL RESPONSE
  !
  !> Read (soft) effective potential response \(\delta^{\bf q}_{\kappa \alpha} V({\bf r})\)
  !> for given \({\bf q}\) point, atom \(\kappa\) and polarization direction \(\alpha\) from file.
  subroutine eph_ephmat_read_dpot( qset, vql, lmaxvr, dpot_mt, Gqset, dpot_ir, &
      directory )
    use dfpt_inout, only: dfpt_io_read_zfun
    use phonons_io_util, only: ph_io_qsap_string
    use phonons_density_potential, only: ph_rhopot_rotate_q_canonical
    use mod_kpointset, only: k_set, G_set
    use mod_atoms, only: nspecies, natoms, idxas
    use mod_muffin_tin, only: nrmt
    use modinput
    !> set of \({\bf q}\)-vectors the potential response \(\delta^{\bf q}_{\kappa \alpha} V({\bf r})\) have been computed for
    type(k_set), intent(in) :: qset
    !> \({\bf q}\)-point (lattice coordinates) to read
    real(dp), intent(in) :: vql(3)
    !> maximum angular momentum \(l\) used in potential response expansion
    integer, intent(in) :: lmaxvr
    !> muffin-tin effective potential response
    complex(dp), intent(out) :: dpot_mt(:,:,:,:,:)
    !> set of \({\bf G}+{\bf q}\) vectors used in potential response expansion
    type(G_set), intent(in) :: Gqset
    !> interstitial effective potential response
    complex(dp), intent(out) :: dpot_ir(:,:,:)
    !> path to directory (default: current directory)
    character(*), optional, intent(in) :: directory

    integer :: is, ia, ias, ip, iq, isym, lmmaxvr
    character(:), allocatable :: fxt, errmsg

    lmmaxvr = (lmaxvr + 1)**2

    ! find equivalent q-point in set
    call findkptinset( vql, qset, isym, iq )

    ! read potential response at equivalent q-point
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        do ip = 1, 3
          fxt = ph_io_qsap_string( qset%ivk(:, iq), qset%ngridk, is, ia, ip )
          call read_potential_response( 'DVEFF_'//trim(fxt)//'.OUT', nspecies, natoms(1:nspecies), nrmt(1:nspecies), lmaxvr, Gqset%ngrid, &
            dpot_mt(:, :, :, ip, ias), [size(dpot_mt, dim=1), size(dpot_mt, dim=2), size(dpot_mt, dim=3)], &
            dpot_ir(:, ip, ias), [size(dpot_ir, dim=1)] )
          ! WARNING! We have to change the sign for Gamma point phonons from super-cell calculations
          ! due to a inconsistency in the super-cell code.
          if (input%phonons%method == 'sc' .and. all( qset%ivk(:, iq) == 0 )) then
            dpot_mt(:, :, :, ip, ias) = - dpot_mt(:, :, :, ip, ias)
            dpot_ir(:, ip, ias) = - dpot_ir(:, ip, ias)
          end if
        end do
      end do
    end do

    ! rotate potential response from equivalent point into requested point if symmetry is not the identity
    if( isym /= 1 ) &
      call ph_rhopot_rotate_q_canonical( qset%vkl(:, iq), vql, isym, lmaxvr, dpot_mt, Gqset, dpot_ir )
  end subroutine eph_ephmat_read_dpot
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! GAUGE INDEPENDENT AVERAGE OF MATRIX ELEMENTS
  !
  !> Compute gauge independent average of matrix elements.
  !>
  !> The gauge independent average of matrix elements is defined as
  !> \[ |g|_{mn,\nu}({\bf k},{\bf q}) = \left[
  !>    \sum_{m' \in \mathcal{D}_{{\bf k}+{\bf q}}} \sum_{n' \in \mathcal{D}_{\bf k}} \sum_{\nu' \in \mathcal{D}_{\bf q}} 
  !>    \frac{|g_{m'n',\nu'}({\bf k},{\bf q})|^2}
  !>         {|\mathcal{D}_{{\bf k}+{\bf q}}|\, |\mathcal{D}_{\bf k}|\, |\mathcal{D}_{\bf q}|}
  !>    \right]^{1/2} \; , \]
  !> where \(\mathcal{D}_{{\bf k}+{\bf q}}\) is the set of electronic states that are degenerate with energy
  !> \(\epsilon_{m{\bf k}+{\bf q}}\), \(\mathcal{D}_{\bf k}\) is the set of electronic states that 
  !> are degenerate with energy \(\epsilon_{n{\bf k}}\), and \(\mathcal{D}_{\bf q}\) is the set of phonon modes
  !> that are degenerate with energy \(\omega_{\nu{\bf q}}\).
  subroutine eph_ephmat_average( el_energy_k, el_energy_kq, eps_el, ph_energy_q, eps_ph, g_kq, g_kq_avg )
    use math_utils, only: get_degeneracies
    !> electron energies \(\epsilon_{n{\bf k}}\)
    real(dp), intent(in) :: el_energy_k(:)
    !> electron energies \(\epsilon_{n{\bf k}+{\bf q}}\)
    real(dp), intent(in) :: el_energy_kq(:)
    !> tolerance for degeneracy of electrons
    real(dp), intent(in) :: eps_el
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:)
    !> tolerance for degeneracy of phonons
    real(dp), intent(in) :: eps_ph
    !> matrix elements \(g_{mn,\nu}({\bf k},{\bf q})\)
    complex(dp), intent(in) :: g_kq(:,:,:)
    !> averaged matrix elements \(|g|_{mn,\nu}({\bf k},{\bf q})\)
    real(dp), allocatable, intent(out) :: g_kq_avg(:,:,:)

    integer :: ist, jst, imode
    integer, allocatable :: deg_k(:,:), deg_kq(:,:), deg_q(:,:)

    ! check input
    call assert( size( el_energy_kq ) == size( g_kq, dim=1 ), &
      'Number of electron energies at k+q and size of 1st dimension of matrix elements must be equal.' )
    call assert( size( el_energy_k ) == size( g_kq, dim=2 ), &
      'Number of electron energies at k and size of 2nd dimension of matrix elements must be equal.' )
    call assert( size( ph_energy_q ) == size( g_kq, dim=3 ), &
      'Number of phonon energies at q and size of 3rd dimension of matrix elements must be equal.' )

    allocate( g_kq_avg, source=abs( g_kq ) )
    deg_kq = get_degeneracies( el_energy_kq, eps_el )
    deg_k = get_degeneracies( el_energy_k, eps_el )
    deg_q = get_degeneracies( ph_energy_q, eps_ph )
    do imode = 1, size( deg_q, dim=2 )
      do jst = 1, size( deg_k, dim=2 )
        do ist = 1, size( deg_kq, dim=2 )
          g_kq_avg(deg_kq(1, ist):deg_kq(2, ist), deg_k(1, jst):deg_k(2, jst), deg_q(1, imode):deg_q(2, imode)) = sqrt( &
            sum( g_kq_avg(deg_kq(1, ist):deg_kq(2, ist), deg_k(1, jst):deg_k(2, jst), deg_q(1, imode):deg_q(2, imode))**2 ) / &
            (deg_kq(3, ist) * deg_k(3, jst) * deg_q(3, imode)) )
        end do
      end do
    end do
  end subroutine eph_ephmat_average
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! LONG-RANGE FRÖHLICH MATRIX ELEMENTS
  !
  !> Compute the \({\bf k}\) independent prefactor of the long-range Fröhlich matrix elements
  !> in atomic Wannier gauge.
  !>
  !> \[ 
  !>    \frac{4\pi}{\Omega} \sum\limits_{{\bf G+q} \neq 0} 
  !>    {\rm e}^{-{\rm i}{\bf G+q} \cdot {\bf \tau}_{\kappa}}
  !>    \frac{{\bf G+q}^\top \cdot {\bf Z}^\ast_{\kappa}}{{\bf G+q}^\top \cdot {\bf \epsilon}^\infty \cdot {\bf G+q}} 
  !>    {\rm e}^{-\frac{{\bf G+q}^\top \cdot {\bf \epsilon}^\infty \cdot {\bf G+q}}{4 \lambda^2}}
  !> \]
  subroutine eph_ephmat_gen_lr_prefactor_atomic( bvec, vql, dielten, borncharge, prefactor )
    use constants, only: zzero, twopi, fourpi
    use mod_atoms, only: natmtot, nspecies, natoms, atposc
    !> reciprocal lattice vectors (columnwise)
    real(dp), intent(in) :: bvec(3, 3)
    !> \({\bf q}\)-point in lattice coordinates
    real(dp), intent(in) :: vql(3)
    !> dielectric tensor \({\bf \epsilon}^\infty\)
    real(dp), intent(in) :: dielten(3, 3)
    !> Born-effective charges \({\bf Z}^\ast_\alpha\)
    real(dp), intent(in) :: borncharge(3, 3, natmtot)
    !> prefactors
    complex(dp), intent(out) :: prefactor(3, natmtot)

    real(dp), parameter :: tol = 1e-12_dp          ! tolerance for terms to include
    real(dp), parameter :: gmax = -log(tol)        ! exp(-gmax) ~ tol
    real(dp), parameter :: lambda = 2.0_dp

    integer :: i, is, ia, ias, ig1, ig2, ig3, ng(3)
    real(dp) :: t1, t2, vqlbz(3), vgql(3), gql, gqegq, beb(3, 3), bz(3, 3), bt(3, natmtot), dotp

    real(8), external :: r3mdet

    prefactor = zzero

    ! map q to 1BZ
    vqlbz = vql
    call r3ws( tol, bvec, vqlbz, ng )
    ! return if q = 0
    if (norm2( vqlbz ) < tol) return
    ! 4pi/Omega
    t1 = fourpi * abs( r3mdet( bvec ) ) / twopi**3
    ! B^T.e.B
    beb = matmul( transpose( bvec ), matmul( dielten, bvec ) )
    ! B^T.tau
    bt = reshape( [( ( matmul( transpose( bvec ), atposc(:, ia, is) ), ia=1, natoms(is) ), is=1, nspecies )], [3, natmtot] )
    ! get limits for reciprocal space sum
    ng = 3 * [( ceiling( sqrt( gmax * 4.0_dp * lambda**2 / beb(i, i) ) ), i=1, 3 )]

    ! reciprocal space sum
    do ig3 = -ng(3), ng(3)
      do ig2 = -ng(2), ng(2)
        do ig1 = -ng(1), ng(1)
          vgql = dble( [ig1, ig2, ig3] ) + vqlbz
          gql = norm2( vgql )
          vgql = vgql / gql
          gqegq = dot_product( vgql, matmul( beb, vgql ) )
          t2 = gql**2 * gqegq / (4.0_dp * lambda**2)
          if ((gqegq < tol) .or. (t2 > gmax)) cycle
          do ias = 1, natmtot
            dotp = - gql * dot_product( vgql, bt(:, ias) )
            prefactor(:, ias) = prefactor(:, ias) + &
              exp( -t2 ) * vgql / (gql * gqegq) * cmplx( cos( dotp ), sin( dotp ), dp )
          end do
        end do
      end do
    end do

    ! product with Born charges and prefactor
    do ias = 1, natmtot
      ! B^T.Z
      bz = matmul( transpose( bvec ), borncharge(:, :, ias) )
      prefactor(:, ias) = t1 * cmplx( matmul( prefactor(:, ias)%re, bz ), &
                                      matmul( prefactor(:, ias)%im, bz ), dp )
    end do
  end subroutine eph_ephmat_gen_lr_prefactor_atomic

  !> See [[eph_ephmat_gen_lr_prefactor_atomic(subroutine)]].
  !>
  !> This version works in Cartesian coordiantes and uses different numerical parameters than [[eph_ephmat_gen_lr_prefactor_atomic(subroutine)]]
  !> in order to be compatible with the codes `elphbolt` and `EPW`.
  subroutine eph_ephmat_gen_lr_prefactor_atomic_elphbolt( bvec, vql, dielten, borncharge, prefactor, g_grid, g_max )
    use constants, only: zzero, twopi, fourpi
    use mod_atoms, only: natmtot, nspecies, natoms, idxas, atposc
    !> reciprocal lattice vectors (columnwise)
    real(dp), intent(in) :: bvec(3, 3)
    !> \({\bf q}\)-point in lattice coordinates
    real(dp), intent(in) :: vql(3)
    !> dielectric tensor \({\bf \epsilon}^\infty\)
    real(dp), intent(in) :: dielten(3, 3)
    !> Born-effective charges \({\bf Z}^\ast_\alpha\)
    real(dp), intent(in) :: borncharge(3, 3, natmtot)
    !> prefactors
    complex(dp), intent(out) :: prefactor(3, natmtot)
    !> grid size for reciprocal space sum
    integer, intent(in) :: g_grid(3)
    !> cut-off for reciprocal space sum
    real(dp), intent(in) :: g_max

    real(dp), parameter :: tol = 1e-12_dp          ! tolerance for terms to include

    integer :: ng(3), ig1, ig2, ig3, is, ia, ias
    real(dp) :: lambda, gmax
    real(dp) :: avec(3, 3), vqlbz(3), vgqc(3), gqc, gqegq, t1, t2, dotp

    real(8), external :: r3mdet

    prefactor = zzero

    ! map q to 1BZ
    vqlbz = vql
    call r3ws( tol, bvec, vqlbz, ng )
    ! return if q = 0
    if (norm2( vqlbz ) < tol) return
    ! 4pi/Omega
    t1 = fourpi * abs( r3mdet( bvec ) ) / twopi**3
    ! A = 2pi B^-T
    call r3minv( bvec, avec )
    avec = twopi * transpose( avec )

    gmax = g_max
    ng = g_grid
    lambda = twopi / norm2( avec(:, 1) )

    ! reciprocal space sum
    do ig3 = -ng(3), ng(3)
      do ig2 = -ng(2), ng(2)
        do ig1 = -ng(1), ng(1)
          vgqc = matmul( bvec, dble( [ig1, ig2, ig3] ) + vqlbz )
          gqc = norm2( vgqc )
          vgqc = vgqc / gqc
          gqegq = dot_product( vgqc, matmul( dielten, vgqc ) )
          t2 = gqc**2 * gqegq / (4.0_dp * lambda**2)
          if ((gqegq < tol) .or. (t2 > gmax)) cycle
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)
              dotp = - gqc * dot_product( vgqc, atposc(:, ia, is) )
              prefactor(:, ias) = prefactor(:, ias) + &
                exp( -t2 ) * vgqc / (gqc * gqegq) * &
                cmplx( cos( dotp ), sin( dotp ), dp )
            end do
          end do
        end do
      end do
    end do

    ! product with Born charges and prefactor
    do ias = 1, natmtot
      prefactor(:, ias) = cmplx( matmul( prefactor(:, ias)%re, borncharge(:, :, ias) ),&
                                 matmul( prefactor(:, ias)%im, borncharge(:, :, ias) ), dp )
    end do
    prefactor = t1 * prefactor
  end subroutine eph_ephmat_gen_lr_prefactor_atomic_elphbolt

  !> Add or subtract long-range Fröhlich matrix elements.
  !>
  !> Computes \(\mathcal{\bf g} \leftarrow \mathcal{\bf g} + \alpha \mathcal{\bf g}^\mathcal{L}\).
  !> This subroutine works in Wannier gauge in which the long-range term is diagonal.
  subroutine eph_ephmat_add_lr( alpha, prefactor, g_aW )
    use mod_atoms, only: natmtot
    !> prefactor \(\alpha\)
    complex(dp), intent(in) :: alpha
    !> long-range matrix element prefactors in atomic coordinates
    complex(dp), intent(in) :: prefactor(3, natmtot)
    !> e-ph matrix elements in atomic Wannier gauge
    complex(dp), intent(inout) :: g_aW(:,:,:)

    integer :: nwfk, nwfkq, ist, ip, ias
    complex(dp) :: ialpha

    call assert( size( g_aW, dim=3 ) == 3*natmtot, &
      '3rd dimension of `g_aW` must equal `3*natmtot`.' )
    
    nwfk = size( g_aW, dim=2 )
    nwfkq = size( g_aW, dim=1 )

    ! i*alpha
    ialpha = cmplx( -aimag( alpha ), dble( alpha ), dp )

    do ias = 1, natmtot
      do ip = 1, 3
        do ist = 1, min( nwfk, nwfkq )
          g_aW(ist, ist, (ias-1)*3+ip) = g_aW(ist, ist, (ias-1)*3+ip) + ialpha * prefactor(ip, ias)
        end do
      end do
    end do
  end subroutine eph_ephmat_add_lr
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! TRANSFORM BETWEEN HAMILTONIAN AND WANNIER GAUGE
  !
  !> Transform electron-phonon matrix elements between Hamitlonian and Wannier gauge.
  !>
  !> Given the Wannier transformation matrices \({\bf U}^({\bf k})\) and \({\bf U}({\bf k+q})\),
  !> the transformation from Hamiltonian to Wannier gauge is given by
  !> \[ \mathcal{\bf g} = {\bf U}({\bf k+q})^\dagger \, {\bf g} \, {\bf U}({\bf k}) \]
  !> and the transformation from Wannier to Hamiltonian gauge is given by
  !> \[ {\bf g} = {\bf U}({\bf k+q}) \, \mathcal{\bf g} \, {\bf U}({\bf k})^\dagger \]
  subroutine eph_ephmat_transform_Hamiltonian_Wannier( Umnk, Umnkq, g_Hamilton, g_Wannier, nmat, dir )
    use constants, only: zzero, zone
    !> Wannier transformation matrix \(U_{mn}({\bf k})\)
    complex(dp), intent(in) :: Umnk(:,:)
    !> Wannier transformation matrix \(U_{mn}({\bf k}+{\bf q})\)
    complex(dp), intent(in) :: Umnkq(:,:)
    !> matrix elements \(g_{mn,:}({\bf k},{\bf q})\) in Hamitlonian gauge
    complex(dp), intent(inout) :: g_Hamilton(:,:,:)
    !> matrix elements \(\mathcal{g}_{mn,:}({\bf k},{\bf q})\) in Wannier gauge
    complex(dp), intent(inout) :: g_Wannier(:,:,:)
    !> number of matrices \(g\)
    integer, intent(in) :: nmat
    !> direction (`>0`: Hamiltonian to Wannier, `<0`: Wannier to Hamiltonian)
    integer, intent(in) :: dir

    integer :: i
    
    complex(dp), allocatable :: aux(:,:)

    if (dir == 0) return

    ! check input
    call assert( size( g_Hamilton, dim=1 ) == size( Umnkq, dim=1 ), &
      'Size of 1st dimension of `g_Hamilton` must equal size of 1st dimension of `Umnkq`.' )
    call assert( size( g_Hamilton, dim=2 ) == size( Umnk, dim=1 ), &
      'Size of 2nd dimension of `g_Hamilton` must equal size of 1st dimension of `Umnk`.' )
    call assert( size( g_Wannier, dim=1 ) == size( Umnkq, dim=2 ), &
      'Size of 1st dimension of `g_Wannier` must equal size of 2nd dimension of `Umnkq`.' )
    call assert( size( g_Wannier, dim=2 ) == size( Umnk, dim=2 ), &
      'Size of 2nd dimension of `g_Wannier` must equal size of 2nd dimension of `Umnk`.' )

    if (dir > 0) then
      allocate( aux(size( g_Hamilton, dim=1 ), size( Umnk, dim=2 )) )
      do i = 1, nmat
        call zgemm( 'n', 'n', size( g_Hamilton, dim=1 ), size( Umnk, dim=2 ), size( Umnk, dim=1 ), zone, &
          g_Hamilton(:, :, i), size( g_Hamilton, dim=1 ), &
          Umnk, size( Umnk, dim=1 ), zzero, &
          aux, size( aux, dim=1 ) )
        call zgemm( 'c', 'n', size( Umnkq, dim=2 ), size( Umnk, dim=2 ), size( Umnkq, dim=1 ), zone, &
          Umnkq, size( Umnkq, dim=1 ), &
          aux, size( aux, dim=1 ), zzero, &
          g_Wannier(:, :, i), size( g_Wannier, dim=1 ) )
      end do
      deallocate( aux )
    else
      allocate( aux(size( Umnkq, dim=1 ), size( g_Wannier, dim=2 )) )
      do i = 1, nmat
        call zgemm( 'n', 'n', size( Umnkq, dim=1 ), size( g_Wannier, dim=2 ), size( Umnkq, dim=2 ), zone, &
          Umnkq, size( Umnkq, dim=1 ), &
          g_Wannier(:, :, i), size( g_Wannier, dim=1 ), zzero, &
          aux, size( aux, dim=1 ) )
        call zgemm( 'n', 'c', size( Umnkq, dim=1), size( Umnk, dim=1 ), size( Umnk, dim=2 ), zone, &
          aux, size( aux, dim=1 ), &
          Umnk, size( Umnk, dim=1 ), zzero, &
          g_Hamilton(:, :, i), size( g_Hamilton, dim=1 ) )
      end do
      deallocate( aux )
    end if
  end subroutine eph_ephmat_transform_Hamiltonian_Wannier
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! TRANSFORM BETWEEN ATOMIC AND PHONON GAUGE
  !
  !> Transform electron-phonon matrix elements between atomic and phonon gauge.
  !>
  !> Given the phonon frequencies \(\omega_{\nu}({\bf q})\) and eigenvectors \({\bf e}_{\nu}({\bf q})\), 
  !> the transformation from atomic to phonon gauge is given by
  !> \[ {\bf g}^{\rm ph}_{\nu}({\bf k},{\bf q}) = \sum\limits_{\kappa, \alpha} \frac{1}{\sqrt{2 \, \omega_{\nu}({\bf q}) \, M_{\kappa}}} {\bf g}^{\rm at}_{\kappa\alpha}({\bf k},{\bf q}) \, e_{\kappa\alpha,\nu}({\bf q}) \]
  !> and the transformation from phonon to atomic gauge is given by
  !> \[ {\bf g}^{\rm at}_{\kappa\alpha}({\bf k},{\bf q}) = \sum\limits_{\nu} \sqrt{2 \, \omega_{\nu}({\bf q}) \, M_{\kappa}} {\bf g}^{\rm ph}_{\nu}({\bf k},{\bf q}) \, e_{\kappa\alpha,\nu}^\ast({\bf q}) \]
  subroutine eph_ephmat_transform_atomic_phonon( ph_energy_q, ph_evec_q, g_atomic, g_phonon, ld, dir )
    use constants, only: zzero, zone
    use mod_atoms, only: natmtot, nspecies, natoms, idxas, spmass
    !> phonon energies \(\omega_{\nu{\bf q}}\)
    real(dp), intent(in) :: ph_energy_q(:)
    !> phonon eigenvectors \(e_{\kappa\alpha,\nu}({\bf q})\)
    complex(dp), intent(in) :: ph_evec_q(:,:)
    !> leading dimension
    integer, intent(in) :: ld
    !> object in atomic coordinates \(g^{\rm at}_{:,\kappa\alpha}\)
    complex(dp), intent(inout) :: g_atomic(ld, *)
    !> object in phonon coordinates \(g^{\rm ph}_{:,\nu}\)
    complex(dp), intent(inout) :: g_phonon(ld, *)
    !> direction (`>0`: atomic to phonon, `<0`: phonon to atomic)
    integer, intent(in) :: dir

    integer :: imode, is, ia, ias
    real(dp) :: t1

    complex(dp), allocatable :: transform(:,:)

    if (dir == 0) return

    ! check input
    call assert( size( ph_energy_q ) == size( ph_evec_q, dim=2 ), &
      'Different number of modes for phonon energies and eigenvectors.' )
    call assert( size( ph_evec_q, dim=1 ) == 3*natmtot, &
      '1st dimension of phonon eigenvectors must be `3*natmtot`.' )

    allocate( transform(3*natmtot, size( ph_energy_q )) )
    if (dir > 0) then
      do imode = 1, size( ph_energy_q )
        if (ph_energy_q(imode) > eph_ph_energy_zero) then
          t1 = 1.0_dp / sqrt( 2.0_dp * ph_energy_q(imode) )
        else
          t1 = 0.0_dp
        end if
        transform(:, imode) = ph_evec_q(:, imode) * t1
      end do
      do is = 1, nspecies
        if (spmass(is) > 1e-12_dp) then
          t1 = 1.0_dp / sqrt( spmass(is) )
        else
          t1 = 0.0_dp
        end if
        do ia = 1, natoms(is)
          ias = (idxas(ia, is) - 1) * 3 + 1
          transform(ias:ias+2, :) = transform(ias:ias+2, :) * t1
        end do
      end do
      call zgemm( 'n', 'n', ld, size( ph_energy_q ), 3*natmtot, zone, &
        g_atomic, ld, &
        transform, 3*natmtot, zzero, &
        g_phonon, ld )
    else
      do imode = 1, size( ph_energy_q )
        if (ph_energy_q(imode) > eph_ph_energy_zero) then
          t1 = sqrt( 2.0_dp * ph_energy_q(imode) )
        else
          t1 = 0.0_dp
        end if
        transform(:, imode) = ph_evec_q(:, imode) * t1
      end do
      do is = 1, nspecies
        if (spmass(is) > 1e-12_dp) then
          t1 = sqrt( spmass(is) )
        else
          t1 = 0.0_dp
        end if
        do ia = 1, natoms(is)
          ias = (idxas(ia, is) - 1) * 3 + 1
          transform(ias:ias+2, :) = transform(ias:ias+2, :) * t1
        end do
      end do
      call zgemm( 'n', 'c', ld, 3*natmtot, size( ph_energy_q ), zone, &
        g_phonon, ld, &
        transform, 3*natmtot, zzero, &
        g_atomic, ld )
    end if
    deallocate( transform )
  end subroutine eph_ephmat_transform_atomic_phonon
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! AUXILIARY PROCEDURES
  !
  !> Generate radial integrals times Gaunt coefficients for potential gradient using 
  !> regularized effective potential.
  subroutine eph_ephmat_gen_grad_pot_mt_basis( lmaxvr, gpot_mt_basis )
    use constants, only: zzero, zone, y00
    use matrix_elements, only: me_mt_alloc, me_mt_prepare
    use mod_potential_and_density, only : pot_mt => veffmt, potxc_mt => vxcmt, potcl_mt => vclmt
    use mod_atoms, only: nspecies, natoms, idxas, spzn, spr
    use mod_muffin_tin, only : nrmtmax, nrmt
    !> maximum angular momentum \(l\) used in potential expansion
    integer, intent(in) :: lmaxvr
    !> radial integrals times Gaunt coefficients for potential gradient
    complex(dp), allocatable, intent(out) :: gpot_mt_basis(:,:,:,:)

    integer :: lmmaxvr, is, ia, ias, ip, ir
    complex(dp), allocatable :: potr_mt(:,:), gpotr_mt(:,:,:)

    lmmaxvr = (lmaxvr + 1)**2

    call me_mt_alloc( gpot_mt_basis, 3 )

    ! contributions from ionic potential gradients
    call rignt_ion( gpot_mt_basis )

    ! contributions from regularized potential gradients
    allocate( potr_mt(lmmaxvr, nrmtmax) )
    allocate( gpotr_mt(lmmaxvr, nrmtmax, 3) )
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        ! get regularized muffin-tin potential (without diverging ionic potential)
        ! and its gradients
        do ir = 1, nrmt(is)
          call rtozflm( lmaxvr, pot_mt(:, ir, ias), potr_mt(:, ir) )
          potr_mt(1, ir) = (potcl_mt(1, ir, ias) * spr(ir, is) - spzn(is) / y00) / spr(ir, is) + potxc_mt(1, ir, ias)
        end do
        call gradzfmt( lmaxvr, nrmt(is), spr(1:nrmt(is), is), lmmaxvr, nrmtmax, potr_mt, gpotr_mt )

        do ip = 1, 3
          call me_mt_prepare( is, ias, lmaxvr, -zone, gpotr_mt(:, :, ip), zone, gpot_mt_basis(:, :, ias, ip) )
        end do
      end do
    end do
    deallocate( potr_mt, gpotr_mt )

  contains
    ! Compute radial integrals times Gaunt coefficients for the gradient of the ionic potential
    ! grad Z/|r| = -Z/|r|^2 * (r/|r|) 
    ! where the normal vector r/|r| is expressed with spherical harmonics.
    subroutine rignt_ion( gion )
      use constants, only: fourpi
      use gaunt
      use matrix_elements_lapw_lo, only: basis
      complex(dp), intent(out) :: gion(:,:,:,:)

      integer :: l1, l2, m1, m2, lm, lm1, lm2, n1, n2, idx1, idx2, lam1, lam2, i
      complex(dp) :: vgion(4, 3)
      type(non_zero_gaunt_real), pointer :: gntr

      real(dp), allocatable :: f1(:), f2(:), f(:), df(:), cf(:,:)
      complex(dp), allocatable :: zri_gion(:,:,:)

      vgion = zzero
      do ip = 1, 3
        if (ip == 1) then
          vgion(2, ip) = -sqrt(1.0_dp / 6.0_dp)
          vgion(4, ip) =  sqrt(1.0_dp / 6.0_dp)
        else if (ip == 2) then
          vgion(2, ip) = -sqrt(1.0_dp / 6.0_dp) * cmplx(0, 1, dp)
          vgion(4, ip) = -sqrt(1.0_dp / 6.0_dp) * cmplx(0, 1, dp)
        else if (ip == 3) then
          vgion(3, ip) = -sqrt(1.0_dp / 3.0_dp)
        end if
      end do
      vgion = vgion * sqrt(fourpi)

      allocate( zri_gion(basis%n_rad_fun_max, basis%n_rad_fun_max, 4) )
      gntr => gaunt_coeff_yyy

      do is = 1, nspecies
        allocate( f(basis%n_rad_grid(is)), df(basis%n_rad_grid(is)), cf(basis%n_rad_grid(is), 3) )
        do ia = 1, natoms(is)
          ias = idxas(ia, is)

          do ip = 1, 3
            do l1 = 0, basis%lmax_basis
              n1 = basis%n_rad_fun(l1, is)
              do l2 = 0, basis%lmax_basis
                n2 = basis%n_rad_fun(l2, is)
                zri_gion = cmplx( 0, 0, dp )
                do lam1 = 1, n1
                  f1 = basis%get_rad_fun( l1, is, ias, lam1 )
                  do lam2 = 1, n2
                    f2 = basis%get_rad_fun( l2, is, ias, lam2 )
                    f = f1 * f2
                    call fderiv( -1, basis%n_rad_grid(is), basis%rad_grid(:, is), f, df, cf )
                    zri_gion(lam1, lam2, :) = df(basis%n_rad_grid(is)) * vgion(:, ip)
                  end do
                end do
                
                do m1 = -l1, l1
                  lm1 = l1 * (l1 + 1) + m1 + 1
                  do m2 = -l2, l2
                    lm2 = l2 * (l2 + 1) + m2 + 1
                    do lam1 = 1, basis%n_rad_fun(l1, is)
                      idx1 = basis%idx_basis_fun(lm1, lam1, is)
                      do lam2 = 1, basis%n_rad_fun(l2, is)
                        idx2 = basis%idx_basis_fun(lm2, lam2, is)

                        do i = 1, gntr%num(lm1, lm2)
                          lm = gntr%lm2(i, lm1, lm2)
                          if (lm > 4) exit
                          gion(idx1, idx2, ias, ip) = gion(idx1, idx2, ias, ip) - &
                            spzn(is) * gntr%val(i, lm1, lm2) * zri_gion(lam1, lam2, lm) 
                        end do

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

        end do
        deallocate( f, df, cf )
      end do
    end subroutine rignt_ion
  end subroutine eph_ephmat_gen_grad_pot_mt_basis

  !> For a given effective potential response \(\delta V_{\rm eff}({\bf r})\), 
  !> this subroutine computes the radial muffin-tin integrals times Gaunt coefficients
  !> and the reciprocal space representation of the effective potential response
  !> times the characteristic function.
  subroutine eph_ephmat_prepare_gmat( lmaxvr, dpot_mt, Gset, dpot_ir, dpot_mt_basis, dpot_cfun_ig, &
      alpha, beta )
    use matrix_elements, only: me_mt_prepare, me_ir_prepare
    use mod_kpointset, only: G_set
    use mod_atoms, only: nspecies, natoms, idxas
    !> maximum angular momentum \(l\) used in potential response expansion
    integer, intent(in) :: lmaxvr
    !> muffin-tin effective potential response
    complex(dp), intent(in) :: dpot_mt(:,:,:)
    !> set of \({\bf G}\) vectors used in potential response expansion
    type(G_set), intent(in) :: Gset
    !> interstitial effective potential response
    complex(dp), intent(in) :: dpot_ir(:)
    !> radial muffin-tin integrals times Gaunt coefficients
    complex(dp), intent(inout) :: dpot_mt_basis(:,:,:)
    !> interstitial effective potential response times characteristic function in reciprocal space
    complex(dp), intent(inout) :: dpot_cfun_ig(:)
    !> prefactor of input (default: 1)
    complex(dp), optional, intent(in) :: alpha
    !> prefactor of output (default: 0)
    complex(dp), optional, intent(in) :: beta

    integer :: is, ia, ias
    complex(dp) :: a, b

    a = cmplx( 1, 0, dp )
    if (present(alpha)) a = alpha
    b = cmplx( 0, 0, dp )
    if (present(beta)) b = beta

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call me_mt_prepare( is, ias, lmaxvr, a, dpot_mt(:, :, ias), b, dpot_mt_basis(:, :, ias) )
      end do
    end do
    
    call me_ir_prepare( a, dpot_ir, b, dpot_cfun_ig, Gset_op=Gset )
  end subroutine eph_ephmat_prepare_gmat

  !> Having prepared the \({\bf k}\)-independent part of the e-ph matrix elements, i.e.,
  !> the radial integrals times the Gaunt coefficients and the interstitial potential response
  !> times the characteristic function, this computes the total e-ph matrix elements with the
  !> given eigenvectors.
  subroutine eph_ephmat_gen_gmat( ik, Gkset, Gkqset, Gqset, fst1, lst1, fst2, lst2, dpot_mt_basis, dpot_cfun_ig, evec1, evec2, &
      gmat )
    use constants, only: zzero, zone
    use matrix_elements, only: me_mt_mat, me_ir_mat
    use mod_kpointset, only: k_set, G_set, Gk_set
    use mod_eigensystem, only: nmatmax_ptr
    use mod_gkvector, only: ngkmax_ptr
    use mod_APW_LO, only: nlotot, apwordmax
    use mod_muffin_tin, only: lmmaxapw
    use mod_atoms, only: natmtot, nspecies, natoms, idxas
    !> index of \({\bf k}\)- and \({\bf k}+{\bf q}\)-point
    integer, intent(in) :: ik
    !> set of \({\bf G}+{\bf k}\) vectors
    type(Gk_set), intent(in) :: Gkset
    !> set of \({\bf G}+{\bf k}+{\bf q}\) vectors
    type(Gk_set), intent(in) :: Gkqset
    !> set of \({\bf G}+{\bf q}\) vectors used in expansion of interstitial potential response
    type(G_set), intent(in) :: Gqset
    !> first and last state on the left for which the matrix elements are calculated
    integer, intent(in) :: fst1, lst1
    !> first and last state on the right for which the matrix elements are calculated
    integer, intent(in) :: fst2, lst2
    !> radial muffin-tin integrals of effective potential response times Gaunt coefficients
    complex(dp), intent(in) :: dpot_mt_basis(:,:,:)
    !> interstitial effective potential response times characteristic function in reciprocal space
    complex(dp), intent(in) :: dpot_cfun_ig(:)
    !> eigenvectors on the left
    complex(dp), intent(in) :: evec1(:, fst1:)
    !> eigenvectors on the right
    complex(dp), intent(in) :: evec2(:, fst2:)
    !> e-ph matrix elements
    complex(dp), intent(out) :: gmat(fst1:, fst2:)

    integer :: is, ia, ias
    integer, target :: ngkmax, ngkqmax, nmatmax

    complex(dp), allocatable :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)

    gmat = zzero

    ngkmax = Gkset%ngkmax
    ngkqmax = Gkqset%ngkmax
    nmatmax = ngkmax + nlotot
    nmatmax_ptr => nmatmax

    ! get matching coefficients at k
    ngkmax_ptr => ngkmax
    allocate( apwalmk(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
    call match( Gkset%ngk(1, ik), Gkset%gkc(:, 1, ik), Gkset%tpgkc(:, :, 1, ik), Gkset%sfacgk(:, :, 1, ik), &
      apwalmk )
    ! get matching coefficients at k+q
    ngkmax_ptr => ngkqmax
    allocate( apwalmkq(ngkmax_ptr, apwordmax, lmmaxapw, natmtot) )
    call match( Gkqset%ngk(1, ik), Gkqset%gkc(:, 1, ik), Gkqset%tpgkc(:, :, 1, ik), Gkqset%sfacgk(:, :, 1, ik), &
      apwalmkq )

    ! compute matrix elements
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call me_mt_mat( is, ias, Gkqset%ngk(1, ik), Gkset%ngk(1, ik), apwalmkq(:, :, :, ias), apwalmk(:, :, :, ias), &
          evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
          zone, dpot_mt_basis(:, :, ias), zone, gmat )
      end do
    end do
    call me_ir_mat( Gkqset, ik, Gkset, ik, &
      evec1(:, fst1:lst1), evec2(:, fst2:lst2), &
      zone, dpot_cfun_ig, zone, gmat, &
      Gset_op=Gqset )

    ! free memory
    deallocate( apwalmk, apwalmkq )
  end subroutine eph_ephmat_gen_gmat
  !-------------------------------------------------------------------------------- 
end module eph_ephmat
