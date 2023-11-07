!> This is the main module for calculating phonons using 
!> density-functional perturbation theory (DFPT).
module phonons
  use dfpt_variables
  use dfpt_inout
  use phonons_variables
  use phonons_inout

  use modmpi
  use precision, only: dp
  use block_data_file, only: block_data_file_type

  implicit none
  private

  ! LOCAL VARIABLES
  !> `.true.` if current \({\bf q}\) point is the Gamma point
  logical :: gamma
  !> eigenvalues at all \({\bf k}\) and \({\bf k+q}\) points
  real(dp), allocatable :: evalk(:,:), evalkq(:,:)
  !> occupation numbers at all \({\bf k}\) and \({\bf k+q}\) points
  real(dp), allocatable :: occk(:,:), occkq(:,:)
  !> eigenvectors at single \({\bf k}\) and \({\bf k+q}\) point
  complex(dp), allocatable :: eveck(:,:), eveckq(:,:)
  !> (L)APW matching coefficients \(A^\kappa_{{\bf G+p},lm,\xi}\) at \({\bf k}\) and \({\bf k+q}\)
  complex(dp), allocatable :: apwalmk(:,:,:,:), apwalmkq(:,:,:,:)
  !> eigenvalue and occupation response at single \({\bf k}\) point
  real(dp), allocatable :: devalk(:,:,:), docck(:,:,:)
  !> eigenvector response at single \({\bf k}\) point
  complex(dp), allocatable :: deveck(:,:)
  !> overlap matrix response at all \({\bf k}\) points and irrep members
  complex(dp), allocatable :: dSmat(:,:,:,:)
  !> constant part of Hamiltonian matrix response at all \({\bf k}\) points and irrep members
  complex(dp), allocatable :: dHmat_const(:,:,:,:)
  !> full Hamiltonian matrix response at single \({\bf k}\) point and irrep member
  complex(dp), allocatable :: dHmat(:,:)
  !> density response matrix for all atoms and irrep members
  complex(dp), allocatable :: drho_mat(:,:,:,:)
  !> muffin-tin and interstitial density response for all irrep members
  complex(dp), allocatable :: drho_mt(:,:,:,:), drho_ir(:,:)
  !> muffin-tin and interstitial effective potential response for all irrep members
  complex(dp), allocatable :: dpot_mt(:,:,:,:), dpot_ir(:,:)
  !> radial integrals of effective potential response times Gaunt coefficients
  complex(dp), allocatable :: dHmat_mt_basis(:,:,:,:)
  !> interstitial potential response times characteristic function in reciprocal space
  complex(dp), allocatable :: dpot_cfun_ig(:,:)
  !> interstitial (scalar relativistic) kinetic energy response times characteristic function in reciprocal space
  complex(dp), allocatable :: dkin_cfun_ig(:,:)
  !> constant part of force response for all canonical displacements
  complex(dp), allocatable :: dforce_const(:,:,:)
  !> force response for all atoms and irrep members
  complex(dp), allocatable :: dforce(:,:,:)
  !> polarization response for all irrep members
  complex(dp), allocatable :: dpol(:,:)
  !> parallel i/o file objects for eigenvectors 
  !> at \({\bf q}\)-dependent \({\bf k}\)- and \({\bf k+q}\)-grid
  type(block_data_file_type) :: feveck, feveckq
  !> parallel i/o file objects for eigenvalue, occupation and eigenvector response 
  !> at \({\bf q}\)-dependent \({\bf k}\)-grid
  type(block_data_file_type) :: fdevalk, fdocck, fdeveck
  !> MPI communicator for all processes working on a part and on the first \({\bf k}\)-point within this part
  type(mpiinfo) :: mpilocal, mpilocalk
  !> `.true.`, if this process is the local master of the current part
  logical :: master
  !> global rank of the local master
  integer :: localmasterrank

  public :: ph_prepare, ph_finalize, ph_dry_run
  public :: ph_part_prepare, ph_part_finalize, ph_part_scf, ph_part_force, ph_part_polarization
  public :: ph_write_dyn_canonical, ph_write_dpot_canonical, ph_write_borncharge_canonical
  public :: ph_dynmat_canonical_from_file, ph_dpot_canonical_from_file, ph_borncharge_canonical_from_file

  contains

    !> Perform a phonons dry run.
    !>
    !> Set up irreps and parallelization settings for a hypothetical parallel run
    !> with `input%phonons%drynumprocs` processes. The information and settings
    !> for the calculation will be written to `PHONON_RUN_INFO.OUT`.
    subroutine ph_dry_run
      use phonons_parallelization, only: ph_par_distribute
      use modinput

      ! initialize global phonon variables
      call ph_var_init
      ph_numprocs = input%phonons%drynumprocs
      ! find independent parts and their computational load
      ! distribute parts according to their load among processes
      call ph_io_find_done_parts( ph_qset, ph_irrep_basis(1:), ph_parts_done )
      call ph_par_distribute( ph_qset, ph_irrep_basis(1:), ph_numprocs, &
             input%phonons%minprocsperpart, input%phonons%maxprocsperpart, ph_parts_done, &
             ph_parts, ph_schedule, all_parts=ph_parts_all )
      ! print run information
      call ph_io_print_calc_info( print_schedule=input%phonons%write_schedule )
      ! free memory
      call ph_var_free
    end subroutine ph_dry_run

    !> This subroutine executes preparative tasks for a DFPT phonon calculation.
    !>
    !> This includes:
    !>
    !> * initialization of global phonon variables
    !> * distribution of independent calcualtion parts among MPI processes
    !> * obtaining eigenvalues and occupation numbers on the \({\bf k}_0\) points
    !>   of the unperturbed system
    !> * calculation of the constant parts for density and potential response
    subroutine ph_prepare( &
        do_force )
      use dfpt_eigensystem, only: dfpt_eig_ks
      use phonons_parallelization, only: ph_par_distribute, ph_par_parts_from_schedule
      use phonons_density_potential, only: ph_rhopot_init
      use phonons_force, only: ph_frc_dpulay_k, ph_frc_dsurf_k, ph_frc_symmetrize

      use constants, only: zzero
      use mod_APW_LO, only: nlotot, apwordmax
      use mod_muffin_tin, only: lmmaxapw
      use mod_atoms, only: natmtot
      use mod_eigenvalue_occupancy, only: nstfv, occmax, efermi
      use mod_charge_and_moment, only: chgval
      use mod_Gkvector, only: ngkmax_ptr
      use modinput
      !> prepare the constant part of the force response (default: `.false.`)
      logical, optional, intent(in) :: do_force

      integer :: ik, ik1, ik2, nmatmaxk, ip
      integer, target :: ngkmaxk
      logical :: force, exists, success

      force = .false.
      if( present( do_force ) ) force = do_force

      ! find independent parts and their computational load
      ! distribute parts according to their load among processes
      call ph_io_find_done_parts( ph_qset, ph_irrep_basis(1:), ph_parts_done )
      call ph_par_distribute( ph_qset, ph_irrep_basis(1:), ph_numprocs, &
             input%phonons%minprocsperpart, input%phonons%maxprocsperpart, ph_parts_done, &
             ph_parts, ph_schedule, all_parts=ph_parts_all )
      ph_parts_per_rank = ph_par_parts_from_schedule( ph_schedule(mpiglobal%rank+1, :) )

      ! print run information
      call ph_io_print_calc_info( print_schedule=input%phonons%write_schedule )

      ! initialize more global variables
      call ph_var_init_q(0)
      ! initialize constant parts for density and potential response
      call ph_rhopot_init
      ! try to read constant part of force response from file
      allocate( dforce_const(3, natmtot, 3), source=zzero )
      call ph_io_read_dforce_const( dforce_const, success )

      ! check if temporary file for eigenvectors and eigenvalues are available
      exists = fevalk0%exists() .and. feveck0%exists()
      exists = exists .and. (fevalk0%is_open() .and. feveck0%is_open())
      call terminate_if_false( exists, &
             '(ph_prepare): Eigenvalues and eigenvectors were not prepared. Call `dfpt_prepare` in advance.' )

      ! compute constant part of force response
      if( force .and. size( ph_parts ) > 0 .and. .not. success ) then
        ! set limits for k-point loops
        ik1 = firstofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
        ik2 = lastofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
        ! allocate local variables
        ngkmaxk = dfpt_Gkset%ngkmax
        nmatmaxk = ngkmaxk + nlotot
        allocate( evalk(nmatmaxk, dfpt_kset%nkpt), source=0.0_dp )
        allocate( occk(nstfv, dfpt_kset%nkpt) )
        allocate( eveck(nmatmaxk, nmatmaxk) )
        allocate( apwalmk(ngkmaxk, apwordmax, lmmaxapw, natmtot) )
        ! read eigenvalues from file
        do ik = ik1, ik2
          call dfpt_eig_ks( ik, dfpt_kset, dfpt_Gset, dfpt_Gkset, nmatmaxk, evalk(:,ik), eveck, &
                 p0set=dfpt_kset, Gp0set=dfpt_Gkset, feval=fevalk0, fevec=feveck0 )
        end do
        ! gather eigenvalues at all processes
        call mpi_allgatherv_ifc( dfpt_kset%nkpt, rlen=nmatmaxk, inplace=.true., comm=mpiglobal, rbuf=evalk )
        ! get occupation numbers
        occk = 0.0_dp
        call find_fermi( dfpt_kset%nkpt, dfpt_kset%wkpt, nstfv, evalk(1:nstfv, :), chgval, occmax, &
               input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
               efermi, occk )
        do ik = ik1, ik2
          ! get matching coefficients at k
          ngkmax_ptr => ngkmaxk
          call match( dfpt_Gkset%ngk(1, ik), dfpt_Gkset%gkc(:, 1, ik), dfpt_Gkset%tpgkc(:, :, 1, ik), dfpt_Gkset%sfacgk(:, :, 1, ik), &
                      apwalmk )
          ! read eigenvectors
          call feveck0%read( ik, eveck )
          ! add k-point contribution to force response
          do ip = 1, 3
            call ph_frc_dpulay_k( ik, dfpt_kset, dfpt_Gkset, dfpt_Gkset, 1, nstfv, &
                   evalk(:, ik), evalk(:, ik), occk(:, ik), occk(:, ik), eveck, eveck, apwalmk, apwalmk, &
                   ph_irrep_basis(0)%irreps(1)%pat(:, :, ip), .false., &
                   dforce_const(:, :, ip), order=1 )
            call ph_frc_dsurf_k( ik, dfpt_kset, dfpt_Gkset, dfpt_Gkset, 1, nstfv, &
                   evalk(:, ik), evalk(:, ik), occk(:, ik), occk(:, ik), eveck, eveck, &
                   ph_irrep_basis(0)%irreps(1)%pat(:, :, ip), .false., &
                   dforce_const(:, :, ip), order=1 )
          end do
        end do
        ! reduce force response over k-points
#ifdef MPI
        call MPI_Allreduce( MPI_IN_PLACE, dforce_const, size( dforce_const ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
               mpiglobal%comm, mpiglobal%ierr )
#endif
        ! symmetrize force response
        call ph_frc_symmetrize( dforce_const, [0.0_dp, 0.0_dp, 0.0_dp], 3, &
               ph_irrep_basis(0)%nsym, ph_irrep_basis(0)%isym, ph_irrep_basis(0)%ivsym, &
               ph_irrep_basis(0)%irreps(1)%symmat )

        ! write dforce_const to file
        if( mpiglobal%rank == 0 ) &
          call ph_io_write_dforce_const( dforce_const, success )

        deallocate( evalk, occk, eveck, apwalmk )
      end if

      ! free unneeded variables
      call ph_var_free_q

      call barrier( mpicom=mpiglobal )
    end subroutine ph_prepare

    !> This subroutine executes finalizing tasks for a DFPT phonon calculation.
    !>
    !> This includes freeing memory from unneeded variables.
    subroutine ph_finalize
      use phonons_density_potential, only: ph_rhopot_free
      if( allocated( dforce_const ) ) deallocate( dforce_const )
      call ph_rhopot_free
      call ph_var_free
    end subroutine ph_finalize

    !> This subroutine executes preparative tasks for an independent part
    !> of a DFPT phonon calculation.
    !>
    !> This includes:
    !> 
    !> * creating and entering a new subdirectory for this part
    !> * allocating module variables
    !> * obtaining eigenvalues, occupation numbers and eigenvectors of 
    !>   \({\bf k}\) and \({\bf k+q}\) for \({\bf k} \in \text{IBZ}({\bf q})\)
    !>   (eigenvalues and eigenvectors are written to temporary files
    !>   that are deleted at the end of the calculation)
    !> * calculation of the overlap matrix response and the part of the 
    !>   Hamiltonian matrix response that remain constant in the 
    !>   self-consistency cycle
    !> * initializing the density and potential response
    subroutine ph_part_prepare( ipart, &
        info_output )
      use dfpt_eigensystem, only: dfpt_eig_ks
      use phonons_eigensystem, only: ph_eig_gen_dSHmat
      use phonons_density_potential, only: ph_rhopot_gen_dpot, ph_rhopot_init_drho

      use constants, only: zzero, zone
      use mod_APW_LO, only: nlotot, apwordmax
      use mod_atoms, only: natmtot
      use mod_muffin_tin, only: lmmaxapw, nrmtmax
      use mod_eigenvalue_occupancy, only: nstfv, occmax, efermi
      use mod_charge_and_moment, only: chgval
      use mod_Gkvector, only: ngkmax_ptr
      use mod_misc, only: scrpath
      use modinput
      !> index of the independent part
      integer, intent(in) :: ipart
      !> create and write info output files (default: `.true.`)
      logical, optional, intent(in) :: info_output

      integer :: iq, iirrep, dirrep, &
                 id, id1, id2, ik, ik1, ik2, &
                 iv(3), nmatmaxk, nmatmaxkq
      integer, target :: ngkmaxk, ngkmaxkq
      real(dp) :: vql(3)
      logical :: write_info, success, from_file
      character(:), allocatable :: fxt

      write_info = .true.
      if( present( info_output ) ) write_info = info_output

      ! return if this rank is not working on this part
      if( .not. ph_parts(ipart)%is_my_rank( mpiglobal%rank ) ) return
      ! get parallelization information
      mpilocal = ph_parts(ipart)%mpi      ! all processes working on this task
      mpilocalk = ph_parts(ipart)%mpik    ! all processes working on the first k-point in this task
      master = (mpilocalk%rank == 0)      ! local master
      localmasterrank = 0                 ! rank of local master within local communicator
      if( master ) localmasterrank = mpilocal%rank
#ifdef MPI
      call MPI_Allreduce( MPI_IN_PLACE, localmasterrank, 1, MPI_INTEGER, MPI_SUM, mpilocal%comm, mpilocal%ierr )
#endif
      call barrier( mpicom=mpilocal )

      ! set q-point and irrep
      iq = ph_parts(ipart)%iq
      iirrep = ph_parts(ipart)%iirrep
      dirrep = ph_parts(ipart)%dirrep

      ! set limits for irrep member and k-point loops
      id1 = ph_parts(ipart)%first_d_of_rank( mpiglobal%rank )
      id2 = ph_parts(ipart)%last_d_of_rank( mpiglobal%rank )
      ik1 = ph_parts(ipart)%first_k_of_rank( mpiglobal%rank )
      ik2 = ph_parts(ipart)%last_k_of_rank( mpiglobal%rank )
      ! check for Gamma point
      vql = ph_qset%vkl(:, iq)
      call r3frac( input%structure%epslat, vql, iv )
      gamma = (norm2( vql ) < input%structure%epslat)
      ! generate subdirectory for part
      call ph_io_genqidir( iq, iirrep, ph_io_qi_dirname, comm=mpilocal )
      scrpath = trim( ph_io_qi_dirname )
      ! initialize info output
      if( master .and. write_info ) call ph_io_info_init( iq, iirrep, directory=ph_io_qi_dirname )
      ! write symmetries and parallelization settings to info file
      if( master .and. write_info ) call ph_io_info_symmetries_and_parallelization( ph_qset, ph_irrep_basis(1:), ph_parts, ipart, ph_parts_all )
      ! write displacement patterns to info file
      if( master .and. write_info ) call ph_io_info_irrep_patterns( ph_irrep_basis(iq)%irreps(iirrep) )
      ! initialize q-dependent variables
      call ph_var_init_q( iq )
      ngkmaxk = ph_Gkset%ngkmax
      ngkmaxkq = ph_Gkqset%ngkmax
      nmatmaxk = ngkmaxk + nlotot
      nmatmaxkq = ngkmaxkq + nlotot
      ! open temporary files for eigenvectors at k and k+q
      feveck = block_data_file_type( ph_io_qi_dirname//'EVECK.TMP', [nmatmaxk, nstfv], zone )
      feveckq = block_data_file_type( ph_io_qi_dirname//'EVECKQ.TMP', [nmatmaxkq, nmatmaxkq], zone )
      call feveck%open( mpilocal, delete_existing=.true. )
      call feveckq%open( mpilocal, delete_existing=.true. )
      ! allocate local variables
      allocate( evalk(nstfv, ph_kset%nkpt), evalkq(nmatmaxkq, ph_kqset%nkpt), devalk(nstfv, ph_kset%nkpt, dirrep) )
      allocate( occk(nstfv, ph_kset%nkpt), occkq(nmatmaxkq, ph_kqset%nkpt), docck(nstfv, ph_kset%nkpt, dirrep) )
      allocate( eveck(nmatmaxk, nstfv), eveckq(nmatmaxkq, nmatmaxkq), deveck(nmatmaxkq, nstfv) )
      allocate( apwalmk(ngkmaxk, apwordmax, lmmaxapw, natmtot), apwalmkq(ngkmaxkq, apwordmax, lmmaxapw, natmtot) )
      allocate( dSmat(nmatmaxkq, nstfv, ik1:ik2, id1:id2) )
      allocate( dHmat_const(nmatmaxkq, nstfv, ik1:ik2, id1:id2) )
      allocate( dHmat(nmatmaxkq, nstfv) )
      allocate( drho_mat(mt_basis%n_basis_fun_max, mt_basis%n_basis_fun_max, natmtot, id1:id2) )
      allocate( dHmat_mt_basis(mt_basis%n_basis_fun_max, mt_basis%n_basis_fun_max, natmtot, id1:id2) )
      allocate( dpot_cfun_ig(ph_Gqset%ngvec, id1:id2) )
      allocate( dkin_cfun_ig(ph_Gqset%ngvec, id1:id2) )
      allocate( dpol(3, dirrep) )
      if( master ) then ! only the master needs to know all irrep members
        allocate( drho_mt(dfpt_lmmaxvr, nrmtmax, natmtot, dirrep) )
        allocate( drho_ir(dfpt_Gset%ngrtot, dirrep) )
        allocate( dpot_mt(dfpt_lmmaxvr, nrmtmax, natmtot, dirrep) )
        allocate( dpot_ir(dfpt_Gset%ngrtot, dirrep) )
        allocate( dforce(3, natmtot, dirrep) )
      else ! the other processes only know their portions
        allocate( drho_mt(dfpt_lmmaxvr, nrmtmax, natmtot, id1:id2) )
        allocate( drho_ir(dfpt_Gset%ngrtot, id1:id2) )
        allocate( dpot_mt(dfpt_lmmaxvr, nrmtmax, natmtot, id1:id2) )
        allocate( dpot_ir(dfpt_Gset%ngrtot, id1:id2) )
        allocate( dforce(3, natmtot, id1:id2) )
      end if
      ! loop over k-points
      do ik = firstofset( mpilocal%rank, ph_kset%nkpt, mpilocal%procs ), lastofset( mpilocal%rank, ph_kset%nkpt, mpilocal%procs )
        ! solve KS equation at k
        call dfpt_eig_ks( ik, ph_kset, dfpt_Gset, ph_Gkset, nstfv, evalk(:, ik), eveck, &
               p0set=dfpt_kset, Gp0set=dfpt_Gkset, feval=fevalk0, fevec=feveck0 )
        ! solve KS equation at k+q
        call dfpt_eig_ks( ik, ph_kqset, dfpt_Gset, ph_Gkqset, nmatmaxkq, evalkq(:, ik), eveckq, &
               p0set=dfpt_kset, Gp0set=dfpt_Gkset, feval=fevalk0, fevec=feveck0 )
        ! write eigenvectors and eigenvalues to file
        call feveck%write( ik, eveck )
        call feveckq%write( ik, eveckq )
      end do
      call mpi_allgatherv_ifc( ph_kset%nkpt, rlen=nstfv, comm=mpilocal, rbuf=evalk )
      call mpi_allgatherv_ifc( ph_kset%nkpt, rlen=nmatmaxkq, comm=mpilocal, rbuf=evalkq )
      ! get occupation numbers
      occk = 0.0_dp
      call find_fermi( ph_kset%nkpt, ph_kset%wkpt, nstfv, evalk, chgval, occmax, &
             input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
             efermi, occk )
      occkq = 0.0_dp
      call find_fermi( ph_kqset%nkpt, ph_kqset%wkpt, nstfv, evalkq(1:nstfv,:), chgval, occmax, &
             input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
             efermi, occkq(1:nstfv,:) )
      ! compute overlap response and constant part of Hamiltonian response
      ! for all k points and irrep members
      dSmat = zzero; dHmat_const = zzero
      do ik = ik1, ik2
        ! get matching coefficients at k
        ngkmax_ptr => ngkmaxk
        call match( ph_Gkset%ngk(1, ik), ph_Gkset%gkc(:, 1, ik), ph_Gkset%tpgkc(:, :, 1, ik), ph_Gkset%sfacgk(:, :, 1, ik), &
                    apwalmk )
        ! get matching coefficients at k+q
        ngkmax_ptr => ngkmaxkq
        call match( ph_Gkqset%ngk(1, ik), ph_Gkqset%gkc(:, 1, ik), ph_Gkqset%tpgkc(:, :, 1, ik), ph_Gkqset%sfacgk(:, :, 1, ik), &
                    apwalmkq )
        ! read eigenvectors
        call feveck%read( ik, eveck )
        call feveckq%read( ik, eveckq )
        do id = id1, id2
          call ph_eig_gen_dSHmat( ik, ph_Gkset, ph_Gkqset, 1, nstfv, &
                 eveck, eveckq, apwalmk, apwalmkq, dSmat(:, :, ik, id), dHmat_const(:, :, ik, id), &
                 pat=ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id) )
        end do
      end do
      ! initialize density and potential response
      from_file = (input%phonons%do == 'fromfile')
      if( master ) then
        docck = 0.0_dp
        id = 0
        do while( id < dirrep )
          id = id + 1
          if( from_file ) then
            call ph_io_irrep_fxt( iq, iirrep, id, fxt )
            call dfpt_io_read_zfun( drho_mt(:, :, :, id), drho_ir(:, id), 'PHONON_DRHO', success, &
              file_extension=fxt, directory=ph_io_qi_dirname )
            from_file = from_file .and. success
            call dfpt_io_read_zfun( dpot_mt(:, :, :, id), dpot_ir(:, id), 'PHONON_DVEFF', success, &
              file_extension=fxt, directory=ph_io_qi_dirname )
            from_file = from_file .and. success
            if( .not. from_file ) id = 0
          else
            call ph_rhopot_init_drho( ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), &
                   drho_mt(:, :, :, id), drho_ir(:, id) )
            call ph_rhopot_gen_dpot( ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), &
                   drho_mt(:, :, :, id), drho_ir(:, id), dpot_mt(:, :, :, id), dpot_ir(:, id) )
          end if
        end do
        if( from_file ) call dfpt_io_info_string( new_line( 'a' )//'Initial density and potential response read from file.' )
      end if
      ! open files for eigenvalue, occupancy and eigenvector response
      call ph_io_irrep_fxt( iq, iirrep, 0, fxt )
      fdevalk = block_data_file_type( ph_io_qi_dirname//'PHONON_DEVAL_'//trim(fxt)//'.OUT', [nstfv], 1.0_dp )
      fdocck = block_data_file_type( ph_io_qi_dirname//'PHONON_DOCC_'//trim(fxt)//'.OUT', [nstfv], 1.0_dp )
      fdeveck = block_data_file_type( ph_io_qi_dirname//'PHONON_DEVEC_'//trim(fxt)//'.OUT', [nmatmaxkq, nstfv], zone )
      call fdevalk%open( mpilocal )
      call fdocck%open( mpilocal )
      call fdeveck%open( mpilocal )

      call barrier( mpicom=mpilocal )
    end subroutine ph_part_prepare

    !> This subroutine executes finalizing tasks for an independent part 
    !> of a DFPT phonon calculation.
    !>
    !> This includes freeing memory from unneeded variables, closing files,
    !> deleting temporary files and leaving the subdirectory of this part.
    subroutine ph_part_finalize( ipart, &
        info_output, delete_files )
      use mod_misc, only: scrpath
      !> index of the independent part
      integer, intent( in) :: ipart
      !> create and write info output files (default: `.true.`)
      logical, optional, intent(in) :: info_output
      !> delete files for eigensystem response (default: `.true.`)
      logical, optional, intent(in) :: delete_files

      logical :: write_info, delete

      write_info = .true.
      if( present( info_output ) ) write_info = info_output
      delete = .true.
      if( present( delete_files ) ) delete = delete_files

      ! return if this rank is not working on this part
      if( .not. ph_parts(ipart)%is_my_rank( mpiglobal%rank) ) return

      ! close file for eigensystem response
      call fdevalk%close( mpilocal )
      call fdocck%close( mpilocal )
      call fdeveck%close( mpilocal )
      ! delete file for eigensystem response
      if( delete ) then
        call fdevalk%delete( mpilocal )
        call fdocck%delete( mpilocal )
        call fdeveck%delete( mpilocal )
      end if
      ! close temporary files for eigenvectors at k and k+q and delete them
      call feveck%close( mpilocal )
      call feveck%delete( mpilocal )
      call feveckq%close( mpilocal )
      call feveckq%delete( mpilocal )
      ! free unneeded variables
      if( allocated( evalk ) ) deallocate( evalk )
      if( allocated( occk ) ) deallocate( occk )
      if( allocated( eveck ) ) deallocate( eveck )
      if( allocated( apwalmk ) ) deallocate( apwalmk )
      if( allocated( evalkq ) ) deallocate( evalkq )
      if( allocated( occkq ) ) deallocate( occkq )
      if( allocated( eveckq ) ) deallocate( eveckq )
      if( allocated( apwalmkq ) ) deallocate( apwalmkq )
      if( allocated( devalk ) ) deallocate( devalk )
      if( allocated( docck ) ) deallocate( docck )
      if( allocated( deveck ) ) deallocate( deveck )
      if( allocated( dSmat ) ) deallocate( dSmat )
      if( allocated( dHmat_const ) ) deallocate( dHmat_const )
      if( allocated( dHmat ) ) deallocate( dHmat )
      if( allocated( drho_mt ) ) deallocate( drho_mt )
      if( allocated( drho_ir ) ) deallocate( drho_ir )
      if( allocated( drho_mat ) ) deallocate( drho_mat )
      if( allocated( dpot_mt ) ) deallocate( dpot_mt )
      if( allocated( dpot_ir ) ) deallocate( dpot_ir )
      if( allocated( dHmat_mt_basis ) ) deallocate( dHmat_mt_basis )
      if( allocated( dpot_cfun_ig ) ) deallocate( dpot_cfun_ig )
      if( allocated( dkin_cfun_ig ) ) deallocate( dkin_cfun_ig )
      if( allocated( dforce ) ) deallocate( dforce )
      if( allocated( dpol ) ) deallocate( dpol )
      call ph_var_free_q
      ! finalize info output
      if( master .and. write_info ) call ph_io_info_finit
      ph_io_qi_dirname = './'
      scrpath = './'
    end subroutine ph_part_finalize

    !> This subroutine runs the self-consistency cycle for obataining 
    !> the density and potential response for an independent part
    !> of a DFPT phonon calculation.
    !>
    !> Each SCF iteration consists of:
    !>
    !> * for each \({\bf k}\) point and irrep member calculate the
    !>   Hamiltonian matrix response
    !> * for each \({\bf k}\) point and irrep member solve the Sternheimer
    !>   equation for the eigenvalue and eigenvector response
    !> * for each irrep member synthesize the new density response
    !>   and compute the corresponding potential response
    !> * mix old and new density/potential response
    subroutine ph_part_scf( ipart, &
        info_output )
      use dfpt_density_potential, only: dfpt_rhopot_mixpack
      use dfpt_eigensystem, only: dfpt_eig_prepare_dHmat
      use phonons_density_potential, only: ph_rhopot_gen_drho_k, ph_rhopot_gen_drho_mt, ph_rhopot_gen_dpot, ph_rhopot_symmetrize
      use phonons_eigensystem, only: ph_eig_gen_dSHmat, ph_eig_sternheimer
      use phonons_parallelization, only: ph_par_zscatter, ph_par_zgather

      use constants, only: zzero
      use mod_potential_and_density, only: pot_mt => veffmt, pot_ir => veffir
      use mod_convergence, only: iscl
      use mod_eigenvalue_occupancy, only: nstfv, occmax, efermi
      use mod_Gkvector, only: ngkmax_ptr
      use modinput
      !> index of the independent part
      integer, intent(in) :: ipart
      !> write info output files (default: `.true.`)
      logical, optional, intent(in) :: info_output

      integer :: iq, iirrep, dirrep, id, id1, id2, ik, ik1, ik2, i
      integer :: nmix, mixermode
      integer, target :: ngkmaxk, ngkmaxkq
      real(dp) :: conv, t_iter, t1, defermi
      logical :: write_info, success, lastiter, converged, skip, to_file
      character(:), allocatable :: fxt

      integer, allocatable :: rlen(:), roff(:)
      real(dp), allocatable :: rvmix(:), vconv(:)

      write_info = .true.
      if( present( info_output ) ) write_info = info_output

      ! return if this rank is not working on this part
      if( .not. ph_parts(ipart)%is_my_rank( mpiglobal%rank ) ) return

      ! set q-point and irrep
      iq = ph_parts(ipart)%iq
      iirrep = ph_parts(ipart)%iirrep
      dirrep = ph_parts(ipart)%dirrep
      ! set limits for irrep member and k-point loops
      id1 = ph_parts(ipart)%first_d_of_rank( mpiglobal%rank )
      id2 = ph_parts(ipart)%last_d_of_rank( mpiglobal%rank )
      ik1 = ph_parts(ipart)%first_k_of_rank( mpiglobal%rank )
      ik2 = ph_parts(ipart)%last_k_of_rank( mpiglobal%rank )
      ! initialize mixer
      iscl = 0
      if( master ) then
        nmix = 2 * (size( dpot_mt ) + size( dpot_ir ))
        mixermode = -1
        allocate( rvmix(nmix) )
        allocate( vconv(input%groundstate%niterconvcheck) )
        call dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, .true., dirrep, nmix, rvmix )
        call mixerifc( input%groundstate%mixernumber, nmix, rvmix, conv, mixermode )
        conv = 1.0_dp
      end if
      ! allocate record length and offsets for MPI communication
      allocate( rlen(0:mpilocalk%procs-1), source=0 )
      allocate( roff(0:mpilocalk%procs-1), source=0 )
      if( mpilocalk%rank >= 0 ) then
        rlen(mpilocalk%rank) = id2 - id1 + 1 ! record length
        roff(mpilocalk%rank) = id1 - 1       ! record offset
#ifdef MPI
        call MPI_Allreduce( MPI_IN_PLACE, rlen, size( rlen ), MPI_INTEGER, MPI_SUM, mpilocalk%comm, mpilocalk%ierr )
        call MPI_Allreduce( MPI_IN_PLACE, roff, size( roff ), MPI_INTEGER, MPI_SUM, mpilocalk%comm, mpilocalk%ierr )
#endif
      end if
      ! initialize scf cycle info output
      if( master .and. write_info ) call dfpt_io_info_scf_init

      !********************************
      !* START SCF LOOP
      iscl = 0
      lastiter = .false.
      converged = .false.
      scfloop: do while( iscl <= input%groundstate%maxscl .or. iscl == 0 )
        call timesec( t1 )
        ! send effective potential response to other processes
        call ph_par_zscatter( dpot_mt, 0, mpilocalk, rlen*size( dpot_mt(:, :, :, 1) ), roff=roff*size( dpot_mt(:, :, :, 1) ) )
        call ph_par_zscatter( dpot_ir, 0, mpilocalk, rlen*size( dpot_ir(:, 1) ),       roff=roff*size( dpot_ir(:, 1) ) )
        do id = 1, dirrep
          if( ph_parts(ipart)%mpid(id)%rank < 0 ) cycle
          i = 0
          if( mpilocalk%rank >= 0 ) i = ph_parts(ipart)%mpid(id)%rank
#ifdef MPI
          call MPI_Allreduce( MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_SUM, ph_parts(ipart)%mpid(id)%comm, ph_parts(ipart)%mpid(id)%ierr )
          call MPI_Bcast( dpot_mt(:, :, :, id), size( dpot_mt(:, :, :, id) ), MPI_DOUBLE_COMPLEX, &
                          i, ph_parts(ipart)%mpid(id)%comm, ph_parts(ipart)%mpid(id)%ierr )
          call MPI_Bcast( dpot_ir(:, id), size( dpot_ir(:, id)), MPI_DOUBLE_COMPLEX, &
                          i, ph_parts(ipart)%mpid(id)%comm, ph_parts(ipart)%mpid(id)%ierr )
#endif
        end do
        ! send occupation response to other processes
#ifdef MPI
        call MPI_Bcast( docck, size( docck ), MPI_DOUBLE, localmasterrank, mpilocal%comm, mpilocal%ierr )
#endif
        ! increment iteration and check for last iteration
        iscl = iscl + 1
        skip = (lastiter .or. iscl > input%groundstate%maxscl)
        if( lastiter ) exit scfloop
        lastiter = converged .or. (iscl >= input%groundstate%maxscl)

        ! perform scf cycle
        if( .not. skip ) then
          drho_mat = zzero
          drho_ir = zzero
          devalk = 0.0_dp
          ngkmaxk = ph_Gkset%ngkmax
          ngkmaxkq = ph_Gkset%ngkmax
          do ik = ik1, ik2
            ! get matching coefficients at k
            ngkmax_ptr => ngkmaxk
            call match( ph_Gkset%ngk(1, ik), ph_Gkset%gkc(:, 1, ik), ph_Gkset%tpgkc(:, :, 1, ik), ph_Gkset%sfacgk(:, :, 1, ik), &
                        apwalmk )
            ! get matching coefficients at k+q
            ngkmax_ptr => ngkmaxkq
            call match( ph_Gkqset%ngk(1, ik), ph_Gkqset%gkc(:, 1, ik), ph_Gkqset%tpgkc(:, :, 1, ik), ph_Gkqset%sfacgk(:, :, 1, ik), &
                        apwalmkq )
            ! read eigenvectors
            call feveck%read( ik, eveck )
            call feveckq%read( ik, eveckq )
            do id = id1, id2
              ! generate k-independent part of Hamiltonian response
              if( ik == ik1 ) &
                call dfpt_eig_prepare_dHmat( pot_mt, pot_ir, dpot_mt(:, :, :, id), dpot_ir(:, id), dHmat_mt_basis(:, :, :, id), dpot_cfun_ig(:, id), dkin_cfun_ig(:, id), &
                                             Gset=ph_Gqset )
              ! generate full Hamiltonian response
              dHmat = dHmat_const(:, :, ik, id)
              call ph_eig_gen_dSHmat( ik, ph_Gkset, ph_Gkqset, 1, nstfv, &
                     eveck, eveckq, apwalmk, apwalmkq, dSmat(:, :, ik, id), dHmat, &
                     dHmat_mt_basis=dHmat_mt_basis(:, :, :, id), dpot_cfun_ig=dpot_cfun_ig(:, id), dkin_cfun_ig=dkin_cfun_ig(:, id) )
              ! solve Sternheimer equation
              call ph_eig_sternheimer( ik, ph_Gkqset, 1, nstfv, &
                     evalk(:, ik), occk(:, ik), evalkq(:, ik), occkq(:, ik), eveckq, dSmat(:, :, ik, id), dHmat, gamma, &
                     devalk(:, ik, id), deveck, projector=.true. )
              ! write eigenvector response to file
              call fdeveck%write( ph_parts(ipart)%get_dk_offset(id, ik)+1, deveck )
              ! add k-point contribution to density response
              call ph_rhopot_gen_drho_k( ik, ph_kset, ph_Gkset, ph_Gkqset, 1, nstfv, &
                     occk(:, ik), docck(:, ik, id), eveck, deveck, apwalmk, apwalmkq, &
                     ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), gamma, &
                     drho_mat(:, :, :, id), drho_ir(:, id) )
            end do
          end do
          call barrier( mpicom=mpilocal )
          do id = id1, id2
            ! reduce density response over k-points
#ifdef MPI
            call MPI_Allreduce( MPI_IN_PLACE, drho_mat(:, :, :, id), size( drho_mat(:, :, :, id) ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                   ph_parts(ipart)%mpid(id)%comm, ph_parts(ipart)%mpid(id)%ierr )
            call MPI_Allreduce( MPI_IN_PLACE, drho_ir(:, id), size( drho_ir(:, id) ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                   ph_parts(ipart)%mpid(id)%comm, ph_parts(ipart)%mpid(id)%ierr )
#endif
            ! update muffin-tin density response
            call ph_rhopot_gen_drho_mt( ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), drho_mat(:, :, :, id), drho_mt(:, :, :, id) )
          end do
          ! reduce eigenvalue response over k-points
#ifdef MPI
          call MPI_Allreduce( MPI_IN_PLACE, devalk, size( devalk ), MPI_DOUBLE, MPI_SUM, &
                 mpilocal%comm, mpilocal%ierr )
#endif
          ! collect density response at local master
          call ph_par_zgather( drho_mt, 0, mpilocalk, rlen*size( drho_mt(:, :, :, 1) ), roff=roff*size( drho_mt(:, :, :, 1) ) )
          call ph_par_zgather( drho_ir, 0, mpilocalk, rlen*size( drho_ir(:, 1) ),       roff=roff*size( drho_ir(:, 1) ) )
          if( master ) then
            ! symmetrize density response
            call ph_rhopot_symmetrize( drho_mt, drho_ir, ph_qset%vkl(:, iq), dirrep, &
                   ph_irrep_basis(iq)%nsym, ph_irrep_basis(iq)%isym, ph_irrep_basis(iq)%ivsym, &
                   ph_irrep_basis(iq)%irreps(iirrep)%symmat )
            ! update effective potential response
            do id = 1, dirrep
              call ph_rhopot_gen_dpot( ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), &
                     drho_mt(:, :, :, id), drho_ir(:, id), dpot_mt(:, :, :, id), dpot_ir(:, id) )
            end do
            ! symmetrize potential response
            call ph_rhopot_symmetrize( dpot_mt, dpot_ir, ph_qset%vkl(:, iq), dirrep, &
                   ph_irrep_basis(iq)%nsym, ph_irrep_basis(iq)%isym, ph_irrep_basis(iq)%ivsym, &
                   ph_irrep_basis(iq)%irreps(iirrep)%symmat )
            ! uppdate occupation response
            ! and write eigenvalue and occupation response to file
            if( gamma ) then
              do id = 1, dirrep
                call find_dfermi( ph_kset%nkpt, ph_kset%wkpt, nstfv, evalk, devalk(:, :, id), 0.0_dp, occmax, efermi, &
                       input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
                       defermi, docck(:, :, id) )
                do ik = 1, ph_kset%nkpt
                  call fdevalk%write( ph_parts(ipart)%get_dk_offset(id, ik) + 1, devalk(:, ik, id) )
                  call fdocck%write( ph_parts(ipart)%get_dk_offset(id, ik) + 1, docck(:, ik, id) )
                end do
              end do
            end if
          end if
        end if ! end of scf step

        ! write density response to file
        if( master ) then
          to_file = .false.
          if( input%groundstate%nwrite > 0 ) to_file = (mod( max( 1, iscl - 1 ), input%groundstate%nwrite ) == 0)
          to_file = to_file .or. lastiter
          if( to_file ) then
            do id = 1, dirrep
              call ph_io_irrep_fxt( iq, iirrep, id, fxt )
              call dfpt_io_write_zfun( drho_mt(:, :, :, id), drho_ir(:, id), 'PHONON_DRHO', success, &
                file_extension=fxt, directory=ph_io_qi_dirname )
              to_file = to_file .and. success
              call dfpt_io_write_zfun( dpot_mt(:, :, :, id), dpot_ir(:, id), 'PHONON_DVEFF', success, &
                file_extension=fxt, directory=ph_io_qi_dirname )
              to_file = to_file .and. success
            end do
          end if
          if( to_file ) call dfpt_io_info_string( 'Density and potential response written to file.' )
        end if

        ! exit on demand
        if( skip .and. lastiter ) exit scfloop

        ! mixing
        if( master ) then
          call dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, .true., dirrep, nmix, rvmix )
          call mixerifc( input%groundstate%mixernumber, nmix, rvmix, conv, mixermode )
          do i = 1, input%groundstate%niterconvcheck - 1
            vconv(i) = vconv(i+1)
          end do
          vconv(i) = conv
          ! check convergence
          converged = all( abs( vconv ) < input%groundstate%epspot )
          if( .not. lastiter ) &
            call dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, .false., dirrep, nmix, rvmix )
        end if
#ifdef MPI
        call MPI_Bcast( lastiter, 1, MPI_LOGICAL, localmasterrank, mpilocal%comm, mpilocal%ierr )
        call MPI_Bcast( converged, 1, MPI_LOGICAL, localmasterrank, mpilocal%comm, mpilocal%ierr )
#endif
        call timesec( t_iter )
        t_iter = t_iter - t1

        ! write iteration info to output file
        if( master .and. write_info ) then
          if( gamma ) then
            call dfpt_io_info_scf( iscl, conv, t_iter, defermi=defermi )
          else
            call dfpt_io_info_scf( iscl, conv, t_iter )
          end if
        end if

        call barrier( mpicom=mpilocal )
      end do scfloop
      !* END SCF LOOP
      !********************************

      ! finalize scf loop output
      if( master .and. write_info ) call dfpt_io_info_scf_finit

      ! deallocate mixer
      deallocate( rlen, roff )
      if( master ) then
        mixermode = -2
        call mixerifc( input%groundstate%mixernumber, nmix, rvmix, conv, mixermode )
        deallocate( rvmix, vconv )
      end if

      call barrier( mpicom=mpilocal )
    end subroutine ph_part_scf

    !> This subroutine calculates the force response (= dynamical matrix column) 
    !> from the converged density and potential response for an independent part
    !> of a DFPT phonon calculation and writes it to file.
    !>
    !> The force response consists of the response to the Hellmann-Feynman force,
    !> the Pulay force and the contribution coming from surface integrals.
    !>
    !> See [[ph_dynmat_canonical_from_file(subroutine)]] for the transformation of the dynamical 
    !> matrix into canonical (Cartesian) coordinates.
    subroutine ph_part_force( ipart )
      use dfpt_eigensystem, only: dfpt_eig_prepare_dHmat
      use phonons_eigensystem, only: ph_eig_gen_dSHmat, ph_eig_sternheimer
      use phonons_parallelization, only: ph_par_zgather
      use phonons_density_potential, only: ph_rhopot_gen_dpot
      use phonons_force
      use phonons_io_util, only: ph_io_write_dynmat_col

      use constants, only: zzero, zone
      use mod_potential_and_density, only: pot_mt => veffmt, pot_ir => veffir
      use mod_eigenvalue_occupancy, only: nstfv
      use mod_Gkvector, only: ngkmax_ptr
      use mod_atoms, only: natmtot
      use modinput
      !> index of the independent part
      integer, intent(in) :: ipart

      integer :: iq, iirrep, dirrep, id, id1, id2, ik, ik1, ik2, i
      integer, target :: ngkmaxk, ngkmaxkq
      logical :: success, to_file
      character(:), allocatable :: fxt

      integer, allocatable :: rlen(:), roff(:)

      ! return if this rank is not working on this part
      if( .not. ph_parts(ipart)%is_my_rank( mpiglobal%rank ) ) return

      ! set q-point and irrep
      iq = ph_parts(ipart)%iq
      iirrep = ph_parts(ipart)%iirrep
      dirrep = ph_parts(ipart)%dirrep
      ! set limits for irrep member and k-point loops
      id1 = ph_parts(ipart)%first_d_of_rank( mpiglobal%rank )
      id2 = ph_parts(ipart)%last_d_of_rank( mpiglobal%rank )
      ik1 = ph_parts(ipart)%first_k_of_rank( mpiglobal%rank )
      ik2 = ph_parts(ipart)%last_k_of_rank( mpiglobal%rank )
      ! allocate record length and offsets for MPI communication
      allocate( rlen(0:mpilocalk%procs-1), source=0 )
      allocate( roff(0:mpilocalk%procs-1), source=0 )
      if( mpilocalk%rank >= 0 ) then
        rlen(mpilocalk%rank) = id2 - id1 + 1 ! record length
        roff(mpilocalk%rank) = id1 - 1       ! record offset
#ifdef MPI
        call MPI_Allreduce( MPI_IN_PLACE, rlen, size( rlen ), MPI_INTEGER, MPI_SUM, mpilocalk%comm, mpilocalk%ierr )
        call MPI_Allreduce( MPI_IN_PLACE, roff, size( roff ), MPI_INTEGER, MPI_SUM, mpilocalk%comm, mpilocalk%ierr )
#endif
      end if

      dforce = zzero
      ngkmaxk = ph_Gkset%ngkmax
      ngkmaxkq = ph_Gkset%ngkmax
      ! k-point contribution
      do ik = ik1, ik2
        ! get matching coefficients at k
        ngkmax_ptr => ngkmaxk
        call match( ph_Gkset%ngk(1, ik), ph_Gkset%gkc(:, 1, ik), ph_Gkset%tpgkc(:, :, 1, ik), ph_Gkset%sfacgk(:, :, 1, ik), &
                    apwalmk )
        ! get matching coefficients at k+q
        ngkmax_ptr => ngkmaxkq
        call match( ph_Gkqset%ngk(1, ik), ph_Gkqset%gkc(:, 1, ik), ph_Gkqset%tpgkc(:, :, 1, ik), ph_Gkqset%sfacgk(:, :, 1, ik), &
                    apwalmkq )
        ! read eigenvectors at k and k+q
        call feveck%read( ik, eveck )
        call feveckq%read( ik, eveckq )
        do id = id1, id2
          ! generate k-independent part of Hamiltonian response
          if( ik == ik1 ) &
            call dfpt_eig_prepare_dHmat( pot_mt, pot_ir, dpot_mt(:, :, :, id), dpot_ir(:, id), dHmat_mt_basis(:, :, :, id), dpot_cfun_ig(:, id), dkin_cfun_ig(:, id), &
                                         Gset=ph_Gqset )
          ! generate full Hamiltonian response
          dHmat = dHmat_const(:, :, ik, id)
          call ph_eig_gen_dSHmat( ik, ph_Gkset, ph_Gkqset, 1, nstfv, &
                 eveck, eveckq, apwalmk, apwalmkq, dSmat(:, :, ik, id), dHmat, &
                 dHmat_mt_basis=dHmat_mt_basis(:, :, :, id), dpot_cfun_ig=dpot_cfun_ig(:, id), dkin_cfun_ig=dkin_cfun_ig(:, id) )
          ! solve Sternheimer equation without projection
          call ph_eig_sternheimer( ik, ph_Gkqset, 1, nstfv, &
                 evalk(:, ik), occk(:, ik), evalkq(:, ik), occkq(:, ik), eveckq, dSmat(:, :, ik, id), dHmat, gamma, &
                 devalk(:, ik, id), deveck, projector=.false. )
          ! add k-point contribution to force response
          call ph_frc_dpulay_k( ik, ph_kset, ph_Gkset, ph_Gkqset, 1, nstfv, &
                 evalk(:, ik), devalk(:, ik, id), occk(:, ik), docck(:, ik, id), eveck, deveck, apwalmk, apwalmkq, &
                 ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), gamma, &
                 dforce(:, :, id), &
                 dHmat_mt_basis=dHmat_mt_basis(:, :, :, id) )
          call ph_frc_dsurf_k( ik, ph_kset, ph_Gkset, ph_Gkqset, 1, nstfv, &
                 evalk(:, ik), devalk(:, ik, id), occk(:, ik), docck(:, ik, id), eveck, deveck, &
                 ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), gamma, &
                 dforce(:, :, id), &
                 dpot_ir=dpot_ir(:, id) )
        end do
      end do
      call barrier( mpicom=mpilocal )
      ! reduce force response over k-points
      do id = id1, id2
#ifdef MPI
        call MPI_Allreduce( MPI_IN_PLACE, dforce(:, :, id), size( dforce(:, :, id) ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
               ph_parts(ipart)%mpid(id)%comm, ph_parts(ipart)%mpid(id)%ierr )
#endif
      end do
      ! collect force response at local master
      call ph_par_zgather( dforce, 0, mpilocalk, rlen*size( dforce(:, :, 1) ), roff=roff*size( dforce(:, :, 1) ) )
      if( master ) then
        do id = 1, dirrep
          ! add integral contribution from Pulay and surface force response
          call ph_frc_dpulay_int( drho_mt(:, :, :, id), dpot_mt(:, :, :, id), dforce(:, :, id) )
          ! get soft Coulomb potential response
          call ph_rhopot_gen_dpot( ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), &
                 drho_mt(:, :, :, id), drho_ir(:, id), dpot_mt(:, :, :, id), dpot_ir(:, id), xc=.false. )
          ! add Hellmann-Feynman force response
          call ph_frc_dhf( drho_mt(:, :, :, id), dpot_mt(:, :, :, id), dforce(:, :, id) )
          ! add constant part from basis function response
          do i = 1, natmtot
            call zgemv( 'n', 3, 3, zone, dforce_const(1, i, 1), 3*natmtot, &
                   ph_irrep_basis(iq)%irreps(iirrep)%pat(1, i, id), 1, zone, &
                   dforce(1, i, id), 1 )
          end do
        end do
        ! symmetrize force response
        call ph_frc_symmetrize( dforce, ph_qset%vkl(:, iq), dirrep, &
               ph_irrep_basis(iq)%nsym, ph_irrep_basis(iq)%isym, ph_irrep_basis(iq)%ivsym, &
               ph_irrep_basis(iq)%irreps(iirrep)%symmat, &
               acoustic_sum_rule=(gamma .and. input%phonons%sumrule) )
        ! write dynamical matrix column in irrep basis to file
        to_file = .true.
        do id = 1, dirrep
          call ph_io_irrep_fxt( iq, iirrep, id, fxt )
          call ph_io_write_dynmat_col( dforce(:, :, id), fxt, success, directory=ph_io_qi_dirname )
          to_file = to_file .and. success
        end do
        if( to_file ) call dfpt_io_info_string( 'Dynamical matrix columns written to file.' )
      end if

      deallocate( rlen, roff )

      call barrier( mpicom=mpilocal )
    end subroutine ph_part_force

    !> This subroutine calculates the polarization response (= Born effective charge column) 
    !> from the converged density and potential response for an independent part
    !> of a DFPT phonon calculation and writes it to file.
    !>
    !> The polarization response is calculated using the Berry phase approach.
    !>
    !> See [[ph_borncharge_canonical_from_file(subroutine)]] for the transformation of the Born effective 
    !> charge into canonical (Cartesian) coordinates.
    subroutine ph_part_polarization( ipart )
      use dfpt_polarization
      use phonons_polarization
      use phonons_symmetry, only: ph_sym_find_kpt, ph_sym_rotate_devec

      use constants, only: zzero
      use mod_eigenvalue_occupancy, only: occmax
      !> index of the independent part
      integer, intent(in) :: ipart

      integer :: iq, iirrep, dirrep, id, ik, ikb, ik1, ik2, ip
      integer :: iknr, ikbnr, isym, nst, ivkb(3)
      real(dp) :: binv(3, 3), tmp(3, 2)
      logical :: success, to_file
      character(:), allocatable :: fxt

      complex(dp), allocatable :: devecknr(:,:,:), eveckbnr(:,:), deveckbnr(:,:,:)

      ! return if this rank is not working on this part
      if( .not. ph_parts(ipart)%is_my_rank( mpiglobal%rank ) ) return

      ! return if q is different from 0
      if( .not. gamma ) return
      
      ! set q-point and irrep
      iq = ph_parts(ipart)%iq
      iirrep = ph_parts(ipart)%iirrep
      dirrep = ph_parts(ipart)%dirrep
      ! set limits for k-point loops
      ik1 = firstofset( mpilocal%rank, ph_kset%nkptnr, nprocs=mpilocal%procs )
      ik2 = lastofset( mpilocal%rank, ph_kset%nkptnr, nprocs=mpilocal%procs )

      allocate( devecknr(size( deveck, dim=1 ), size( deveck, dim=2 ), dirrep) )
      allocate( eveckbnr(size( eveck, dim=1 ), size( eveck, dim=2 )) )
      allocate( deveckbnr(size( deveck, dim=1 ), size( deveck, dim=2 ), dirrep) )

      dpol = zzero
      ! initialize plane wave matrix elements
      call ph_pol_init( ph_irrep_basis(iq)%irreps(iirrep) )
      do ip = 1, 3
        ! initialize plane wave matrix elements
        call dfpt_pol_init_p( - ph_kset%bvec(:, ip) / ph_kset%ngridk(ip) )
        ! we have to loop over the full k-set
        do iknr = ik1, ik2
          ! find equivalent point in reduced set and corresponding symmetry
          call ph_sym_find_kpt( ph_kset%vklnr(:, iknr), ph_kset%vkl, ph_kset%nkpt, ph_irrep_basis(iq), isym, ik )
          call terminate_if_false( ik > 0, '(ph_part_polarization) &
            Equivalent k-point not found.' )
          ! get number of occupied states
          nst = count( occk(:, ik) > occmax/2 )
          ! read eigenvector at k
          call feveck%read( ik, eveck )
          call rotate_evecfv( ph_irrep_basis(iq)%isym(isym), ph_kset%vkl(:, ik), ph_kset%vklnr(:, iknr), &
            ph_Gkset%ngk(1, ik), ph_Gkset%vgkl(:, :, 1, ik), ph_Gkset%vgknrl(:, :, 1, iknr), &
            eveck, size( eveck, dim=1 ), nst )
          ! read eigenvector response at k
          do id = 1, dirrep
            call fdeveck%read( ph_parts(ipart)%get_dk_offset(id, ik) + 1, devecknr(:, :, id) )
          end do
          call ph_sym_rotate_devec( ph_irrep_basis(iq), iirrep, isym, ph_kset%vkl(:, ik), ph_kset%vklnr(:, iknr), &
            ph_Gkset%ngk(1, ik), ph_Gkset%vgkl(:, :, 1, ik), ph_Gkset%vgknrl(:, :, 1, iknr), &
            devecknr, size( devecknr, dim=1 ), size( devecknr, dim=2 ), nst )
          ! find vector k+b and its equivalent in reduced set
          ivkb = ph_kset%ivknr(:, iknr)
          ivkb(ip) = modulo( ivkb(ip) + 1, ph_kset%ngridk(ip) )
          ikbnr = ph_kset%ikmapnr(ivkb(1), ivkb(2), ivkb(3))
          call ph_sym_find_kpt( ph_kset%vklnr(:, ikbnr), ph_kset%vkl, ph_kset%nkpt, ph_irrep_basis(iq), isym, ikb )
          call terminate_if_false( ikb > 0, '(ph_part_polarization) &
            Equivalent k-point not found.' )
          ! read eigenvector at k+b
          call feveck%read( ikb, eveckbnr )
          call rotate_evecfv( ph_irrep_basis(iq)%isym(isym), ph_kset%vkl(:, ikb), ph_kset%vklnr(:, ikbnr), &
            ph_Gkset%ngk(1, ikb), ph_Gkset%vgkl(:, :, 1, ikb), ph_Gkset%vgknrl(:, :, 1, ikbnr), &
            eveckbnr, size( eveckbnr, dim=1 ), nst )
          ! read eigenvector response at k+b
          do id = 1, dirrep
            call fdeveck%read( ph_parts(ipart)%get_dk_offset(id, ikb) + 1, deveckbnr(:, :, id) )
          end do
          call ph_sym_rotate_devec( ph_irrep_basis(iq), iirrep, isym, ph_kset%vkl(:, ikb), ph_kset%vklnr(:, ikbnr), &
            ph_Gkset%ngk(1, ikb), ph_Gkset%vgkl(:, :, 1, ikb), ph_Gkset%vgknrl(:, :, 1, ikbnr), &
            deveckbnr, size( deveckbnr, dim=1 ), size( deveckbnr, dim=2 ), nst )

          ! add k-point contribution to polarization response
          do id = 1, dirrep
            call ph_pol_dpol_k( iknr, ikbnr, ph_kset, ph_Gkset, 1, nst, &
              eveck, eveckbnr, devecknr(:, :, id), deveckbnr(:, :, id), &
              ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id), id, &
              dpol(ip, id) )
          end do

        end do
      end do
      ! free memory
      call ph_pol_free
      call barrier( mpicom=mpilocal )
      ! reduce polarization response over k-points
#ifdef MPI
      call MPI_Allreduce( MPI_IN_PLACE, dpol, size( dpol ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
             mpilocal%comm, mpilocal%ierr )
#endif
      ! transform into Cartesian coordinates
      call r3minv( ph_kset%bvec, binv )
      do id = 1, dirrep
        dpol(:, id) = dpol(:, id) * ph_kset%ngridk
        call r3mtv( binv, dble(  dpol(:, id) ), tmp(:, 1) )
        call r3mtv( binv, aimag( dpol(:, id) ), tmp(:, 2) )
        dpol(:, id) = cmplx( tmp(:, 1), tmp(:, 2), dp )
      end do
      ! add ionic contribution
      call ph_pol_dpol_ion( dirrep, ph_irrep_basis(iq)%irreps(iirrep)%pat, dpol )
      ! symmetrize
      call ph_pol_symmetrize( dpol, dirrep, &
             ph_irrep_basis(iq)%nsym, ph_irrep_basis(iq)%isym, ph_irrep_basis(iq)%ivsym, &
             ph_irrep_basis(iq)%irreps(iirrep)%symmat )
      ! write to file
      if( master ) then
        to_file = .true.
        do id = 1, dirrep
          call ph_io_irrep_fxt( iq, iirrep, id, fxt )
          call ph_io_write_borncharge_col( dpol(:, id), fxt, success, directory=ph_io_qi_dirname )
          to_file = to_file .and. success
        end do
        if( to_file ) call dfpt_io_info_string( 'Born effective charge columns written to file.' )
      end if

      deallocate( devecknr, eveckbnr, deveckbnr )

      call barrier( mpicom=mpilocal )
    end subroutine ph_part_polarization

    !> This subroutine tries to read the dynamical matrices in (half) irrep coordinates
    !> from a previous SCF run and transforms them into the canonical (Cartesian) basis
    !> and writes them to file.
    !>
    !> The dynamical matrices in canonical (Cartesian) coordinates are obtained
    !> from the irrep displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\) by
    !> \[ D_{\kappa\alpha,\kappa'\beta}({\bf q}) = 
    !>    \sum_{I, \mu} {p^{I \mu}_{\kappa'\beta}}^\ast ({\bf q}) \, D_{\kappa\alpha,I \mu}({\bf q}) \;. \]
    subroutine ph_write_dyn_canonical
      use phonons_io_util, only: ph_io_write_dynmat_col
      use mod_atoms, only: natmtot, nspecies, natoms, idxas

      integer :: iq, is, ia, ias, ip, i
      character(:), allocatable :: fxt
      logical :: success

      complex(dp), allocatable :: dyn(:,:)

      if( mpiglobal%rank /= 0 ) return

      do iq = 1, ph_qset%nkpt
        call ph_dynmat_canonical_from_file( ph_qset%vkl(:, iq), dyn )
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            do ip = 1, 3
              i = (ias - 1) * 3 + ip
              ! get file extension
              call ph_io_canon_fxt( iq, is, ia, ip, fxt )
              ! try to write dynamical matrix to file
              call ph_io_write_dynmat_col( reshape( dyn(:, i), [3, natmtot] ), fxt, success )
              call terminate_if_false( success, '(ph_write_dyn_canonical) &
                Was not able to write dynamical matrix column to file '//&
                new_line( 'a' )//'DYN_'//trim( fxt )//'.OUT' )
            end do
          end do
        end do
      end do

      if( allocated( dyn ) ) deallocate( dyn )
    end subroutine ph_write_dyn_canonical

    !> For a given phonon wavevector \({\bf q}\), this subroutine tries to read the dynamical matrix
    !> in irrep coordinates from a previous SCF run and transforms it into the canonical (Cartesian) basis.
    !>
    !> First, for the given \({\bf q}\) its symmetry equivalent point \({\bf q}_0\) in the set of
    !> \({\bf q}\)-vectors and the connecting symmetry operation is found. Then, the dynamical matrix
    !> in irrep coordinates at \({\bf q}_0\) is read from file, transformed into canonical coordinates via
    !> \[ D_{\kappa\alpha,\kappa'\beta}({\bf q}_0) = 
    !>    \sum_{I, \mu} {p^{I \mu}_{\kappa'\beta}}^\ast ({\bf q}_0) \, D_{\kappa\alpha,I \mu}({\bf q}_0) \;, \]
    !> where \(p^{I \mu}_{\kappa\alpha}({\bf q})\) are the irrep displacement patterns,
    !> and possibly rotated into \({\bf q}\) using [[ph_util_symapp_dyn(subroutine)]].
    subroutine ph_dynmat_canonical_from_file( vql, dyn )
      use constants, only: zzero
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use phonons_util, only: ph_util_symapp_dyn
      use phonons_io_util, only: ph_io_read_dynmat_col
      !> phonon wavevector \({\bf q}\) in lattice coordinates
      real(dp), intent(in) :: vql(3)
      !> dynamical matrix
      complex(dp), allocatable, intent(out) :: dyn(:,:)

      integer :: iq, isym, iirrep, dirrep, id, &
                 is, ia, ias, ip, i
      complex(dp) :: pat(3, natmtot)
      character(:), allocatable :: qidirname, fxt
      logical :: success
      
      complex(dp), allocatable :: dyn_i(:,:)

      ! allocate output arrays
      if( allocated( dyn ) ) deallocate( dyn )
      allocate( dyn(3*natmtot, 3*natmtot), source=zzero )
      ! allocate local arrays
      allocate( dyn_i(3, natmtot) )

      ! find equivalent q-point q0 in set and connecting symmetry
      call findkptinset( vql, ph_qset, isym, iq )

      ! read dyamical matrix in irrep basis at q0
      ! and transform it into canonical basis
      do iirrep = 1, ph_irrep_basis(iq)%nirrep
        dirrep = ph_irrep_basis(iq)%irreps(iirrep)%dim
        ! get subdirectory of part
        call ph_io_irrep_fxt( iq, iirrep, 0, qidirname )
        do id = 1, dirrep
         ! get file extension of part
          call ph_io_irrep_fxt( iq, iirrep, id, fxt )
          ! try to read dynamical matrix from file
          call ph_io_read_dynmat_col( dyn_i, fxt, success, directory=qidirname )
          call terminate_if_false( success, '(ph_dynmat_canonical_from_file) &
            Was not able to read dynamical matrix column supposed to be in file '//&
            new_line( 'a' )//trim( qidirname )//'/DYN_'//trim(fxt)//'.OUT' )
          ! store irrep displacement pattern
          pat = ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id)
          ! add contribution to dynamical matrix in canonical basis
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)
              do ip = 1, 3
                i = (ias - 1) * 3 + ip
                dyn(:, i) = dyn(:, i) + conjg( pat(ip, ias) ) * reshape( dyn_i, [3*natmtot] )
              end do
            end do
          end do
        end do
      end do
      deallocate( dyn_i )

      ! rotate dynamical matrix from q0 into q if symmetry is not the identity
      if( isym == 1 ) return
      allocate( dyn_i, source=dyn )
      dyn = zzero
      call ph_util_symapp_dyn( isym, ph_qset%vkl(:, iq), dyn_i, dyn )
      deallocate( dyn_i )
    end subroutine ph_dynmat_canonical_from_file

    !> This subroutine tries to read the effective potential response in irrep coordinates
    !> from a previous SCF run and transforms it into the canonical (Cartesian) basis
    !> and writes it to file.
    !>
    !> The effective potential response in canonical (Cartesian) coordinates is obtained
    !> from the irrep displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\) by
    !> \[ \delta^{\bf q}_{\kappa\alpha}V_{\rm eff}({\bf r}) = 
    !>    \sum_{I, \mu} {p^{I \mu}_{\kappa\alpha}}^\ast ({\bf q}) \, \delta^{\bf q}_{I \mu}V_{\rm eff}({\bf r}) \;. \]
    subroutine ph_write_dpot_canonical
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_phonon, only: natoms0, natmtot0, ngrid0, ngrtot0

      integer :: iq, is, ia, ias, ip

      complex(dp), allocatable :: dpot_mt(:,:,:,:,:), dpot_ir(:,:,:)

      if( mpiglobal%rank /= 0 ) return

      ! for needed super cell phonon globals
      natoms0(1:nspecies) = natoms(1:nspecies)
      natmtot0 = natmtot
      ngrid0 = dfpt_Gset%ngrid
      ngrtot0 = dfpt_Gset%ngrtot
      call init2

      do iq = 1, ph_qset%nkpt
        call ph_dpot_canonical_from_file( ph_qset%vkl(:, iq), dpot_mt, dpot_ir )
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            do ip = 1, 3
              call writedveff( iq, is, ia, ip, dpot_mt(:, :, :, ip, ias), dpot_ir(:, ip, ias) )
            end do
          end do
        end do
      end do

      if( allocated( dpot_mt ) ) deallocate( dpot_mt )
      if( allocated( dpot_ir ) ) deallocate( dpot_ir )
    end subroutine ph_write_dpot_canonical

    !> For a given phonon wavevector \({\bf q}\), this subroutine tries to read the effective potential
    !> response in irrep coordinates from a previous SCF run and transforms it into the canonical (Cartesian) basis.
    !>
    !> First, for the given \({\bf q}\) its symmetry equivalent point \({\bf q}_0\) in the set of
    !> \({\bf q}\)-vectors and the connecting symmetry operation is found. Then, the potential response
    !> in irrep coordinates at \({\bf q}_0\) is read from file, transformed into canonical coordinates via
    !> \[ \delta^{{\bf q}_0}_{\kappa\alpha}V_{\rm eff}({\bf r}) = 
    !>    \sum_{I, \mu} {p^{I \mu}_{\kappa\alpha}}^\ast ({\bf q}_0) \, \delta^{{\bf q}_0}_{I \mu}V_{\rm eff}({\bf r}) \;, \]
    !> where \(p^{I \mu}_{\kappa\alpha}({\bf q})\) are the irrep displacement patterns,
    !> and possibly rotated into \({\bf q}\) using [[ph_rhopot_rotate_q_canonical(subroutine)]].
    subroutine ph_dpot_canonical_from_file( vql, dpot_mt, dpot_ir )
      use phonons_density_potential, only: ph_rhopot_rotate_q_canonical
      use constants, only: zzero
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_muffin_tin, only: nrmtmax
      !> phonon wavevector \({\bf q}\) in lattice coordinates
      real(dp), intent(in) :: vql(3)
      !> muffin-tin potential response
      complex(dp), allocatable, intent(out) :: dpot_mt(:,:,:,:,:)
      !> interstitial potential response
      complex(dp), allocatable, intent(out) :: dpot_ir(:,:,:)

      integer :: iq, isym, iirrep, dirrep, id, &
                 is, ia, ias, ip
      complex(dp) :: pat(3, natmtot)
      character(:), allocatable :: qidirname, fxt
      logical :: success
      
      complex(dp), allocatable :: dpot_mt_i(:,:,:), dpot_ir_i(:)

      ! allocate output arrays
      if( allocated( dpot_mt ) ) deallocate( dpot_mt )
      allocate( dpot_mt(dfpt_lmmaxvr, nrmtmax, natmtot, 3, natmtot), source=zzero )
      if( allocated( dpot_ir ) ) deallocate( dpot_ir )
      allocate( dpot_ir(dfpt_Gset%ngrtot, 3, natmtot), source=zzero )
      ! allocate local arrays
      allocate( dpot_mt_i(dfpt_lmmaxvr, nrmtmax, natmtot) )
      allocate( dpot_ir_i(dfpt_Gset%ngrtot) )

      ! find equivalent q-point q0 in set and connecting symmetry
      call findkptinset( vql, ph_qset, isym, iq )

      ! read potential response in irrep basis at q0
      ! and transform it into canonical basis
      do iirrep = 1, ph_irrep_basis(iq)%nirrep
        dirrep = ph_irrep_basis(iq)%irreps(iirrep)%dim
        ! get subdirectory of part
        call ph_io_irrep_fxt( iq, iirrep, 0, qidirname )
        do id = 1, dirrep
          ! get file extension of part
          call ph_io_irrep_fxt( iq, iirrep, id, fxt )
          ! try to read potential response from file
          call dfpt_io_read_zfun( dpot_mt_i, dpot_ir_i, 'PHONON_DVEFF', success, &
                 file_extension=fxt, directory=qidirname )
          call terminate_if_false( success, '(ph_dpot_canonical_from_file) &
            Was not able to read potential response supposed to be in file '//&
            new_line( 'a' )//trim( qidirname )//'/PHONON_DVEFF_'//trim(fxt)//'.OUT' )
          ! store irrep displacement pattern
          pat = ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id)
          ! add contribution to potential response in canonical basis
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)
              do ip = 1, 3
                dpot_mt(:, :, :, ip, ias) = dpot_mt(:, :, :, ip, ias) + conjg( pat(ip, ias) ) * dpot_mt_i
                dpot_ir(:, ip, ias) = dpot_ir(:, ip, ias) + conjg( pat(ip, ias) ) * dpot_ir_i
              end do
            end do
          end do
        end do
      end do

      deallocate( dpot_mt_i, dpot_ir_i )

      ! rotate potential response from q0 into q if symmetry is not the identity
      if( isym /= 1 ) &
        call ph_rhopot_rotate_q_canonical( ph_qset%vkl(:, iq), vql, isym, dpot_mt, dpot_ir )
    end subroutine ph_dpot_canonical_from_file

    !> This subroutine tries to read the Born effective charges in irrep coordinates
    !> from a previous SCF run and transforms them into the canonical (Cartesian) basis
    !> and writes them to file.
    !>
    !> The Born effective charges in canonical (Cartesian) coordinates are obtained
    !> from the irrep displacement patterns \(p^{I \mu}_{\kappa\alpha}({\bf q})\) by
    !> \[ Z^{\ast}_{i,\kappa'\beta} = 
    !>    \sum_{I, \mu} {p^{I \mu}_{\kappa'\beta}}^\ast ({\bf 0}) \, Z^{\ast}_{i,I \mu} \;. \]
    subroutine ph_write_borncharge_canonical
      use phonons_io_util, only: ph_io_write_borncharge
      use phonons_util, only: ph_util_sumrule_borncharge
      use modinput

      integer :: iq
      real(dp) :: correction(3, 3)
      logical :: success

      real(dp), allocatable :: zstar(:,:,:)

      if( mpiglobal%rank /= 0 ) return

      do iq = 1, ph_qset%nkpt
        if( norm2( ph_qset%vkl(:, iq) ) > 1e-6_dp ) cycle
        ! get Born effective charge in canonical coordinates
        call ph_borncharge_canonical_from_file( zstar )
        if( input%phonons%sumrule ) then
          ! impose acoustic sum rule
          call ph_util_sumrule_borncharge( zstar, correction )
          ! try to write Born effective charge to file
          call ph_io_write_borncharge( zstar, 'ZSTAR.OUT', success, &
            sumrule_correction=correction )
        else
          ! try to write Born effective charge to file
          call ph_io_write_borncharge( zstar, 'ZSTAR.OUT', success )
        end if
        if( success ) then
          write( *, * )
          write( *, '("Info (ph_write_borncharge_canonical):")' )
          write( *, '(" Born-effective charge tensors written to ZSTAR.OUT")' )
        end if
        call terminate_if_false( success, '(ph_write_borncharge_canonical) &
          Was not able to write Born effective charge tensors to file '//&
          new_line( 'a' )//'ZSTAR.OUT' )
      end do

      if( allocated( zstar ) ) deallocate( zstar )
    end subroutine ph_write_borncharge_canonical

    !> For a given phonon wavevector \({\bf q}\), this subroutine tries to read the dynamical matrix
    !> in irrep coordinates from a previous SCF run and transforms it into the canonical (Cartesian) basis.
    !>
    !> First, for the given \({\bf q}\) its symmetry equivalent point \({\bf q}_0\) in the set of
    !> \({\bf q}\)-vectors and the connecting symmetry operation is found. Then, the dynamical matrix
    !> in irrep coordinates at \({\bf q}_0\) is read from file, transformed into canonical coordinates via
    !> \[ D_{\kappa\alpha,\kappa'\beta}({\bf q}_0) = 
    !>    \sum_{I, \mu} {p^{I \mu}_{\kappa'\beta}}^\ast ({\bf q}_0) \, D_{\kappa\alpha,I \mu}({\bf q}_0) \;, \]
    !> where \(p^{I \mu}_{\kappa\alpha}({\bf q})\) are the irrep displacement patterns,
    !> and possibly rotated into \({\bf q}\) using [[ph_util_symapp_dyn(subroutine)]].
    subroutine ph_borncharge_canonical_from_file( zstar )
      use constants, only: zzero
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      !> dynamical matrix
      real(dp), allocatable, intent(out) :: zstar(:,:,:)

      integer :: iq, isym, iirrep, dirrep, id, &
                 is, ia, ias, ip
      complex(dp) :: pat(3, natmtot), zstar_i(3)
      character(:), allocatable :: qidirname, fxt
      logical :: success
      
      ! allocate output arrays
      if( allocated( zstar ) ) deallocate( zstar )
      allocate( zstar(3, 3, natmtot), source=0.0_dp )

      ! find equivalent Gamma point in set
      call findkptinset( [0._dp, 0._dp, 0._dp], ph_qset, isym, iq )

      ! read Born effective charge in irrep basis
      ! and transform it into canonical basis
      do iirrep = 1, ph_irrep_basis(iq)%nirrep
        dirrep = ph_irrep_basis(iq)%irreps(iirrep)%dim
        ! get subdirectory of part
        call ph_io_irrep_fxt( iq, iirrep, 0, qidirname )
        do id = 1, dirrep
          ! get file extension of part
          call ph_io_irrep_fxt( iq, iirrep, id, fxt )
          ! try to read Born effective charge from file
          call ph_io_read_borncharge_col( zstar_i, fxt, success, directory=qidirname )
          call terminate_if_false( success, '(ph_borncharge_canonical_from_file) &
            Was not able to read Born effective charge column supposed to be in file '//&
            new_line( 'a' )//trim( qidirname )//'/ZSTAR_'//trim(fxt)//'.OUT' )
          ! store irrep displacement pattern
          pat = ph_irrep_basis(iq)%irreps(iirrep)%pat(:, :, id)
          ! add contribution to dynamical matrix in canonical basis
          do is = 1, nspecies
            do ia = 1, natoms(is)
              ias = idxas(ia, is)
              do ip = 1, 3
                zstar(:, ip, ias) = zstar(:, ip, ias) + dble( conjg( pat(ip, ias) ) * zstar_i )
              end do
            end do
          end do
        end do
      end do
    end subroutine ph_borncharge_canonical_from_file
end module phonons
