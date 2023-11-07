!> This is the main moudle for calculating the response to an
!> external electric field using density-functional perturbation theory (DFPT).
module efield
  use dfpt_variables
  use dfpt_inout
  use efield_variables

  use modmpi
  use precision, only: dp
  use block_data_file, only: block_data_file_type

  implicit none
  private

  ! LOCAL VARIABLES
  !> eigenvalues at all \({\bf k}\) point
  real(dp), allocatable :: evalk(:,:)
  !> occupation numbers at all \({\bf k}\) point
  real(dp), allocatable :: occk(:,:)
  !> eigenvectors at single \({\bf k}\) point
  complex(dp), allocatable :: eveck(:,:)
  !> (L)APW matching coefficients \(A^\alpha_{{\bf G+p},lm,\xi}\) at \({\bf k}\)
  complex(dp), allocatable :: apwalmk(:,:,:,:)
  !> eigenvalue and occupation response at single \({\bf k}\) point
  real(dp), allocatable :: devalk(:,:,:), docck(:,:,:)
  !> eigenvector response at single \({\bf k}\) point
  complex(dp), allocatable :: deveck(:,:)
  !> constant part of Hamiltonian matrix response at all \({\bf k}\) points and polarization directions
  complex(dp), allocatable :: dHmat_const(:,:,:,:)
  !> full Hamiltonian matrix response at single \({\bf k}\) point and polarization direction
  complex(dp), allocatable :: dHmat(:,:)
  !> density response matrix for all atoms and polarization directions
  complex(dp), allocatable :: drho_mat(:,:,:,:)
  !> muffin-tin and interstitial density response for all polarization directions
  complex(dp), allocatable :: drho_mt(:,:,:,:), drho_ir(:,:)
  !> muffin-tin and interstitial effective potential response for all polarization directions
  complex(dp), allocatable :: dpot_mt(:,:,:,:), dpot_ir(:,:)
  !> radial integrals of effective potential response times Gaunt coefficients
  complex(dp), allocatable :: dHmat_mt_basis(:,:,:,:)
  !> interstitial potential response times characteristic function in reciprocal space
  complex(dp), allocatable :: dpot_cfun_ig(:,:)
  !> interstitial (scalar relativistic) kinetic energy response times characteristic function in reciprocal space
  complex(dp), allocatable :: dkin_cfun_ig(:,:)
  !> polarization response for all polarization directions
  complex(dp), allocatable :: dpol(:,:)
  !> parallel i/o file objects for eigenvalue, occupation and eigenvector response on \({\bf k}\)-grid
  type(block_data_file_type) :: fdevalk, fdocck, fdeveck
  !> `.true.`, if this process is the local master of the current part
  logical :: master

  public :: ef_prepare, ef_finalize
  public :: ef_scf, ef_polarization

  contains
    
    !> This subroutine executes preparative tasks for a DFPT electric field calculation.
    !>
    !> This includes:
    !>
    !> * initialization of global electric field variables
    !> * obtaining eigenvalues and occupation numbers on the \({\bf k}_0\) points
    !>   of the unperturbed system
    !> * calculation of the constant parts for density and potential response
    subroutine ef_prepare( &
        info_output )
      use dfpt_eigensystem, only: dfpt_eig_ks
      use efield_eigensystem, only: ef_eig_init, ef_eig_gen_dHmat

      use constants, only: zzero, zone
      use mod_APW_LO, only: nlotot, apwordmax
      use mod_atoms, only: natmtot
      use mod_muffin_tin, only: lmmaxapw, nrmtmax
      use mod_eigenvalue_occupancy, only: nstfv, occmax, efermi
      use mod_charge_and_moment, only: chgval
      use mod_Gkvector, only: ngkmax_ptr
      use modinput
      !> create and write info output files (default: `.true.`)
      logical, optional, intent(in) :: info_output

      integer :: ik, ik1, ik2, ip, nmatmax
      integer, target :: ngkmax
      logical :: write_info, exists, from_file, success
      character(1) :: fxt

      write_info = .true.
      if( present( info_output ) ) write_info = info_output

      master = (mpiglobal%rank == 0)
      ngkmax = dfpt_Gkset%ngkmax
      nmatmax = ngkmax + nlotot
      
      ! check if temporary file for eigenvectors and eigenvalues are available
      exists = fevalk0%exists() .and. feveck0%exists()
      exists = exists .and. (fevalk0%is_open() .and. feveck0%is_open())
      call terminate_if_false( exists, '(ef_prepare) &
             Eigenvalues and eigenvectors were not prepared. Call `dfpt_prepare` in advance.' )
      ! initialize info output
      if( master .and. write_info ) call dfpt_io_info_init
      ! set limits for k-point loops
      ik1 = firstofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
      ik2 = lastofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
      ! allocate local variables
      allocate( evalk(nmatmax, dfpt_kset%nkpt), devalk(nstfv, dfpt_kset%nkpt, 3) )
      allocate( occk(nmatmax, dfpt_kset%nkpt), docck(nstfv, dfpt_kset%nkpt, 3) )
      allocate( eveck(nmatmax, nmatmax), deveck(nmatmax, nstfv) )
      allocate( apwalmk(ngkmax, apwordmax, lmmaxapw, natmtot) )
      allocate( dHmat_const(nmatmax, nstfv, ik1:ik2, 3) )
      allocate( dHmat(nmatmax, nstfv) )
      allocate( dHmat_mt_basis(mt_basis%n_basis_fun_max, mt_basis%n_basis_fun_max, natmtot, 3) )
      allocate( dpot_cfun_ig(dfpt_Gset%ngvec, 3) )
      allocate( dkin_cfun_ig(dfpt_Gset%ngvec, 3) )
      allocate( drho_mat(mt_basis%n_basis_fun_max, mt_basis%n_basis_fun_max, natmtot, 3) )
      allocate( drho_mt(dfpt_lmmaxvr, nrmtmax, natmtot, 3) )
      allocate( drho_ir(dfpt_Gset%ngrtot, 3) )
      allocate( dpot_mt(dfpt_lmmaxvr, nrmtmax, natmtot, 3) )
      allocate( dpot_ir(dfpt_Gset%ngrtot, 3) )
      allocate( dpol(3, 3) )
      ! loop over k-points
      do ik = ik1, ik2
        ! solve KS equation at k
        call dfpt_eig_ks( ik, dfpt_kset, dfpt_Gset, dfpt_Gkset, nmatmax, evalk(:, ik), eveck, &
               p0set=dfpt_kset, Gp0set=dfpt_Gkset, feval=fevalk0, fevec=feveck0 )
      end do
      call mpi_allgatherv_ifc( dfpt_kset%nkpt, rlen=nmatmax, comm=mpiglobal, rbuf=evalk )
      ! get occupation numbers
      occk = 0.0_dp
      call find_fermi( dfpt_kset%nkpt, dfpt_kset%wkpt, nstfv, evalk(1:nstfv, :), chgval, occmax, &
             input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
             efermi, occk(1:nstfv, :) )
      ! initialize eigensystem response
      call ef_eig_init
      ! compute constant part of Hamiltonian response
      ! for all k points and polarization directions
      dHmat_const = zzero
      do ik = ik1, ik2
        ! get matching coefficients at k
        ngkmax_ptr => ngkmax
        call match( dfpt_Gkset%ngk(1, ik), dfpt_Gkset%gkc(:, 1, ik), dfpt_Gkset%tpgkc(:, :, 1, ik), dfpt_Gkset%sfacgk(:, :, 1, ik), &
                    apwalmk )
        ! read eigenvectors
        call feveck0%read( ik, eveck )
        do ip = 1, 3
          call ef_eig_gen_dHmat( ik, dfpt_Gkset, 1, nstfv, evalk(:, ik), eveck, apwalmk, dHmat_const(:, :, ik, ip), &
            ip=ip )
        end do
      end do
      ! initialize density and potential response
      from_file = (input%phonons%do == 'fromfile')
      if( master ) then
        docck = 0.0_dp
        ip = 0
        do while( ip < 3 )
          ip = ip + 1
          if( from_file ) then
            write( fxt, '(i1.1)' ) ip
            call dfpt_io_read_zfun( drho_mt(:, :, :, ip), drho_ir(:, ip), 'EFIELD_DRHO', success, &
              file_extension=fxt )
            from_file = from_file .and. success
            call dfpt_io_read_zfun( dpot_mt(:, :, :, ip), dpot_ir(:, ip), 'EFIELD_DVEFF', success, &
              file_extension=fxt )
            from_file = from_file .and. success
            if( .not. from_file ) ip = 0
          else
            drho_mt(:, :, :, ip) = zzero; drho_ir(:, ip) = zzero
            dpot_mt(:, :, :, ip) = zzero; dpot_ir(:, ip) = zzero
          end if
        end do
        if( from_file ) call dfpt_io_info_string( 'Density and potential response read from file.' )
      end if
      ! open files for eigenvalue, occupancy and eigenvector response
      fdevalk = block_data_file_type( 'EFIELD_DEVAL.OUT', [nstfv], 1.0_dp )
      fdocck = block_data_file_type( 'EFIELD_DOCC.OUT', [nstfv], 1.0_dp )
      fdeveck = block_data_file_type( 'EFIELD_DEVEC.OUT', [nmatmax, nstfv], zone )
      call fdevalk%open( mpiglobal )
      call fdocck%open( mpiglobal )
      call fdeveck%open( mpiglobal )

      call barrier( mpicom=mpiglobal )
    end subroutine ef_prepare

    !> This subroutine executes finalizing tasks for a DFPT electric field calculation.
    !>
    !> This includes freeing memory from unneeded variables, closing files and 
    !> deleting temporary files.
    subroutine ef_finalize( &
        info_output, delete_files )
      use efield_eigensystem, only: ef_eig_free
      !> create and write info output files (default: `.true.`)
      logical, optional, intent(in) :: info_output
      !> delete files for eigensystem response (default: `.true.`)
      logical, optional, intent(in) :: delete_files

      logical :: write_info, delete

      write_info = .true.
      if( present( info_output ) ) write_info = info_output
      delete = .true.
      if( present( delete_files ) ) delete = delete_files

      ! close file for eigensystem response
      call fdevalk%close( mpiglobal )
      call fdocck%close( mpiglobal )
      call fdeveck%close( mpiglobal )
      ! delete file for eigensystem response
      if( delete ) then
        call fdevalk%delete( mpiglobal )
        call fdocck%delete( mpiglobal )
        call fdeveck%delete( mpiglobal )
      end if
      ! free unneeded variables
      if( allocated( evalk ) ) deallocate( evalk )
      if( allocated( occk ) ) deallocate( occk )
      if( allocated( eveck ) ) deallocate( eveck )
      if( allocated( apwalmk ) ) deallocate( apwalmk )
      if( allocated( devalk ) ) deallocate( devalk )
      if( allocated( docck ) ) deallocate( docck )
      if( allocated( deveck ) ) deallocate( deveck )
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
      if( allocated( dpol ) ) deallocate( dpol )
      call ef_eig_free
      ! finalize info output
      if( master .and. write_info ) call dfpt_io_info_finit
    end subroutine ef_finalize

    !> This subroutine runs the self-consistency cycle for obataining 
    !> the density and potential response for a DFPT electric field calculation.
    !>
    !> Each SCF iteration consists of:
    !>
    !> * for each \({\bf k}\) point and polarization direction calculate the
    !>   Hamiltonian matrix response
    !> * for each \({\bf k}\) point and polarization direction solve the Sternheimer
    !>   equation for the eigenvalue and eigenvector response
    !> * for each polarization direction synthesize the new density response
    !>   and compute the corresponding potential response
    !> * mix old and new density/potential response
    subroutine ef_scf( &
        info_output )
      use dfpt_density_potential, only: dfpt_rhopot_mixpack, dfpt_rhopot_drho_k, dfpt_rhopot_gen_drho_mt
      use dfpt_eigensystem, only: dfpt_eig_prepare_dHmat
      use efield_density_potential, only: ef_rhopot_gen_dpot, ef_rhopot_symmetrize
      use efield_eigensystem, only: ef_eig_gen_dHmat, ef_eig_sternheimer

      use constants, only: zzero
      use mod_potential_and_density, only: pot_mt => veffmt, pot_ir => veffir
      use mod_convergence, only: iscl
      use mod_eigenvalue_occupancy, only: nstfv, occmax, efermi
      use mod_Gkvector, only: ngkmax_ptr
      use mod_symmetry, only: nsymcrys
      use modinput
      !> write info output files (default: `.true.`)
      logical, optional, intent(in) :: info_output

      integer :: ik, ik1, ik2, ip, i
      integer :: nmix, mixermode
      integer, target :: ngkmax
      real(dp) :: conv, t_iter, t1, defermi
      logical :: write_info, success, lastiter, converged, skip, to_file
      character(1) :: fxt

      real(dp), allocatable :: rvmix(:), vconv(:)

      write_info = .true.
      if( present( info_output ) ) write_info = info_output

      ! set limits for k-point loops
      ik1 = firstofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
      ik2 = lastofset( mpiglobal%rank, dfpt_kset%nkpt, mpiglobal%procs )
      ! initialize mixer
      iscl = 0
      if( master ) then
        nmix = 2 * (size( dpot_mt ) + size( dpot_ir ))
        mixermode = -1
        allocate( rvmix(nmix) )
        allocate( vconv(input%groundstate%niterconvcheck) )
        call dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, .true., 3, nmix, rvmix )
        call mixerifc( input%groundstate%mixernumber, nmix, rvmix, conv, mixermode )
        conv = 1.0_dp
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
#ifdef MPI
        call MPI_Bcast( dpot_mt, size( dpot_mt ), MPI_DOUBLE_COMPLEX, 0, mpiglobal%comm, mpiglobal%ierr )
        call MPI_Bcast( dpot_ir, size( dpot_ir ), MPI_DOUBLE_COMPLEX, 0, mpiglobal%comm, mpiglobal%ierr )
#endif
        ! send occupation response to other processes
#ifdef MPI
        call MPI_Bcast( docck, size( docck ), MPI_DOUBLE, 0, mpiglobal%comm, mpiglobal%ierr )
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
          ngkmax = dfpt_Gkset%ngkmax
          do ik = ik1, ik2
            ! get matching coefficients at k
            ngkmax_ptr => ngkmax
            call match( dfpt_Gkset%ngk(1, ik), dfpt_Gkset%gkc(:, 1, ik), dfpt_Gkset%tpgkc(:, :, 1, ik), dfpt_Gkset%sfacgk(:, :, 1, ik), &
                        apwalmk )
            ! read eigenvectors
            call feveck0%read( ik, eveck )
            do ip = 1, 3
              ! generate k-independent part of Hamiltonian response
              if( ik == ik1 ) &
                call dfpt_eig_prepare_dHmat( pot_mt, pot_ir, dpot_mt(:, :, :, ip), dpot_ir(:, ip), dHmat_mt_basis(:, :, :, ip), dpot_cfun_ig(:, ip), dkin_cfun_ig(:, ip) )
              ! generate full Hamiltonian response
              dHmat = dHmat_const(:, :, ik, ip)
              call ef_eig_gen_dHmat( ik, dfpt_Gkset, 1, nstfv, &
                evalk(:, ik), eveck, apwalmk, dHmat, &
                dHmat_mt_basis=dHmat_mt_basis(:, :, :, ip), dpot_cfun_ig=dpot_cfun_ig(:, ip), dkin_cfun_ig=dkin_cfun_ig(:, ip) )
              ! solve Sternheimer equation
              call ef_eig_sternheimer( ik, dfpt_Gkset, 1, nstfv, &
                evalk(:, ik), occk(:, ik), eveck, dHmat, &
                devalk(:, ik, ip), deveck, projector=.true. )
              ! write eigenvector response to file
              call fdeveck%write( (ik-1)*3+ip, deveck )
              ! add k-point contribution to density response
              call dfpt_rhopot_drho_k( ik, dfpt_kset, dfpt_Gkset, dfpt_Gkset, 1, nstfv, &
                occk(:, ik), eveck, deveck, apwalmk, apwalmk, &
                drho_mat(:, :, :, ip), drho_ir(:, ip), &
                docck=docck(:, ik, ip) )
            end do
          end do
          call barrier( mpicom=mpiglobal )
          ! reduce density response over k-points
#ifdef MPI
          call MPI_Allreduce( MPI_IN_PLACE, drho_mat, size( drho_mat ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                 mpiglobal%comm, mpiglobal%ierr )
          call MPI_Allreduce( MPI_IN_PLACE, drho_ir, size( drho_ir ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
                 mpiglobal%comm, mpiglobal%ierr )
#endif
          ! update muffin-tin density response
          do ip = 1, 3
            call dfpt_rhopot_gen_drho_mt( drho_mat(:, :, :, ip), drho_mt(:, :, :, ip) )
          end do
          ! reduce eigenvalue response over k-points
#ifdef MPI
          call MPI_Allreduce( MPI_IN_PLACE, devalk, size( devalk ), MPI_DOUBLE, MPI_SUM, &
                 mpiglobal%comm, mpiglobal%ierr )
#endif
          if( master ) then
            ! symmetrize density response
            call ef_rhopot_symmetrize( drho_mt, drho_ir, nsymcrys, [(i, i=1, nsymcrys)] )
            ! update effective potential response
            do ip = 1, 3
              call ef_rhopot_gen_dpot( drho_mt(:, :, :, ip), drho_ir(:, ip), dpot_mt(:, :, :, ip), dpot_ir(:, ip) )
            end do
            ! symmetrize potential response
            call ef_rhopot_symmetrize( dpot_mt, dpot_ir, nsymcrys, [(i, i=1, nsymcrys)] )
            ! uppdate occupation response
            ! and write eigenvalue and occupation response to file
            do ip = 1, 3
              call find_dfermi( dfpt_kset%nkpt, dfpt_kset%wkpt, nstfv, evalk(1:nstfv, :), devalk(:, :, ip), 0.0_dp, occmax, efermi, &
                     input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
                     defermi, docck(:, :, ip) )
              do ik = 1, dfpt_kset%nkpt
                call fdevalk%write( (ik-1)*3+ip, devalk(:, ik, ip) )
                call fdocck%write( (ik-1)*3+ip, docck(:, ik, ip) )
              end do
            end do
          end if
        end if ! end of scf step

        ! write density response to file
        if( master ) then
          to_file = .false.
          if( input%groundstate%nwrite > 0 ) to_file = (mod( max( 1, iscl - 1 ), input%groundstate%nwrite ) == 0)
          to_file = to_file .or. lastiter
          if( to_file ) then
            do ip = 1, 3
              write( fxt, '(i1.1)' ) ip
              call dfpt_io_write_zfun( drho_mt(:, :, :, ip), drho_ir(:, ip), 'EFIELD_DRHO', success, &
                file_extension=fxt )
              to_file = to_file .and. success
              call dfpt_io_write_zfun( dpot_mt(:, :, :, ip), dpot_ir(:, ip), 'EFIELD_DVEFF', success, &
                file_extension=fxt )
              to_file = to_file .and. success
            end do
          end if
          if( to_file ) call dfpt_io_info_string( 'Density and potential response written to file.' )
        end if

        ! exit on demand
        if( skip .and. lastiter ) exit scfloop

        ! mixing
        if( master ) then
          call dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, .true., 3, nmix, rvmix )
          call mixerifc( input%groundstate%mixernumber, nmix, rvmix, conv, mixermode )
          do i = 1, input%groundstate%niterconvcheck - 1
            vconv(i) = vconv(i+1)
          end do
          vconv(i) = conv
          ! check convergence
          converged = all( abs( vconv ) < input%groundstate%epspot )
          if( .not. lastiter ) &
            call dfpt_rhopot_mixpack( drho_mt, drho_ir, dpot_mt, dpot_ir, .false., 3, nmix, rvmix )
        end if
#ifdef MPI
        call MPI_Bcast( lastiter, 1, MPI_LOGICAL, 0, mpiglobal%comm, mpiglobal%ierr )
        call MPI_Bcast( converged, 1, MPI_LOGICAL, 0, mpiglobal%comm, mpiglobal%ierr )
#endif
        call timesec( t_iter )
        t_iter = t_iter - t1

        ! write iteration info to output file
        if( master .and. write_info ) then
          call dfpt_io_info_scf( iscl, conv, t_iter, defermi=defermi )
        end if

        call barrier( mpicom=mpiglobal )
      end do scfloop
      !* END SCF LOOP
      !********************************

      ! finalize scf loop output
      if( master .and. write_info ) call dfpt_io_info_scf_finit

      ! deallocate mixer
      if( master ) then
        mixermode = -2
        call mixerifc( input%groundstate%mixernumber, nmix, rvmix, conv, mixermode )
        deallocate( rvmix, vconv )
      end if

      call barrier( mpicom=mpiglobal )
    end subroutine ef_scf

    !> This subroutine calculates the polarization response (= Born effective charge column) 
    !> from the converged density and potential response for an independent part
    !> of a DFPT phonon calculation and writes it to file.
    !>
    !> The polarization response is calculated using the Berry phase approach.
    !>
    !> See [[ph_borncharge_canonical_from_file(subroutine)]] for the transformation of the Born effective 
    !> charge into canonical (Cartesian) coordinates.
    subroutine ef_polarization
      use dfpt_polarization
      use efield_polarization
      use efield_eigensystem, only: ef_eig_rotate_devec
      use phonons_io_util, only: ph_io_write_dielten

      use constants, only: zzero, zone, twopi
      use mod_eigenvalue_occupancy, only: occmax
      use mod_symmetry, only: nsymcrys

      integer :: ik, ikb, ik1, ik2, ip1, ip2, i
      integer :: iknr, ikbnr, isym, nst, ivkb(3)
      real(dp) :: binv(3, 3), tmp(3, 2), t1
      logical :: success, to_file

      complex(dp), allocatable :: devecknr(:,:,:), eveckbnr(:,:), deveckbnr(:,:,:)

      real(dp), external :: r3mdet

      ! set limits for k-point loops
      ik1 = firstofset( mpiglobal%rank, dfpt_kset%nkptnr, nprocs=mpiglobal%procs )
      ik2 = lastofset( mpiglobal%rank, dfpt_kset%nkptnr, nprocs=mpiglobal%procs )

      allocate( devecknr(size( deveck, dim=1 ), size( deveck, dim=2 ), 3) )
      allocate( eveckbnr(size( eveck, dim=1 ), size( eveck, dim=2 )) )
      allocate( deveckbnr(size( deveck, dim=1 ), size( deveck, dim=2 ), 3) )

      dpol = zzero
      do ip1 = 1, 3
        ! initialize plane wave matrix elements
        call dfpt_pol_init_p( - dfpt_kset%bvec(:, ip1) / dfpt_kset%ngridk(ip1) )
        ! we have to loop over the full k-set
        do iknr = ik1, ik2
          ! find equivalent point in reduced set and corresponding symmetry
          call findkptinset( dfpt_kset%vklnr(:, iknr), dfpt_kset, isym, ik )
          call terminate_if_false( ik > 0, '(ef_polarization) &
            Equivalent k-point not found.' )
          ! get number of occupied states
          nst = count( occk(:, ik) > occmax/2 )
          ! read eigenvector at k
          call feveck0%read( ik, eveck )
          call rotate_evecfv( isym, dfpt_kset%vkl(:, ik), dfpt_kset%vklnr(:, iknr), &
            dfpt_Gkset%ngk(1, ik), dfpt_Gkset%vgkl(:, :, 1, ik), dfpt_Gkset%vgknrl(:, :, 1, iknr), &
            eveck, size( eveck, dim=1 ), nst )
          ! read eigenvector response at k
          do ip2 = 1, 3
            call fdeveck%read( (ik-1)*3+ip2, devecknr(:, :, ip2) )
          end do
          call ef_eig_rotate_devec( isym, dfpt_kset%vkl(:, ik), dfpt_kset%vklnr(:, iknr), &
            dfpt_Gkset%ngk(1, ik), dfpt_Gkset%vgkl(:, :, 1, ik), dfpt_Gkset%vgknrl(:, :, 1, iknr), &
            devecknr, size( devecknr, dim=1 ), size( devecknr, dim=2 ), nst )
          ! find vector k+b and its equivalent in reduced set
          ivkb = dfpt_kset%ivknr(:, iknr)
          ivkb(ip1) = modulo( ivkb(ip1) + 1, dfpt_kset%ngridk(ip1) )
          ikbnr = dfpt_kset%ikmapnr(ivkb(1), ivkb(2), ivkb(3))
          call findkptinset( dfpt_kset%vklnr(:, ikbnr), dfpt_kset, isym, ikb )
          call terminate_if_false( ikb > 0, '(ef_polarization) &
            Equivalent k-point not found.' )
          ! read eigenvector at k+b
          call feveck0%read( ikb, eveckbnr )
          call rotate_evecfv( isym, dfpt_kset%vkl(:, ikb), dfpt_kset%vklnr(:, ikbnr), &
            dfpt_Gkset%ngk(1, ikb), dfpt_Gkset%vgkl(:, :, 1, ikb), dfpt_Gkset%vgknrl(:, :, 1, ikbnr), &
            eveckbnr, size( eveckbnr, dim=1 ), nst )
          ! read eigenvector response at k+b
          do ip2 = 1, 3
            call fdeveck%read( (ikb-1)*3+ip2, deveckbnr(:, :, ip2) )
          end do
          call ef_eig_rotate_devec( isym, dfpt_kset%vkl(:, ikb), dfpt_kset%vklnr(:, ikbnr), &
            dfpt_Gkset%ngk(1, ikb), dfpt_Gkset%vgkl(:, :, 1, ikb), dfpt_Gkset%vgknrl(:, :, 1, ikbnr), &
            deveckbnr, size( deveckbnr, dim=1 ), size( deveckbnr, dim=2 ), nst )

          ! add k-point contribution to polarization response
          do ip2 = 1, 3
            call ef_pol_dpol_k( iknr, ikbnr, dfpt_kset, dfpt_Gkset, 1, nst, &
              eveck, eveckbnr, devecknr(:, :, ip2), deveckbnr(:, :, ip2), &
              dpol(ip1, ip2) )
          end do

        end do
      end do
      ! free memory
      call dfpt_pol_free
      call barrier( mpicom=mpiglobal )
      ! reduce polarization response over k-points
#ifdef MPI
      call MPI_Allreduce( MPI_IN_PLACE, dpol, size( dpol ), MPI_DOUBLE_COMPLEX, MPI_SUM, &
             mpiglobal%comm, mpiglobal%ierr )
#endif
      ! transform into Cartesian coordinates
      call r3minv( dfpt_kset%bvec, binv )
      t1 = - abs( r3mdet( dfpt_kset%bvec ) ) / twopi**2 * 2 ! = 4pi / omega
      do ip2 = 1, 3
        dpol(:, ip2) = t1 * dpol(:, ip2) * dfpt_kset%ngridk
        call r3mtv( binv, dble(  dpol(:, ip2) ), tmp(:, 1) )
        call r3mtv( binv, aimag( dpol(:, ip2) ), tmp(:, 2) )
        dpol(:, ip2) = cmplx( tmp(:, 1), tmp(:, 2), dp )
      end do
      ! add ionic contribution
      do ip1 = 1, 3
        dpol(ip1, ip1) = dpol(ip1, ip1) + zone
      end do
      ! symmetrize
      call ef_pol_symmetrize( dpol, nsymcrys, [(i, i=1, nsymcrys)] )
      ! write to file
      if( master ) then
        to_file = .true.
        call ph_io_write_dielten( dble( dpol ), 'EPSINF.OUT', success )
        if( to_file ) then
          call dfpt_io_info_string( 'Dielectric tensor written to file.' )
          write( *, * )
          write( *, '("Info (ef_polarization):")' )
          write( *, '(" dielectric tensor written to EPSINF.OUT")' )
        end if
      end if

      deallocate( devecknr, eveckbnr, deveckbnr )

      call barrier( mpicom=mpiglobal )
    end subroutine ef_polarization
end module efield
