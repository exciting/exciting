!> Module handling phonon part of an EPH calculation.
module eph_phonons
  use eph_variables

  use precision, only: dp
  use asserts, only: assert
  use modmpi
  use matrix_fourier_interpolation, only: mfi_type
  use block_data_file, only: block_data_file_type

  implicit none
  private

  !> phonon frequencies on phonon \({\bf q}\)-grid
  real(dp), allocatable, public :: eph_ph_energy_q(:,:)
  !> phonon eigenvectors \(e_{\kappa\alpha,\nu}({\bf q})\) on phonon \({\bf q}\)-grid
  complex(dp), allocatable, public :: eph_ph_evec_q(:,:,:)
  !> object for matrix Fourier interpolation on phonon \({\bf q}\)-grid
  type(mfi_type), public :: eph_ph_mfi
  !> phonon dynamical matrices in real space atomic gauge, \(\mathcal{D}_{\kappa\alpha,\lambda\beta}({\bf R})\), 
  !> (interatomic force constants IFCs)
  complex(dp), allocatable :: eph_ph_DR(:,:,:)
  !> name for binary file to save dynamical matrices in real space atomic gauge for later access
  character(*), parameter :: eph_ph_DR_filename = "EPH_DR.OUT"

  public :: eph_ph_free, eph_ph_setup_interpolation, eph_ph_interpolate

contains

  !> Free memory from module variables.
  subroutine eph_ph_free
    if (allocated(eph_ph_energy_q)) deallocate( eph_ph_energy_q )
    if (allocated(eph_ph_evec_q)) deallocate( eph_ph_evec_q )
    if (allocated(eph_ph_DR)) deallocate( eph_ph_DR )
    call eph_ph_mfi%destroy
  end subroutine eph_ph_free

  !================================================================================ 
  ! SETUP FOURIER INTERPOLATION OF PHONONS
  !
  !> Set up the Fourier interpolation of the dynamical matrix.
  !>
  !> This includes
  !> 
  !>   * checking for availability of phonons
  !>   * preparation of matrix Fourier interpolation on phonon \({\bf q}\)-grid
  !>   * computation of localized dynamical matrices in real space atomic gauge \(\mathcal{D}_{\kappa\alpha,\lambda\beta}({\bf R})\)
  !>     and writing to file
  !>   * or reading \(\mathcal{D}_{\kappa\alpha,\lambda\alpha}({\bf R})\) from file, if possible
  subroutine eph_ph_setup_interpolation( write_localization )
    use phonons_util, only: ph_util_setup_interpolation, ph_util_interpolate, ph_util_diag_dynmat
    use mod_atoms, only: nspecies, natoms, natmtot, atposc
    use modinput
    !> write spatial localization of \(\mathcal{\bf D}({\bf R})\) to file (default: `.false.`)
    logical, optional, intent(in) :: write_localization

    integer :: iq, ir, ia, is, un, stat
    logical :: write_loc
    type(block_data_file_type) :: DR_file

    complex(dp), allocatable :: Dq(:,:,:)

    write_loc = .false.
    if (present(write_localization)) write_loc = write_localization

    ! delete existing module variables
    call eph_ph_free

    ! set up binary file
    DR_file = block_data_file_type( eph_ph_DR_filename, [3*natmtot, 3*natmtot], cmplx( 0, 0, dp ) )

    ! check if there is already a file with D_a(R)
    if (DR_file%exists()) then
      ! setup matrix fourier interpolation
      eph_ph_mfi = mfi_type( eph_qset_ph%ngridk, eph_qset_ph%bvec, eph_qset_ph%vklnr )
      call eph_ph_mfi%set_localization_centers( &
        left_centers=reshape( [((atposc(:, ia, is), ia=1, natoms(is)), is=1, nspecies)], [3, natmtot] ), &
        right_centers=reshape( [((atposc(:, ia, is), ia=1, natoms(is)), is=1, nspecies)], [3, natmtot] ), &
        coordinates='c' )
      ! allocate real space atomic dynamical matrix D_a(R)
      allocate( eph_ph_DR(3*natmtot, 3*natmtot, eph_ph_mfi%nr) )
      ! try to read D_a(R) from file
      call DR_file%open( mpiglobal )
      do ir = 1, eph_ph_mfi%nr
        call DR_file%read( ir, eph_ph_DR(:, :, ir) )
      end do
      call DR_file%close( mpiglobal )
    ! otherwise compute D_a(R) and write to file
    else
      ! setup matrix fourier interpolation and real space atomic dynamical matrix D_a(R)
      if (eph_polar) then
        call ph_util_setup_interpolation( eph_qset_ph%bvec, eph_qset_ph%ngridk, eph_qset_ph%nkpt, eph_qset_ph%ivk, eph_qset_ph%vkl, &
          eph_ph_mfi, eph_ph_DR, &
          sumrule=input%phonons%sumrule, &
          dielten=eph_dielten, borncharge=eph_borncharge, &
          elphbolt_compatible=input%eph%elphbolt )
      else
        call ph_util_setup_interpolation( eph_qset_ph%bvec, eph_qset_ph%ngridk, eph_qset_ph%nkpt, eph_qset_ph%ivk, eph_qset_ph%vkl, &
          eph_ph_mfi, eph_ph_DR, &
          sumrule=input%phonons%sumrule, &
          elphbolt_compatible=input%eph%elphbolt )
      end if
      ! write D_a(R) to file
      call DR_file%open( mpiglobal )
      do ir = 1, eph_ph_mfi%nr
        call DR_file%write( ir, eph_ph_DR(:, :, ir) )
      end do
      call DR_file%close( mpiglobal )
    end if
    
    ! get phonon frequencies and eigenvectors on support grid
    allocate( eph_ph_energy_q(eph_nmode_tot, eph_qset_ph%nkpt) )
    allocate( eph_ph_evec_q(3*natmtot, eph_nmode_tot, eph_qset_ph%nkpt) )
    if (eph_polar) then
      call ph_util_interpolate( eph_qset_ph%nkpt, eph_qset_ph%vkl(:, 1:eph_qset_ph%nkpt), eph_ph_mfi, eph_ph_DR, Dq, & 
        dielten=eph_dielten, borncharge=eph_borncharge, &
        elphbolt_compatible=input%eph%elphbolt )
    else
      call ph_util_interpolate( eph_qset_ph%nkpt, eph_qset_ph%vkl(:, 1:eph_qset_ph%nkpt), eph_ph_mfi, eph_ph_DR, Dq, &
        elphbolt_compatible=input%eph%elphbolt ) 
    end if
    do iq = 1, eph_qset_ph%nkpt
      call ph_util_diag_dynmat( Dq(:, :, iq), eph_ph_energy_q(:, iq), eph_ph_evec_q(:, :, iq) )
    end do
    if (allocated(Dq)) deallocate( Dq )

    ! write spatial localization to file
    if (mpiglobal%rank == 0 .and. write_loc) then
      open( newunit=un, file='eph_ph_loc.dat', action='write', form='formatted', iostat=stat )
      call terminate_if_false( stat == 0, '(eph_ph_setup_interpolation) &
        Failed to open file `eph_ph_loc.dat`.' )
      write( un, '("#",a5,a26,a26)' ) 'iR', '|R|', 'max(|D(R)|)'
      do ir = 1, eph_ph_mfi%nr
        write( un, '(i6,2g26.16)' ) ir, eph_ph_mfi%rlen(ir), maxval( abs( eph_ph_DR(:, :, ir) ) )
      end do
      close( un )
    end if
  end subroutine eph_ph_setup_interpolation
  !-------------------------------------------------------------------------------- 

  !================================================================================ 
  ! INTERPOLATE PHONONS
  !
  !> Get phonon frequencies \(\omega_{\nu}({\bf q}')\) and eigenvectors \(e_{\kappa\alpha,\nu}({\bf q}')\)
  !> for a given set of wave vectors \({\bf q}'\) by Fourier interpolation.
  !>
  !> This is done by Fourier interpolating the dynamical matrix to \({\bf q}'\) by
  !> \[ \mathcal{D}({\bf q}') = \sum_{\bf R} {\rm e}^{{\rm i} {\bf q}'\cdot{\bf R}}\, \mathcal{D}({\bf R}) \]
  !> and diagonalizing \(\mathcal{D}({\bf q}')\) such that
  !> \[ D({\bf q}') = \operatorname{diag}\left(\omega^2_{\nu}({\bf q}')\right) 
  !>    = e^\dagger({\bf q}')\, M^{-1/2}\, \mathcal{D}({\bf q}')\, M^{-1/2}\, e({\bf q}') \;, \]
  !> where \(M^{-1/2}_{\kappa\alpha,\lambda\beta} = M_{\kappa}^{-1/2}\, \delta_{\kappa\lambda}\, \delta_{\alpha\beta}\)
  !> with \(M_\kappa\) being the mass of nucleus \(\kappa\).
  !> See also [[ph_util_interpolate(subroutine)]].
  !>
  !> MPI parallelization is over \({\bf q}'\) points, but each process can specify a different mode range `mrange`.
  subroutine eph_ph_interpolate( vql, phfreq, phevec, &
      mrange, mpicomm )
    use phonons_util, only: ph_util_interpolate, ph_util_diag_dynmat
    use exciting_mpi, only: xmpi_allreduce
#ifdef MPI
    use mpi_f08, only: MPI_Send, MPI_Recv, MPI_DOUBLE, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, MPI_COMM
#endif
    use mod_atoms, only: natmtot
    use modinput
    !> set of wave vectors \({\bf q}'\) in lattice coordinates
    real(dp), intent(in) :: vql(:,:)
    !> phonon frequencies \(\omega_{\nu}({\bf q}')\)
    real(dp), allocatable, intent(out) :: phfreq(:,:)
    !> phonon eigenvectors \(e_{\kappa\alpha,\nu}({\bf q}')\)
    complex(dp), allocatable, intent(out) :: phevec(:,:,:)
    !> range of phonon modes for which to interpolate (default: `[eph_fmode, eph_lmode]`)
    integer, optional, intent(in) :: mrange(2)
    !> MPI communicator (default: global MPI communicator)
    type(mpiinfo), optional, intent(inout) :: mpicomm

    integer :: nq, iq, iq1, iq2, nmode, mrng(2), irank, jrank
    type(mpiinfo) :: mpi

    integer, allocatable :: mrng_list(:,:), iq_range(:,:), ix_range(:,:)
    real(dp), allocatable :: eval(:,:)
    complex(dp), allocatable :: dynmat(:,:,:), evec(:,:)

    mrng = [eph_fmode, eph_lmode]
    if (present(mrange)) mrng = mrange

    mpi = mpiglobal
    if (present(mpicomm)) mpi = mpicomm

    nmode = mrng(2) - mrng(1) + 1
    nq = size( vql, dim=2 )

    call assert( size( vql, dim=1 ) == 3, &
      '`vql` must be a set of vectors of length 3.' )

    ! communicate information on mode distribution
    allocate( mrng_list(2, mpi%procs) )
    mrng_list = 0; mrng_list(:, mpi%rank+1) = mrng
    call xmpi_allreduce( mrng_list, mpi )

    ! distribute q among processes
    call patchwork_distribution( nq, 1, mpi%procs, iq_range, ix_range )
    iq1 = iq_range(1, mpi%rank+1)
    iq2 = iq_range(2, mpi%rank+1)

    ! allocate output
    if (allocated(phfreq)) deallocate( phfreq )
    allocate( phfreq(mrng(1):mrng(2), nq) )
    if (allocated(phevec)) deallocate( phevec )
    allocate( phevec(3*natmtot, mrng(1):mrng(2), nq) )

    ! interpolate dynamical matrix
    if (eph_polar) then
      call ph_util_interpolate( iq2-iq1+1, vql(:, iq1:iq2), eph_ph_mfi, eph_ph_DR, dynmat, & 
        dielten=eph_dielten, borncharge=eph_borncharge, minimal_distances=.true., &
        elphbolt_compatible=input%eph%elphbolt )
    else
      call ph_util_interpolate( iq2-iq1+1, vql(:, iq1:iq2), eph_ph_mfi, eph_ph_DR, dynmat, &
        minimal_distances=.true., &
        elphbolt_compatible=input%eph%elphbolt ) 
    end if

    ! diagonalize dynamical matrix
    allocate( eval(eph_nmode_tot, iq1:iq2), evec(3*natmtot, eph_nmode_tot) )
    do iq = iq1, iq2
      call ph_util_diag_dynmat( dynmat(:, :, iq-iq1+1), eval(:, iq), evec )
      ! store eigenvectors in dynmat
      dynmat(:, :, iq-iq1+1) = evec
      phfreq(:, iq) = eval(mrng(1):mrng(2), iq)
      phevec(:, :, iq) = evec(:, mrng(1):mrng(2))
    end do

    ! distribute results
#ifdef MPI
    do jrank = 0, mpi%procs-1
      if (mpi%rank == jrank) then
        do irank = 0, mpi%procs-1
          if (irank == jrank) cycle
          do iq = iq_range(1, irank+1), iq_range(2, irank+1)
            call MPI_Recv( phfreq(:, iq), nmode, MPI_DOUBLE, irank, irank*nq+iq, MPI_COMM(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
            call MPI_Recv( phevec(:, :, iq), nmode*3*natmtot, MPI_DOUBLE_COMPLEX, irank, (mpi%procs+irank)*nq+iq, MPI_COMM(mpi%comm), MPI_STATUS_IGNORE, mpi%ierr )
          end do
        end do
      else
        do iq = iq1, iq2
          evec = dynmat(mrng_list(1, jrank+1):mrng_list(2, jrank+1), :, iq-iq1+1)
          call MPI_Send( eval(mrng_list(1, jrank+1), iq), size( evec, dim=1 ), MPI_DOUBLE, jrank, mpi%rank*nq+iq, MPI_COMM(mpi%comm), mpi%ierr )
          call MPI_Send( evec, size( evec ), MPI_DOUBLE_COMPLEX, jrank, (mpi%procs+mpi%rank)*nq+iq, MPI_COMM(mpi%comm), mpi%ierr )
        end do
      end if
    end do 
#endif
    deallocate( dynmat, eval, evec, mrng_list, iq_range, ix_range )
  end subroutine eph_ph_interpolate
  !-------------------------------------------------------------------------------- 
end module eph_phonons
