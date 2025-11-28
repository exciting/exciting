
module mod_vxc
  use asserts, only: assert
  use constants, only: zzero, zone
  use gw_io, only: file_format_binary, file_format_text
  use modmpi, only: terminate_if_false
  use m_getunit, only: getunit
  use precision, only: dp, i32

  implicit none

  private

  !> APW-APW exchange-correlation  integrals
  real(dp), public, allocatable :: vxcraa(:,:,:,:,:,:)
      
  !> local-orbital-APW exchange-correlation  integrals
  real(dp), public, allocatable :: vxcrloa(:,:,:,:,:)
      
  !> local-orbital-local-orbital exchange-correlation  integrals
  real(dp), public, allocatable :: vxcrlolo(:,:,:,:)

  !> G-space interstitial exchange-correlation potential
  complex(dp), public, allocatable :: vxcig(:)
  
  character(len=*), parameter, private :: file_name_vxcnn = 'VXCNN'
  character(len=*), parameter, private :: extension_binary_format = '.OUT'
  character(len=*), parameter, private :: extension_text_format = '.DAT'
  
  !> Type to encapsulate the information about diagonal elements of VXC
  !> and the corresponding k-points
  type, public :: vxc_diagonal_elements
    !> Array with the diagonal elements of VXC
    complex(dp), allocatable :: diag_elements(:, :)
    !> Indexes of the k-points
    integer(i32), allocatable :: kpt_indexes(:)
    !> Lattice coordinates of the k-points
    real(dp), allocatable :: kpt_lattice_coord(:, :)
  contains
    procedure, private :: init_components, deallocate_components
  end type

  !> Singleton with the diagonal matrix elements of the exchange-correlation potential
  type(vxc_diagonal_elements), public, protected :: vxcnn

  public :: calcvxcnn, read_vxcnn, write_vxcnn, deallocate_vxcnn

contains
  !> Initialize the components of an object with type `vxc_diagonal_elements`
  subroutine init_components( this, first_band, last_band, kpt_indexes, kpt_lattice_coord, vxcnn_elements )
    class(vxc_diagonal_elements), intent(out) :: this
    !> Index of the first Kohn-Sham band
    integer(i32), intent(in) :: first_band
    !> Index of the last Kohn-Sham band
    integer(i32), intent(in) :: last_band
    !> Indexes of k-points
    integer(i32), intent(in) :: kpt_indexes(:)
    !> Lattice coordinates of the k-points
    real(dp), intent(in) :: kpt_lattice_coord(:, :)
    !> When present, initialize `diag_elements` with this array
    complex(dp), intent(in), optional :: vxcnn_elements(first_band:, :)

    call assert( size(kpt_indexes) == size(kpt_lattice_coord, 2), &
      'kpt_indexes and kpt_lattice_coord have incompatible sizes' )
    call assert( size(kpt_lattice_coord, 1)==3, 'kpt_lattice_coord must have size 3 along 1st dim.' )
    if( present( vxcnn_elements ) ) then 
      call assert( size( vxcnn_elements, 2 ) == size(kpt_indexes), 'vxcnn_elements and kpt_indexes must have compatible sizes')
      call assert( ubound( vxcnn_elements, 1 ) == last_band, 'vxcnn_elements must have ubound along 1st dim. equal to last_band')
      allocate( this%diag_elements, source=vxcnn_elements )
    else 
      allocate( this%diag_elements(first_band:last_band, size(kpt_indexes) ), source=zzero )
    end if
    allocate( this%kpt_indexes, source=kpt_indexes )
    allocate( this%kpt_lattice_coord, source=kpt_lattice_coord )

  end subroutine

  
  subroutine deallocate_components( this )
    class(vxc_diagonal_elements), intent(inout) :: this

    deallocate( this%diag_elements, this%kpt_indexes, this%kpt_lattice_coord )
  end subroutine


  subroutine deallocate_vxcnn 
    call vxcnn%deallocate_components()
  end subroutine


  !> Write the diagonal part of VXC into a file
  subroutine write_vxcnn( file_format, first_band_to_write, last_band_to_write )
    !> Format of the output file
    character(len=*), intent(in) :: file_format
    !> Index of the first KS band to write to the output file
    integer(i32), intent(in) :: first_band_to_write
    !> Index of the last KS band to write to the output file
    integer(i32), intent(in) :: last_band_to_write

    call assert( allocated( vxcnn%diag_elements ), 'vxcnn%diag_elements not allocated' )
    call assert( first_band_to_write >= lbound( vxcnn%diag_elements, 1 ), 'first_band_to_write out of bounds' )
    call assert( last_band_to_write <= ubound( vxcnn%diag_elements, 1 ), 'last_band_to_write out of bounds' )

    select case( trim(file_format) )
    case( file_format_text )
      call write_vxcnn_text( first_band_to_write, last_band_to_write )
    case( file_format_binary )
      call write_vxcnn_binary( first_band_to_write, last_band_to_write )
    case default
      call terminate_if_false( .false., '(write_vxcnn): Unrecognized file_format: ' // trim(file_format) )
    end select

  end subroutine


  !> Write the diagonal part of VXC into a file with text format
  subroutine write_vxcnn_text( first_band_to_write, last_band_to_write )
    !> Index of the first KS band to write to the output file
    integer(i32), intent(in) :: first_band_to_write
    !> Index of the last KS band to write to the output file
    integer(i32), intent(in) :: last_band_to_write

    integer(i32) :: fid, ik, i, j

    call getunit( fid )
    open( fid, file=file_name_vxcnn//extension_text_format, form='FORMATTED', status='UNKNOWN' )
    write( fid, '(3I10,A)' ) first_band_to_write, last_band_to_write, &
      size( vxcnn%diag_elements, 2 ), ' : index of first band, index of last band, number of k-points '
    do j = 1, size( vxcnn%kpt_indexes )
      ik = vxcnn%kpt_indexes(j)
      write( fid, '("ik=",I6,"    vkl=",3F15.8)' ) ik, vxcnn%kpt_lattice_coord(:, j)
      do i = first_band_to_write, last_band_to_write
        write( fid, '(I4,2F16.8)' ) i, vxcnn%diag_elements(i, j)
      end do
      write( fid, * )
    end do
    close( fid )

  end subroutine


  !> Write VXC into a file with binary format
  subroutine write_vxcnn_binary( first_band_to_write, last_band_to_write )
    !> Index of the first KS band to write to the output file
    integer(i32), intent(in) :: first_band_to_write
    !> Index of the last KS band to write to the output file
    integer(i32), intent(in) :: last_band_to_write

    integer(i32) :: fid, j

    call getunit(fid)

    open( fid, file=file_name_vxcnn//extension_binary_format, Form='UNFORMATTED', Status='UNKNOWN' )
    write( fid ) first_band_to_write, last_band_to_write, size( vxcnn%diag_elements, 2 )
    do j = 1, size( vxcnn%kpt_indexes )
      write( fid ) vxcnn%kpt_indexes(j), vxcnn%kpt_lattice_coord(:, j), &
                   vxcnn%diag_elements(first_band_to_write:last_band_to_write, j)
    end do
    close( fid )

  end subroutine


  !> Read VXC stored in a file
  subroutine read_vxcnn( file_format )
    character(len=*), intent(in) :: file_format

    select case( trim(file_format) )
    case ('text')
      call read_vxcnn_text
    case( 'binary' )
      call read_vxcnn_binary
    case default
      call terminate_if_false( .true., '(write_vxcnn): Unrecognized file_format: ' // trim(file_format) )
    end select
  end subroutine


  !> Read VXC from a file with text format
  subroutine read_vxcnn_text()
    
    integer(i32) :: first_band, last_band, n_kpts
    integer(i32) :: fid, ik, ib, integer_ignore
    character(len=1) :: char_to_ignore
    integer(i32), allocatable :: kpt_indexes(:)
    real(dp) :: x, y
    real(dp), allocatable :: kpt_lattice_coord(:, :)
    complex(dp), allocatable :: vxcnn_elements(:, :)


    call getunit(fid)
    open( fid, file=file_name_vxcnn//extension_text_format, form='FORMATTED', status='UNKNOWN', action="READ" )
    read( fid, * ) first_band, last_band, n_kpts
    allocate( kpt_indexes(n_kpts), kpt_lattice_coord(3, n_kpts) )
    allocate( vxcnn_elements(first_band:last_band, n_kpts) )
    do ik = 1, n_kpts
      read( fid, * ) char_to_ignore, kpt_indexes(ik), char_to_ignore, kpt_lattice_coord(:, ik)
      do ib = first_band, last_band
        read( fid, * ) integer_ignore, x, y 
        vxcnn_elements(ib, ik) = cmplx( x, y )
      end do
      read( fid, * )
    end do
    close( fid )
    call vxcnn%init_components( first_band, last_band, kpt_indexes, kpt_lattice_coord, vxcnn_elements )
      
  end subroutine
  

  !> Read VXC from a file with binary format
  subroutine read_vxcnn_binary()
    
    integer(i32) :: first_band, last_band, n_kpts
    integer(i32) :: fid, ik
    integer(i32), allocatable :: kpt_indexes(:)
    real(dp), allocatable :: kpt_lattice_coord(:, :)
    complex(dp), allocatable :: vxcnn_elements(:, :)


    call getunit(fid)
    open( fid, file=file_name_vxcnn//extension_binary_format, form='UNFORMATTED', status='UNKNOWN', action="READ" )
    read( fid ) first_band, last_band, n_kpts
    allocate( kpt_indexes(n_kpts), kpt_lattice_coord(3, n_kpts) )
    allocate( vxcnn_elements(first_band:last_band, n_kpts) )
    do ik = 1, n_kpts
      read( fid ) kpt_indexes(ik), kpt_lattice_coord(:, ik), vxcnn_elements(:, ik)
    end do
    close( fid )
    call vxcnn%init_components( first_band, last_band, kpt_indexes, kpt_lattice_coord, vxcnn_elements )
      
  end subroutine

  !> This subroutine calculates the diagonal matrix elements of
  !> the exchange correlation potential (only for valence states).
  subroutine calcvxcnn( first_band, last_band, kpt_indexes, kpt_lattice_coord, mpi_env )
    use exciting_mpi, only: mpiinfo, xmpi_allgatherv
    use genvxcig, only: generate_vxcig
    use modinput, only: input
    use mod_APW_LO, only: apwordmax, nlomax, nlotot
    use mod_LDA_LU, only: ldapu, llu
    use mod_atoms, only: natoms, natmtot, nspecies
    use mod_eigensystem, only: nmatmax
    use mod_eigenvalue_occupancy, only: nstfv
    use mod_gw_degeneracies, only: get_degenerate_limits_qp_interval_ikp, &
                                   degenerate_subspaces                               
    use mod_hybrids, only: hybridhf, vxnl
    use mod_misc, only: filext
    use mod_muffin_tin, only: lmmaxapw, lmmaxvr
    use mod_potential_and_density, only: vhalfir, vhalfmt, vxcir, vxcmt
    use modgw, only: Gkqset, Gkset, Gset, kset, kqset, time_vxc
    use modmpi, only: distribute_loop, mpiglobal
    use modxs, only: isreadstate0
    use vector_multiplication, only: dot_multiply
    use vxcrad, only: obtain_vxc_radial
    
    !> Index of the first KS band to calculate the matrix elements of vxc
    integer(i32), intent(in) :: first_band
    !> Index of the last KS band to calculate the matrix elements of vxc
    integer(i32), intent(in) :: last_band
    !> Indexes of k-points to be calculated
    integer(i32), intent(in) :: kpt_indexes(:)
    !> Lattice coordinates of the k-points
    real(dp), intent(in) :: kpt_lattice_coord(:, :)
    !> The MPI environment type, for distribution over MPI processes and 
    !> summation over MPI-distributed arrays
    type(mpiinfo), intent(in) :: mpi_env
    
    integer(i32) :: ikp, i, i_first, i_last, n_kpt, ik
    integer(i32) :: ib, dim, ia, is, ngp
    ! For the averaging over degenerate states
    integer(i32) :: ispace_init, ispace_final, ispace, lowband, upband, size_deg
    real(dp) :: tstart, tend
    real(dp), allocatable :: vxc_eff(:, :, :)
    complex(dp), allocatable :: apwalm(:,:,:,:), evecfv(:,:), h(:)
    character(80) :: filext_save
    logical :: isreadstate0_save, is_dft_half
    
    
    call timesec(tstart)
    
    if (hybridhf) then
      filext_save = filext
      isreadstate0_save = isreadstate0
      filext = '.OUT'
      isreadstate0 = .false.
      call readstate()
      filext = filext_save
      isreadstate0 = isreadstate0_save
      ! read the non-local potential used in hybrids
      call read_vxnl()
    end if
    
    ! Global array to store <n|Vxc|n>
    n_kpt = size( kpt_indexes )
    
    call vxcnn%init_components( first_band, last_band, kpt_indexes, kpt_lattice_coord )

    ! allocate exchange-correlation integral arrays
    if ( allocated( vxcraa ) ) deallocate( vxcraa )
    allocate( vxcraa(apwordmax, 0:input%groundstate%lmaxmat, apwordmax, &
      0:input%groundstate%lmaxapw, 0:lmmaxvr, natmtot) )
    if ( allocated( vxcrloa ) ) deallocate( vxcrloa )
    allocate(vxcrloa(nlomax, apwordmax, 0:input%groundstate%lmaxmat, 0:lmmaxvr, natmtot))
    if ( allocated( vxcrlolo ) ) deallocate( vxcrlolo )
    allocate( vxcrlolo(nlomax, nlomax, 0:lmmaxvr, natmtot) )
    
    is_dft_half = associated( input%groundstate%dfthalf )

    ! Here an auxiliary array is introduced to keep vxcmt unchanged in the case of DFT-1/2
    vxc_eff = vxcmt
    if( is_dft_half ) vxc_eff = vxc_eff + vhalfmt
    ! Calculate radial integrals
    call obtain_vxc_radial( vxc_eff, vxcraa, vxcrloa, vxcrlolo )

    allocate(vxcig(Gset%ngvec))
    ! Here, avoid creating a new array, as done for the MT part: typically, vxcir is very large.
    ! Instead, adopt the in-place–modify-and-restore strategy
    if( is_dft_half ) vxcir = vxcir + vhalfir
    ! Fourier transform the interstitial part of Vxc
    call generate_vxcig( vxcir, vxcig )
    ! Restore vxcir
    if( is_dft_half ) vxcir = vxcir - vhalfir
    
    allocate( apwalm(Gkset%ngkmax,apwordmax,lmmaxapw,natmtot) )
    allocate( evecfv(nmatmax,nstfv) )
    allocate( h(nmatmax) )

    call distribute_loop( mpi_env, n_kpt, i_first, i_last )
    do i = i_first, i_last
      ikp = kpt_indexes(i)
      ik = kset%ikp2ik(ikp)
      ngp = Gkqset%ngk(1,ik)
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
      call match(ngp, Gkqset%gkc(:,1,ik), Gkqset%tpgkc(:,:,1,ik), &
                  Gkqset%sfacgk(:,:,1,ik), apwalm)
    
      do ib = first_band, last_band
        h(:) = zzero
        ! muffin-tin contributions
        do is = 1, nspecies
          do ia = 1, natoms(is)
            call vxcaa(is, ia, ngp, apwalm, evecfv(:, ib), h)
            call vxcalo(is, ia, ngp, apwalm, evecfv(:, ib), h)
            call vxclolo(is, ia, ngp, evecfv(:, ib), h)
            !----------
            ! LDA+U
            !----------
            if ((ldapu /= 0) .and. (llu(is) >= 0)) then
              call vmat_ldapu(is, ia, ngp, apwalm, evecfv(:, ib), vxcnn%diag_elements(ib, i))
            end if
          end do
        end do
        ! interstitial contribution
        call vxcistl(ngp, Gkqset%igkig(:,1,ik), evecfv(:,ib), h)
        dim = ngp+nlotot
        vxcnn%diag_elements(ib, i) = vxcnn%diag_elements(ib, i) + dot_multiply(evecfv(1:dim, ib), h(1:dim), conjg_a=.true.)
      end do ! i
    
      if (hybridhf) then
        ! setup the hybrid Vxc
        do ib = first_band, last_band
          vxcnn%diag_elements(ib,i) = vxcnn%diag_elements(ib,i)+ &
            input%groundstate%Hybrid%excoeff*vxnl(ib,ib,ikp)
        end do
      end if
    end do ! ikp
    
    deallocate(vxcraa)
    deallocate(vxcrloa)
    deallocate(vxcrlolo)
    if (hybridhf) deallocate(vxnl)
    
    call xmpi_allgatherv( mpiglobal, vxcnn%diag_elements, &
      (last_band - first_band + 1) * (i_last - i_first + 1) )
    ! Here we enforce degeneracies in the VXCNN (the degeneracy lifting is a numerical artifact here)
    do i = i_first, i_last
      ikp = kpt_indexes(i)
      ! First we compute indeces of the subspaces we are interested in
      call get_degenerate_limits_qp_interval_ikp( ikp, ispace_init, ispace_final )
      ! Averaging
      do ispace = ispace_init, ispace_final
        lowband  = degenerate_subspaces(1, ispace, ikp)
        upband   = degenerate_subspaces(2, ispace, ikp)
        call assert( lowband >= first_band, 'lowband is smaller than first_band' )
        call assert( upband <= last_band, 'upband is larger than last_band' )
        size_deg = degenerate_subspaces(3, ispace, ikp)
        vxcnn%diag_elements(lowband:upband, i) = sum( vxcnn%diag_elements(lowband:upband, i) ) / size_deg
      end do
    end do
    
    call timesec(tend)
    time_vxc = time_vxc + tend - tstart
  end subroutine


end module
