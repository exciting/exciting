module matrix_fourier_interpolation
  use precision, only: dp
  use modmpi, only: terminate_if_false

  implicit none
  private

  !> effective zero for real/reciprocal space
  real(dp), parameter :: epslat = 1e-12_dp

  !> object to set up and perform a matrix Fourier interpolation (MFI)
  type mfi_type
    !> number of reciprocal space grid points in each dimension, \(N_1,N_2,N_3\)
    !> (equivalent to corresponding real space super cell size)
    integer :: ngrid(3)

    ! reciprocal space grid
    !> reciprocal space lattice vectors, \({\bf b}_1,{\bf b}_2,{\bf b}_3\) (columnwise)
    real(dp) :: bvec(3, 3)
    !> number of reciprocal space grid points, \(N_{\bf p}\)
    integer :: np
    !> list of reciprocal space grid points \({\bf p}\) in lattice coordinates
    real(dp), allocatable :: vpl(:,:)
    !> list of reciprocal space grid points \({\bf p}\) in Cartesian coordinates
    real(dp), allocatable :: vpc(:,:)

    ! real space grid
    !> real space lattice vectors, \({\bf a}_1,{\bf a}_2,{\bf a}_3\) (columnwise)
    real(dp) :: avec(3, 3)
    !> number of real space grid points, \(N_{\bf R}\)
    integer :: nr
    !> list of real space grid points \({\bf R}\) in lattice coordinates
    integer, allocatable :: vrl(:,:)
    !> list of real space grid points \({\bf R}\) in Cartesian coordinates
    real(dp), allocatable :: vrc(:,:)
    !> length of each real space grid point
    real(dp), allocatable :: rlen(:)
    !> multiplicity \(M_{\bf R}\) of each real space grid point
    integer, allocatable :: rmul(:)
    !> index of the opposite point \(-{\bf R}\) for each real space grid point
    integer, allocatable :: rneg(:)

    ! minimum distance interpolation
    !> localization centers \({\bf \tau}_m\) of orbitals on the left of matrix elements (bra) in lattice coordinates
    real(dp), allocatable :: left_centers(:,:)
    !> localization centers \({\bf \tau}_n\) of orbitals on the right of matrix elements (ket) in lattice coordinates
    real(dp), allocatable :: right_centers(:,:)
    !> minimal distance real space lattice vectors for each pair of centers in lattice coordinates,
    !> \({\bf R}_{mn,i}({\bf R}) = {\bf R} + {\bf S}_{mn,i}({\bf R}) \)
    integer, allocatable :: mindist_vrl(:,:,:,:,:)
    !> multiplicity \(M_{mn}({\bf R})\) of real space lattice vectors for each pair of centers
    integer, allocatable :: mindist_rmul(:,:,:)

    contains

      !> destructor, free memory
      procedure :: destroy
      !> set localization centers and find minimal distance vectors
      procedure :: set_localization_centers
      !> Fourier transform from \({\bf p}\)-grid (coarse grid) to R-grid
      procedure :: transform_p2R
      !> Fourier transform from \({\bf R}\)-grid to any reciprocal space point \({\bf p}\)
      procedure :: transform_R2p
      !> print set of p-vectors to text file
      procedure :: print_pset
      !> print set of R-vectors to text file
      procedure :: print_rset
  end type mfi_type
  !> constructor
  interface mfi_type
    module procedure :: generate_mfi_type
  end interface

  public :: mfi_type

contains

  !> Destructor. Free memory from type-bound arrays.
  subroutine destroy( this )
    !> matrix fourier interpolation object
    class(mfi_type), intent(inout) :: this

    if( allocated( this%vpl ) ) deallocate( this%vpl )
    if( allocated( this%vpc ) ) deallocate( this%vpc )
    if( allocated( this%vrl ) ) deallocate( this%vrl )
    if( allocated( this%vrc ) ) deallocate( this%vrc )
    if( allocated( this%rlen ) ) deallocate( this%rlen )
    if( allocated( this%rmul ) ) deallocate( this%rmul )
    if( allocated( this%rneg ) ) deallocate( this%rneg )
    if( allocated( this%left_centers ) ) deallocate( this%left_centers )
    if( allocated( this%right_centers ) ) deallocate( this%right_centers )
    if( allocated( this%mindist_vrl ) ) deallocate( this%mindist_vrl )
    if( allocated( this%mindist_rmul ) ) deallocate( this%mindist_rmul )
  end subroutine destroy

  !> Constructor. Set up a matrix fourier interpolation
  function generate_mfi_type( ngrid, bvec, vpl ) result( mfi )
    use constants, only: twopi
    !> number of reciprocal space grid points in each dimension, \(N_1,N_2,N_3\)
    !> (equivalent to corresponding real space super cell size)
    integer, intent(in) :: ngrid(3)
    !> reciprocal space lattice vectors, \({\bf b}_1,{\bf b}_2,{\bf b}_3\) (columnwise)
    real(dp), intent(in) :: bvec(3, 3)
    !> list of reciprocal space grid points \({\bf p}\) in lattice coordinates
    real(dp), intent(in) :: vpl(:,:)
    !> matrix fourier interpolation setup
    type(mfi_type) :: mfi

    integer :: ip, iv(3)
    real(dp) :: vploff(3), rv(3)

    logical, allocatable :: grid_point_covered(:,:,:)

    real(8), external :: r3mdet

    ! set grid dimensions
    mfi%ngrid = ngrid
    call terminate_if_false( all( mfi%ngrid > 0 ), '(setup_mfi): &
      All grid dimensions must be positive.' )
    ! set reciprocal space lattice vectors
    mfi%bvec = bvec
    call terminate_if_false( abs( r3mdet( bvec ) ) > epslat, '(setup_mfi): &
      Reciprocal space lattice vectors are not linearly independent.' )
    ! set number of p-vectors
    mfi%np = product( mfi%ngrid )
    ! set real space lattice vectors
    call r3minv( mfi%bvec, mfi%avec )
    mfi%avec = twopi * transpose( mfi%avec )
    ! set reciprocal space grid points
    call terminate_if_false( size( vpl, dim=1 ) == 3, '(setup_mfi): &
      Reciprocal space grid points must be vectors of length 3.' )
    call terminate_if_false( size( vpl, dim=2 ) == mfi%np, '(setup_mfi): &
      Number of given vectors does not equal number of reciprocal space grid points.' )
    vploff = vpl(:, 1) * mfi%ngrid
    vploff = (vploff - floor( vploff )) / mfi%ngrid
    allocate( mfi%vpl(3, mfi%np), mfi%vpc(3, mfi%np) )
    allocate( grid_point_covered(0:mfi%ngrid(1)-1, 0:mfi%ngrid(2)-1, 0:mfi%ngrid(3)-1), source=.false. )
    do ip = 1, mfi%np
      rv = (vpl(:, ip) - vploff) * mfi%ngrid
      iv = nint( rv )
      call terminate_if_false( norm2( rv - iv ) < epslat, '(setup_mfi): &
        Given reciprocal space grid points are not compatible with given grid size.' )
      iv = modulo( iv, mfi%ngrid )
      grid_point_covered(iv(1), iv(2), iv(3)) = .true.
      mfi%vpl(:, ip) = vpl(:, ip)
      call r3mv( mfi%bvec, mfi%vpl(:, ip), mfi%vpc(:, ip) )
    end do
    call terminate_if_false( all( grid_point_covered ), '(setup_mfi): &
      The given reciprocal space vectors do not cover the full grid.' )
    ! find Wigner-Seitz real space lattice vectors
    call generate_r_vectors( mfi%avec, mfi%ngrid, epslat, &
      mfi%nr, mfi%vrl, mfi%vrc, mfi%rlen, mfi%rmul, mfi%rneg )
  end function generate_mfi_type

  !> Set localization centers of orbitals, \({\bf \tau}_m\) and \({\bf \tau}_n\), on the left
  !> and on the right, respectively, and for each real space lattice vector \({\bf R}\) find all
  !> super cell lattice vectors \({\bf S}_{i,mn}({\bf R})\) such that 
  !> \({\bf R}-{\bf \tau}_m+{\bf \tau}_n+{\bf S}_{i,mn}({\bf R})\) is inside the Wigner-Seitz cell
  !> of the super cell, where \(i=1,\dots,M_{mn}({\bf R})\) labels all possible such vectors.
  subroutine set_localization_centers( this, left_centers, right_centers, coordinates )
    use constants, only: twopi
    !> matrix fourier interpolation object
    class(mfi_type), intent(inout) :: this
    !> localization centers of orbitals on the left of matrix elements (bra)
    real(dp), optional, intent(in) :: left_centers(:,:)
    !> localization centers of orbitals on the right of matrix elements (ket)
    real(dp), optional, intent(in) :: right_centers(:,:)
    !> coordinates in which centers are given 
    !> (`'l'` or `'lattice'` for lattice (default), `'c'` or `'cartesian'` for Cartesian)
    character(*), optional, intent(in) :: coordinates

    integer :: ncl, ncr, icl, icr, ir, im
    character :: coord
    real(dp) :: vl(3)

    integer, allocatable :: wrap(:,:), vrl(:,:,:,:,:)

    coord = 'l'
    if( present( coordinates ) ) then
      call terminate_if_false( &
        any( trim( coordinates ) == ['l', 'L', 'c', 'C'] ) .or. &
        any( trim( coordinates ) == ['lattice', 'Lattice', 'LATTICE'] ) .or. &
        any( trim( coordinates ) == ['cartesian', 'Cartesian', 'CARTESIAN'] ), '(set_localization_centers): &
        Invalid value for `coordinates`.' )
      if( any( trim( coordinates ) == ['c', 'C'] ) .or. &
          any( trim( coordinates ) == ['cartesian', 'Cartesian', 'CARTESIAN'] ) ) coord = 'c'
    end if

    ! set left centers
    if( present( left_centers ) ) then
      call terminate_if_false( size( left_centers, dim=1 ) == 3, '(set_localization_centers): &
        Left localization centers must be vectors of length 3.' )
      if( coord == 'l' ) then
        allocate( this%left_centers, source=left_centers )
      else
        allocate( this%left_centers, source=matmul( transpose( this%bvec ), left_centers )/twopi )
      end if
    else
      allocate( this%left_centers(3, 1), source=0.0_dp )
    end if

    ! set right centers
    if( present( right_centers ) ) then
      call terminate_if_false( size( right_centers, dim=1 ) == 3, '(set_localization_centers): &
        Right localization centers must be vectors of length.' )
      if( coord == 'l' ) then
        allocate( this%right_centers, source=right_centers )
      else
        allocate( this%right_centers, source=matmul( transpose( this%bvec ), right_centers )/twopi )
      end if
    else
      allocate( this%right_centers(3, 1), source=0.0_dp )
    end if

    ! find multiplicities and minimal distance vectors
    ncl = size( this%left_centers, dim=2 )
    ncr = size( this%right_centers, dim=2 )
    allocate( this%mindist_rmul(ncl, ncr, this%nr), source=0 )
    allocate( vrl(3, 16, ncl, ncr, this%nr), source=0 )
!$omp parallel default( shared ) private( ir, icr, icl, vl, wrap, im )
!$omp do collapse(3)
    do ir = 1, this%nr
      do icr = 1, ncr
        do icl = 1, ncl
          vl = this%vrl(:, ir) - this%left_centers(:, icl) + this%right_centers(:, icr)
          vl = vl / this%ngrid
          wrap = ws_wrapping_vectors( this%avec, vl, epslat )
          this%mindist_rmul(icl, icr, ir) = size( wrap, dim=2 )
          do im = 1, this%mindist_rmul(icl, icr, ir)
            vrl(:, im, icl, icr, ir) = &
              nint( (vl + wrap(:, im)) * this%ngrid + this%left_centers(:, icl) - this%right_centers(:, icr) )
          end do
        end do
      end do
    end do
!$omp end do
!$omp end parallel
    im = maxval( this%mindist_rmul )
    allocate( this%mindist_vrl(3, im, ncl, ncr, this%nr), source=vrl(:, 1:im, :, :, :) )
  end subroutine set_localization_centers

  !> Perform the Fourier transform from the coarse reciprocal space grid \({\bf p}\)
  !> to the real space grid \({\bf R}\), i.e., for a set of matrices \(M_{mn,i}({\bf p})\)
  !> with \(i=1,\dots,N_M\), compute
  !> \[ M_{mn,i}({\bf R}) = \frac{1}{N_{\bf p}} \sum_{\bf p} {\rm e}^{-{\rm i}{\bf p}\cdot{\bf R}}
  !>    M_{mn,i}({\bf p}) \;. \]
  subroutine transform_p2R( this, matsize, nummat, Mp, ldp1, ldp2, MR, ldr1, ldr2 )
    use constants, only: zzero, zone
    !> matrix fourier interpolation object
    class(mfi_type), intent(in) :: this
    !> size of each matrix \(\mathrm{\bf M}\)
    integer, intent(in) :: matsize(2)
    !> number of matrices \(\mathrm{\bf M}\), \(N_M\)
    integer, intent(in) :: nummat
    !> matrices \(\mathrm{\bf M}({\bf p})\) on coarse recirpocal space grid
    complex(dp), intent(in) :: Mp(ldp1, ldp2, *)
    !> leading dimensions of array `Mp`
    integer, intent(in) :: ldp1, ldp2
    !> matrices \(\mathrm{\bf M}({\bf R})\) on real space grid
    complex(dp), intent(out) :: MR(ldr1, ldr2, *)
    !> leading dimensions of array `MR`
    integer, intent(in) :: ldr1, ldr2

    real(dp), allocatable :: phases(:,:)
    complex(dp), allocatable :: weights(:,:)

    ! sanity check
    call terminate_if_false( all( matsize > 0 ), '(transform_p2R): &
      Matrix dimensions must be positive.' )
    call terminate_if_false( nummat > 0, '(transform_p2R): &
      Number of matrices must be positive.' )
    call terminate_if_false( product( matsize ) == ldp1 .and. product( matsize ) == ldr1, '(transform_p2R): &
      Number of elements per matrix must be equal to first leading dimension.' )
    call terminate_if_false( nummat <= ldp2 .and. nummat <= ldr2, '(transform_p2R): &
      Number of matrices must not be larger than second leading dimension.' )

    ! generate Fourier weights
    allocate( phases(this%np, this%nr) )
    allocate( weights(this%np, this%nr) )
    call dgemm( 't', 'n', this%np, this%nr, 3, 1.0_dp, &
           this%vpc, 3, &
           this%vrc, 3, 0.0_dp, &
           phases, this%np )
    weights = cmplx( cos( phases ), -sin( phases ), dp )

    ! perform Fourier transform
    call zgemm( 'n', 'n', product(matsize)*nummat, this%nr, this%np, zone/this%np, &
           Mp, ldp1*ldp2, &
           weights, this%np, zzero, &
           MR, ldr1*ldr2 )
  end subroutine transform_p2R

  !> Perform the Fourier transform from the real space grid \({\bf R}\)
  !> to an arbitrary set of reciprocal space points \({\bf p}\), i.e., for a set of matrices \(M_{mn,i}({\bf R})\)
  !> with \(i=1,\dots,N_M\), compute
  !> \[ M_{mn,i}({\bf p}) = \sum_{\bf R} \frac{1}{M_{\bf R}} {\rm e}^{{\rm i}{\bf p}\cdot{\bf R}}
  !>    M_{mn,i}({\bf R}) \;. \]
  !> If minimal distance interpolation is employed, instead, 
  !> \[ M_{mn,i}({\bf p}) = \sum_{\bf R} \frac{1}{M_{\bf R}} w_{mn}({\bf R},{\bf p})\, M_{mn,i}({\bf R}) \]
  !> is computed with Fourier weights
  !> \[ w_{mn}({\bf R},{\bf p}) = \frac{1}{M_{mn}({\bf R})} \sum_{i=1}^{M_{mn}({\bf R})}
  !>    {\rm e}^{{\rm i}{\bf p}\cdot({\bf R}+{\bf S}_{mn,i}({\bf R}))} \;. \]
  !> Note, that the orbital localization centers \({\bf \tau}_m\) and \({\bf \tau}_n\) must be set
  !> in advance using [[set_localization_centers(subroutine)]].
  subroutine transform_R2p( this, matsize, nummat, MR, ldr1, ldr2, Mp, ldp1, ldp2, vpl, minimal_distances )
    use constants, only: zzero, zone, twopi
    !> matrix fourier interpolation object
    class(mfi_type), intent(in) :: this
    !> size of each matrix \(\mathrm{\bf M}\)
    integer, intent(in) :: matsize(2)
    !> number of matrices \(\mathrm{\bf M}\), \(N_M\)
    integer, intent(in) :: nummat
    !> matrices \(\mathrm{\bf M}({\bf R})\) on real space grid
    complex(dp), intent(in) :: MR(ldr1, ldr2, *)
    !> leading dimensions of array `MR`
    integer, intent(in) :: ldr1, ldr2
    !> matrices \(\mathrm{\bf M}({\bf p})\) on recirpocal space interpolation points
    complex(dp), intent(out) :: Mp(ldp1, ldp2, *)
    !> leading dimensions of array `Mp`
    integer, intent(in) :: ldp1, ldp2
    !> set of reciprocal space interpolation points \({\bf p}\) in lattice coordinates
    real(dp), intent(in) :: vpl(:,:)
    !> `.true.` if minimal distance interpolation should be used (default: `.false.`)
    logical, optional, intent(in) :: minimal_distances

    integer :: np, nelem, ip, ir, i, m, n, mn, icl, icr, &
               left_skip, right_skip, ncl, ncr
    logical :: mindist

    integer, allocatable :: wgt_idx_map(:,:)
    real(dp), allocatable :: phases(:,:)
    complex(dp), allocatable :: weights(:,:), weights1d(:)

    ! sanity check
    call terminate_if_false( all( matsize > 0 ), '(transform_R2p): &
      Matrix dimensions must be positive.' )
    call terminate_if_false( nummat > 0, '(transform_R2p): &
      Number of matrices must be positive.' )
    call terminate_if_false( product( matsize ) == ldp1 .and. product( matsize ) == ldr1, '(transform_R2p): &
      Number of elements per matrix must be equal to first leading dimension.' )
    call terminate_if_false( nummat <= ldp2 .and. nummat <= ldr2, '(transform_R2p): &
      Number of matrices must not be larger than second leading dimension.' )
    call terminate_if_false( size(vpl, dim=1) == 3, '(transform_R2p): &
      Interpolation points must be vectors of length 3.' )

    np = size( vpl, dim=2 )
    nelem = product( matsize )

    mindist = .false.
    if( present( minimal_distances ) ) mindist = minimal_distances

    ! simple transform
    if( .not. mindist ) then
      ! generate Fourier weights
      allocate( phases(np, this%nr) )
      allocate( weights(np, this%nr) )
      call dgemm( 't', 'n', np, this%nr, 3, twopi, &
             vpl, 3, &
             dble(this%vrl), 3, 0.0_dp, &
             phases, np )
      weights = cmplx( cos( phases ), sin( phases ), dp )
!$omp parallel default(shared) private(ir)
!$omp do
      do ir = 1, this%nr
        weights(:, ir) = weights(:, ir) / this%rmul(ir)
      end do
!$omp end do
!$omp end parallel

      ! perform Fourier transform
      call zgemm( 'n', 't', nelem*nummat, np, this%nr, zone, &
             MR, ldr1*ldr2, &
             weights, np, zzero, &
             Mp, ldp1*ldp2 )
    ! minimal distance transform
    else
      ! check localization centers
      call terminate_if_false( allocated( this%left_centers ) .and. allocated( this%right_centers ), '(transform_R2p): &
        Localization centers must be set for minimal distance interpolation. &
        Call `set_localization_centers` first.' )
      ncl = size( this%left_centers, dim=2 )
      ncr = size( this%right_centers, dim=2 )
      left_skip = matsize(1) / ncl
      right_skip = matsize(2) / ncr
      call terminate_if_false( left_skip*ncl == matsize(1), '(transform_R2p): &
        First matrix dimension must be an integer multiple of number of left localization centers.' )
      call terminate_if_false( right_skip*ncr == matsize(2), '(transform_R2p): &
        Second matrix dimension must be an integer multiple of number of right localization centers.' )

      ! build 1d -> 2d index map for weights
      allocate( weights1d(nelem), wgt_idx_map(nelem, 2) )
      do n = 1, matsize(2)
        icr = (n - 1) / right_skip + 1
        do m = 1, matsize(1)
          icl = (m - 1) / left_skip + 1
          mn = (n - 1) * matsize(1) + m
          wgt_idx_map(mn, :) = [icl, icr]
        end do
      end do

      ! perform Fourier transform
      Mp(:, :, 1:np) = zzero
!$omp parallel default( shared ) private( ir, weights, weights1d, i )
!$omp do
      do ip = 1, np
        do ir = 1, this%nr
          ! get Fourier weights
          weights = generate_mindist_weights( vpl(:, ip), this%mindist_vrl(:, :, :, :, ir), this%mindist_rmul(:, :, ir) ) / this%rmul(ir)
          weights1d = [(weights(wgt_idx_map(i, 1), wgt_idx_map(i, 2)), i=1, nelem)]
          ! add up result
          do i = 1, nummat
            Mp(1:nelem, i, ip) = Mp(1:nelem, i, ip) + weights1d * MR(1:nelem, i, ir)
          end do
        end do
      end do
!$omp end do
!$omp end parallel

      deallocate( weights1d, wgt_idx_map )
    end if
  end subroutine transform_R2p

  !> Generate the set of real space lattice vectors \({\bf R}\) that are within the
  !> Wigner-Seitz cell of the \(N_1 \times N_2 \times N_3\) super cell. Their multiplicity
  !> counts the number of equivalent vectors, i.e., lattice vectors that differ by
  !> a super cell lattice vector.
  subroutine generate_r_vectors( avec, ngrid, eps, nr, vrl, vrc, rlen, rmul, rneg )
    use sorting, only: sort_index_1d
    !> real space lattice vectors (columnwise)
    real(dp), intent(in) :: avec(3, 3)
    !> real space super cell size \(N_1,N_2,N_3\)
    integer, intent(in) :: ngrid(3)
    !> tolerance
    real(dp), intent(in) :: eps
    !> number of real space lattice vectors
    integer, intent(out) :: nr
    !> real space lattice vectors in lattice coordinates
    integer, allocatable, intent(out) :: vrl(:,:)
    !> real space lattice vectors in Cartesian coordinates
    real(dp), allocatable, intent(out) :: vrc(:,:)
    !> length of each real space lattice vector
    real(dp), allocatable, intent(out) :: rlen(:)
    !> multiplicity of each real space lattice vector
    integer, allocatable, intent(out) :: rmul(:)
    !> index of the negative vector for each real space lattice vector
    integer, allocatable, intent(out) :: rneg(:)

    integer :: i, j, k, l, mul
    real(dp) :: vl(3)

    integer, allocatable :: wrappers(:,:), idx(:)

    ! find number of R-vectors
    nr = 0
    do k = 0, ngrid(3) - 1
      do j = 0, ngrid(2) - 1
        do i = 0, ngrid(1) - 1
          vl = dble( [i,j,k] ) / ngrid
          wrappers = ws_wrapping_vectors( avec, vl, eps )
          nr = nr + size( wrappers, dim=2 )
        end do
      end do
    end do

    ! allocate arrays
    allocate( vrl(3, nr), vrc(3, nr), rlen(nr), rmul(nr), rneg(nr) )

    ! find R-vectors
    nr = 0
    do k = 0, ngrid(3) - 1
      do j = 0, ngrid(2) - 1
        do i = 0, ngrid(1) - 1
          vl = dble( [i,j,k] ) / ngrid
          wrappers = ws_wrapping_vectors( avec, vl, eps )
          mul = size( wrappers, dim=2 )
          do l = 1, mul
            nr = nr + 1
            vrl(:, nr) = nint( (vl + wrappers(:, l)) * ngrid )
            call r3mv( avec, dble( vrl(:, nr) ), vrc(:, nr) )
            rlen(nr) = norm2( vrc(:, nr) )
            rmul(nr) = mul
          end do
        end do
      end do
    end do

    ! consistency check
    call terminate_if_false( nr == size( rmul ), '(generate_r_vectors): &
      Inconsistent number of real space lattice vectors.' )

    ! sort R-vectors with increasing length
    allocate( idx(nr) )
    idx = sort_index_1d( nr, rlen, 1 )
    vrl = vrl(:, idx)
    vrc = vrc(:, idx)
    rlen = rlen(idx)
    rmul = rmul(idx)

    ! find indices of opposite vectors
    rneg = 0
    do i = 1, nr
      if( rneg(i) > 0 ) cycle
      inner: do j = 1, nr
        if( all( vrl(:, i) + vrl(:, j) == 0 ) ) then
          rneg(i) = j
          rneg(j) = i
          exit inner
        end if
      end do inner
    end do
    call terminate_if_false( all( rneg > 0 ), '(generate_r_vectors): &
      Not all opposite vectors contained in set.' )
  end subroutine generate_r_vectors

  !> For a given reciprocal space point \({\bf p}\) and a set of minimal distance real space
  !> lattice vectors \({\bf R}_{mn,i}({\bf R}) = {\bf R}+{\bf S}_{mn,i}({\bf R})\) with \(i=1,\dots,M_{mn}({\bf R})\)
  !> compute the Fourier weights 
  !> \[ w_{mn}({\bf R},{\bf p}) = \frac{1}{M_{mn}({\bf R})} \sum_{i=1}^{M_{mn}({\bf R})} 
  !>    {\rm e}^{{\rm i}{\bf p} \cdot {\bf R}_{mn,i}({\bf R})} \;. \]
  function generate_mindist_weights( vpl, vrl, rmul ) result( weights )
    use constants, only: zzero, twopi
    !> reciprocal space point \({\bf p}\) in lattice coordinates
    real(dp), intent(in) :: vpl(3)
    !> set of minimal distance real space lattice vectors \({\bf R}_{mn,i}({\bf R})\) in lattice coordinates
    integer, intent(in) :: vrl(:,:,:,:)
    !> multiplicity \(M_{mn}({\bf R})\) of minimal distance real space lattice vectors
    integer, intent(in) :: rmul(:,:)
    !> minimal distance Fourier weights \(w_{mn}({\bf R},{\bf p})\)
    complex(dp), allocatable :: weights(:,:)

    integer :: ncl, ncr, icl, icr, ir
    real(dp) :: phase, vpl2(3)

    ncl = size( vrl, dim=3 )
    ncr = size( vrl, dim=4 )
    vpl2 = twopi * vpl

    if( allocated( weights ) ) deallocate( weights )
    allocate( weights(ncl, ncr), source=zzero )

    do icr = 1, ncr
      do icl = 1, ncl
        do ir = 1, rmul(icl, icr)
          phase = vpl2(1) * vrl(1, ir, icl, icr) + vpl2(2) * vrl(2, ir, icl, icr) + vpl2(3) * vrl(3, ir, icl, icr)
          weights(icl, icr) = weights(icl, icr) + cmplx( cos( phase ), sin( phase ), dp )
        end do
        weights(icl, icr) = weights(icl, icr) / rmul(icl, icr)
      end do
    end do
  end function generate_mindist_weights

  !> For a given vector \({\bf s}\), this function returns all the lattice vectors \({\bf S}\) 
  !> such that \({\bf s} + {\bf S}\) is in the Wigner-Seitz cell of the lattice defined by the
  !> given lattice vectors. There might be multiple such wrapping vectors when \({\bf s}\) gets
  !> wrapped onto the boundary of the Wigner-Seitz cell.
  function ws_wrapping_vectors( lvec, vl, eps ) result( wrappers )
    !> lattice vectors (columnwise)
    real(dp), intent(in) :: lvec(3, 3)
    !> input vector \({\bf s}\) in lattice coordinates
    real(dp), intent(in) :: vl(3)
    !> tolerance
    real(dp), intent(in) :: eps
    !> list of wrapping vectors \({\bf S}\)
    integer, allocatable :: wrappers(:,:)

    ! super cell search size
    integer, parameter :: sc_size = 3

    integer :: i, j, k, n, iv0(3)
    real(dp) :: d, wl(3), wc(3)

    wl = vl
    call r3frac( epslat, wl, iv0 )

    ! find length and number of shortest vectors
    n = 0
    d = huge( 0.0_dp )
    do k = -sc_size, sc_size
      do j = -sc_size, sc_size
        do i = -sc_size, sc_size
          call r3mv( lvec, wl+[i,j,k], wc )
          if( abs( norm2( wc ) - d ) < eps ) then
            n = n + 1
          else if( norm2( wc ) < d ) then
            n = 1
            d = norm2( wc )
          end if
        end do
      end do
    end do

    ! find all wrapping vectors
    allocate( wrappers(3, n) )
    n = 0
    do k = -sc_size, sc_size
      do j = -sc_size, sc_size
        do i = -sc_size, sc_size
          call r3mv( lvec, wl+[i,j,k], wc )
          if( abs( norm2( wc ) - d ) < eps ) then
            n = n + 1
            wrappers(:, n) = [i,j,k] - iv0
          end if
        end do
      end do
    end do
  end function ws_wrapping_vectors

  !> Print set of \({\bf p}\)-vectors to text file.
  subroutine print_pset( this, un, decimals )
    !> matrix fourier interpolation object
    class(mfi_type), intent(in) :: this
    !> output unit
    integer, intent(in) :: un
    !> number of decimals to print (default: `16`)
    integer, optional, intent(in) :: decimals

    integer :: d, ip
    character(6) :: fmt

    d = 16
    if( present( decimals ) ) d = decimals
    write( fmt, '("g",i2.2,".",i2.2)' ) d+8, d

    write( un, '("# reciprocal lattice vectors (in rows)")' )
    write( un, '(3'//fmt//')' ) this%bvec
    write( un, '("# reciprocal space vectors (lattice / Cartesian coordinates)")' )
    do ip = 1, this%np
      write( un, '(i6,6x,3'//fmt//',6x,3'//fmt//')' ) ip, this%vpl(:, ip), this%vpc(:, ip)
    end do
  end subroutine print_pset

  !> Print set of \({\bf R}\)-vectors to text file.
  subroutine print_rset( this, un, decimals )
    !> matrix fourier interpolation object
    class(mfi_type), intent(in) :: this
    !> output unit
    integer, intent(in) :: un
    !> number of decimals to print (default: `16`)
    integer, optional, intent(in) :: decimals

    integer :: d, ir
    character(6) :: fmt

    d = 16
    if( present( decimals ) ) d = decimals
    write( fmt, '("g",i2.2,".",i2.2)' ) d+8, d

    write( un, '("# direct lattice vectors (in rows)")' )
    write( un, '(3'//fmt//')' ) this%avec
    write( un, '("# real space lattice vectors (lattice / Cartesian coordinates) and multiplicity")' )
    do ir = 1, this%nr
      write( un, '(i6,6x,3i6,6x,3'//fmt//',6x,i6)' ) ir, this%vrl(:, ir), this%vrc(:, ir), this%rmul(ir)
    end do
  end subroutine print_rset

end module matrix_fourier_interpolation
