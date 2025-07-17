
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!

!> symmetry variables  
module mod_symmetry

      use precision, only: i32, dp

      implicit none
      
      private
      public  :: symmetrize_real_mt, symmetrize_real_ir, symapp_zfig, &
                 find_equivalent_wavevectors, what_maps_q_2_qp, get_equivalent_qpairs

      !> nosym is .true. if no symmetry information should be used
      !> replaced by inputstructurelogical::nosym
      !> number of Bravais lattice point group symmetries
      integer(i32), public :: nsymlat
      !> Bravais lattice point group symmetries
      integer(i32), public :: symlat (3, 3, 48)
      !> determinants of lattice symmetry matrices (1 or -1)
      integer(i32), public :: symlatd (48)
      !> index to inverses of the lattice symmetries
      integer(i32), public :: isymlat (48)
      !> lattice point group symmetries in Cartesian coordinates
      real(dp), public :: symlatc (3, 3, 48)
      !> tshift is .true. if atomic basis is allowed to be shifted
      !> replaced by inputstructurelogical::tshift
      !> maximum of symmetries allowed
      integer(i32), public, parameter :: maxsymcrys = 192
      !> number of crystal symmetries
      integer(i32), public :: nsymcrys
      !> crystal symmetry translation vector in lattice coordinates
      real(dp), public :: vtlsymc (3, maxsymcrys)
      !> spatial rotation element in lattice point group for each crystal symmetry
      integer(i32), public :: lsplsymc (maxsymcrys)
      !> global spin rotation element in lattice point group for each crystal symmetry
      integer(i32), public :: lspnsymc (maxsymcrys)
      !> equivalent atom index for each crystal symmetry
      integer(i32), allocatable, public :: ieqatom (:, :, :)
      !> eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
      logical, allocatable, public :: eqatoms (:, :, :)
      !> number of site symmetries
      integer(i32), allocatable, public :: nsymsite (:)
      !> site symmetry spatial rotation element in lattice point group
      integer(i32), allocatable, public :: lsplsyms (:, :)
      !> site symmetry global spin rotation element in lattice point group
      integer(i32), allocatable, public :: lspnsyms (:, :)

  contains

    !> Symmetrize a real function in the muffin-tin spheres
    !> given as a real spherical harmonics expansion.
    subroutine symmetrize_real_mt( &
        f, lmax, nr, isym, &
        rstep, rstart)
      
      use constants, only: zzero
      use mod_atoms, only: nspecies, natoms, natmtot, natmmax, idxas

      !> on input: real function; 
      !> on output: symmetrized real function
      real(dp), intent(inout) :: f(:,:,:)
      !> maximum l to use in the expansion
      integer(i32), intent(in) :: lmax
      !> number of radial points per species
      integer(i32), intent(in) :: nr(:)
      !> global indices of symmetries
      integer(i32), intent(in) :: isym(:)
      !> radial step (default: 1)
      integer(i32), optional, intent(in) :: rstep
      !> first radial point per species (default: 1)
      integer(i32), optional, intent(in) :: rstart(:)

      integer(i32) :: is, ia, ias, ja, jas, &
                      ir, ir_start, ir_step, irc, nrc, &
                      nsym, i, lspl, lmmax
      real(dp) :: sc(3,3)
      
      real(dp),    allocatable :: ft(:)
      complex(dp), allocatable :: zf(:,:,:), szf(:,:)

      ir_step = 1
      if( present( rstep)) ir_step = rstep

      lmmax = (lmax + 1)**2
      nrc = maxval(nr) / ir_step
      nsym = size(isym)

      allocate( ft(lmmax))
      allocate( zf(lmmax,nrc,natmmax))
      allocate( szf(lmmax,nrc))

      do is = 1, nspecies
        ir_start = 1
        if( present( rstart)) ir_start = rstart(is)
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          ! make a complex copy of the input and delete input
          irc = 0
          do ir = ir_start, nr(is), ir_step
            irc = irc + 1
            call rtozflm( lmax, f(:,ir,ias), zf(:,irc,ia))
            f(:,ir,ias) = 0._dp
          end do
        end do
        nrc = irc
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          ! loop over symmetries
          do i = 1, nsym
            ! rotate function
            lspl = lsplsymc(isym(i))
            sc = symlatc(:,:,lspl)
            ja = ieqatom(ia,is,isym(i))
            jas = idxas(ja,is)
            call rotzflm( sc, lmax, nrc, lmmax, zf(:,:,ja), szf)
            ! add to output
            irc = 0
            do ir = ir_start, nr(is), ir_step
              irc = irc + 1
              call ztorflm( lmax, szf(:,irc), ft)
              f(:lmmax,ir,ias) = f(:lmmax,ir,ias) + ft/nsym
            end do
          end do
        end do
      end do

      deallocate( ft, zf, szf)
    end subroutine symmetrize_real_mt

    !> symmetrize a real function in the interstitial region
    !> given as a Fourier series on a real-space grid
    subroutine symmetrize_real_ir( &
        f, ng, ivg, intgv, ivgig, igfft, isym)

      use constants, only: zzero
      use m_zfftifc, only: zfftifc

      !> on input: real function; 
      !> on output: symmetrized real function
      real(dp), intent(inout) :: f(:)
      !> number of G-vectors to consider in the expansion
      integer(i32), intent(in) :: ng
      !> integer(i32) components of G-vectors the function is defined on
      integer(i32), intent(in) :: ivg(:,:)
      !> range of integer(i32) components of G-vectors
      integer(i32), intent(in) :: intgv(3,2)
      !> map from integer(i32) components of G-vector to its index in the list
      integer(i32), intent(in) :: ivgig(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),intgv(3,1):intgv(3,2))
      !> map from G-vector list to FFT grid
      integer(i32), intent(in) :: igfft(:)
      !> global indices of symmetries
      integer(i32), intent(in) :: isym(:)

      integer(i32) :: ngrid(3), ngrtot, nsym, i, lspl

      complex(dp), allocatable :: zf(:), szf(:)

      ngrid = intgv(:,2) - intgv(:,1) + 1
      ngrtot = product( ngrid)
      nsym = size(isym)

      allocate( zf(ngrtot), szf(ngrtot))

      ! transform function to reciprocal space
      zf(1:ngrtot) = cmplx( f(1:ngrtot), 0._dp, kind=dp)
      call zfftifc( 3, ngrid, -1, zf)

      ! loop over symmetries
      szf = zzero
      do i = 1, nsym
        lspl = lsplsymc(isym(i))
        ! rotate the function
        call symapp_zfig( symlat(:,:,lspl), vtlsymc(:,isym(i)), [0._dp,0._dp,0._dp], &
               zf, ng, ivg, igfft, .true., &
               szf, intgv, ivgig, igfft, .true.)
      end do

      ! tranform function to real space and normalize
      call zfftifc( 3, ngrid, 1, szf)
      f = dble(szf)*(1._dp/dble(nsym))

      deallocate( zf, szf)
    end subroutine symmetrize_real_ir

    !> Applies a crystal symmetry operation to an interstitial region function
    !> given by its Fourier components
    !> and adds the result to another function.
    !> The symmetry operation is defined by the translation followed by the rotation.
    !> Both functions can be given on different sets of G-vectors.
    subroutine symapp_zfig( rotl, vtl, vpl, zfig1, ng, ivg1, igfft1, fft1, zfig2, intgv2, ivgig2, igfft2, fft2)

      use constants, only: twopi
    
      !> symmetry rotation matrix in lattice coordinates
      integer(i32), intent(in) :: rotl(3,3)
      !> symmetry translation vector in lattice coordinates
      real(dp), intent(in) :: vtl(3)
      !> Bloch wavevector of function in lattice coordinates
      real(dp), intent(in) :: vpl(3)
      !> function 1 to which the symmetry operation is applied
      complex(dp), intent(in) :: zfig1(:)
      !> number of G-vectors in the expansion of the function 1
      integer(i32), intent(in) :: ng
      !> integer(i32) components of G-vectors function 1 is defined on
      integer(i32), intent(in) :: ivg1(:,:)
      !> map from G-vector list to FFT grid for function 1
      !> (not referenced if `fft1=.false.`)
      integer(i32), intent(in) :: igfft1(:)
      !> if `.true.` function 1 is given on the FFT grid; 
      !> if `.false.` function 1 is given on the G-vector grid
      logical, intent(in) :: fft1
      !> function 2 to which the result is added
      complex(dp), intent(inout) :: zfig2(:)
      !> range of integer(i32) components of G-vectors for function 2
      integer(i32), intent(in) :: intgv2(3,2)
      !> map from integer(i32) components of G-vector to 
      !> its index in the list for function 2
      integer(i32), intent(in) :: ivgig2(intgv2(1,1):intgv2(1,2),intgv2(2,1):intgv2(2,2),intgv2(3,1):intgv2(3,2))
      !> map from G-vector list to FFT grid for function 2
      !> (not referenced if `fft2=.false.`)
      integer(i32), intent(in) :: igfft2(:)
      !> if `.true.` function 2 is given on the FFT grid; 
      !> if `.false.` function 2 is given on the G-vector grid
      logical, intent(in) :: fft2
    
      integer(i32) :: irotl(3,3), ig, igf, jg, jgf, ivg(3), shift(3), ngrid(3)
      real(dp) :: v(3), phase
    
      ! get inverse of rotation matrix
      call i3minv( rotl, irotl)
      ! apply inverse rotation to p from the left
      call r3mtv( dble( irotl), vpl, v)
      ! get reciprocal lattice vector that maps R^-T.p back to the 1st BZ
      call r3frac( 1e-6_dp, v, shift)
    
      ngrid = intgv2(:,2) - intgv2(:,1) + 1
    
      do ig = 1, ng
        igf = ig
        if( fft1) igf = igfft1(ig)
        ! apply inverse rotation to G+p from the left and save integer(i32) part
        ivg(1) = irotl(1,1)*ivg1(1,ig) + irotl(2,1)*ivg1(2,ig) + irotl(3,1)*ivg1(3,ig) + shift(1)
        ivg(2) = irotl(1,2)*ivg1(1,ig) + irotl(2,2)*ivg1(2,ig) + irotl(3,2)*ivg1(3,ig) + shift(2)
        ivg(3) = irotl(1,3)*ivg1(1,ig) + irotl(2,3)*ivg1(2,ig) + irotl(3,3)*ivg1(3,ig) + shift(3)
        ! get phase factor from translation
        phase = -twopi*dot_product( dble( ivg1(:,ig)) + vpl, vtl)
        ! get index of rotated G+p vector in output G set
        ivg = modulo( ivg-intgv2(:,1), ngrid) + intgv2(:,1)
        jg = ivgig2( ivg(1), ivg(2), ivg(3))
        jgf = jg
        if( fft2) jgf = igfft2(jg)
        ! update entry in output function
        zfig2(jgf) = zfig2(jgf) + cmplx( cos(phase), sin(phase), dp)*zfig1(igf)
      end do    
    end subroutine symapp_zfig

    !> Find all pairs of points and rotations from a list of points and rotations that 
    !> are symmetry equivalent to a given target point.
    !>
    !> Points \({\bf p}\) and \({\bf p}'\) are symmetry equivalent, if there is
    !> a rotation \({\bf R}\) such that
    !> \[ {\bf R}({\bf p} + {\bf G}) = {\bf p}' \;, \]
    !> where \({\bf G}\) is an integer(i32) vector.
    pure subroutine find_equivalent_wavevectors( &
        dim, point, pointlist, npt, rotations, nrot, &
        equiv_point_idx, rotation_idx, &
        integer_vectors, tolerance, first_only, unique_only )
      use precision, only: dp
      !> spatial dimension
      integer(i32), intent(in) :: dim
      !> wavevector \({\bf p}'\) whos equivalents to find (in lattice coordiantes)
      real(dp), intent(in) :: point(dim)
      !> number of wavevectors in list
      integer(i32), intent(in) :: npt
      !> list of wavevectors \({\bf p}\) to search in (in lattice coordiantes)
      real(dp), intent(in) :: pointlist(dim, *)
      !> number of allowed rotations
      integer(i32), intent(in) :: nrot
      !> allowed rotation matrices \({\bf R}\) (in lattice coordinates)
      integer(i32), intent(in) :: rotations(dim, dim, *)
      !> list of indices of equivalent points in list
      integer(i32), allocatable, intent(out) :: equiv_point_idx(:)
      !> list of rotations that rotate quivalent point into target point
      integer(i32), allocatable, intent(out) :: rotation_idx(:)
      !> list of integer(i32) integer(i32) vectors \({\bf G}\)
      integer(i32), allocatable, optional, intent(out) :: integer_vectors(:, :)
      !> tolerance for two points beeing equivalent (default: `1e-12`)
      real(dp), optional, intent(in) :: tolerance
      !> return only first equivalent point and rotation (default: `.false.`)
      logical, optional, intent(in) :: first_only
      !> Save unique points only
      logical, optional, intent(in) :: unique_only 
    
      integer(i32) :: nequiv, ipt, irot, idiff(dim)
      real(dp) :: tol, rot_point(dim), diff(dim)
      logical :: first, unique 
    
      integer(i32), allocatable :: tmp(:, :)
    
      ! set tolerance
      tol = 1e-12_dp
      if( present( tolerance ) ) tol = tolerance
    
      ! return after first pair was found
      first = .false.
      if( present( first_only ) ) first = first_only

      ! return only unique points
      unique = .false.
      if (present(unique_only)) unique = unique_only
    
      allocate( tmp(npt*nrot, 2+dim), source=-1 )
    
      nequiv = 0
      outer: do irot = 1, nrot
        ! p' = R.p in Cartesian <==> p = R^T.p' in lattice
        ! multiply matrix from right == multiply transpose from left
        rot_point = matmul( point, rotations(:, :, irot) )
        do ipt = 1, npt
          diff =  rot_point - pointlist(:, ipt)
          idiff = nint( diff )
          ! check if points differ by an integer(i32) vector
          if( any( abs( diff - idiff ) > tol ) ) cycle
          ! check if the equivalent point is already mapped
          if( unique .and. any(tmp(:,1) == ipt )) cycle
          ! add to list of equivalent points
          nequiv = nequiv + 1
          tmp(nequiv, 1:2) = [ipt, irot]
          tmp(nequiv, 3:) = idiff
          if( first ) exit outer
        end do
      end do outer
    
      allocate( equiv_point_idx, source=tmp(1:nequiv, 1) )
      allocate( rotation_idx, source=tmp(1:nequiv, 2) )
      if( present( integer_vectors ) ) &
        allocate( integer_vectors, source=transpose( tmp(1:nequiv, 3:) ) )
    
      deallocate( tmp )
    end subroutine find_equivalent_wavevectors


    !> Find the symmetry operation that maps the q-point with index `iq` 
    !> to the q-point with index `iqp`.
    !> Rotation matrices for q-points are constructed using Cartesian rotations, 
    !> as in [spglib](https://github.com/spglib/spglib)
    integer(i32) function what_maps_q_2_qp(iq, iqp, points_list, number_symmetries, symmetry_list, &
                              rotations_cartesian, bvec, binv, qmesh_size, qmesh_offset)

      use math_utils,  only: get_integer_indexes
      use modmpi,      only: terminate_if_false

      !> first q point (FBZ)
      integer(i32), intent(in) :: iq
      !> second q point (FBZ)
      integer(i32), intent(in) :: iqp
      !> reciprocal points in reduced coordinates (FBZ)
      real(dp), intent(in)     :: points_list(:,:)
      !> number of crystal symmetries
      integer(i32), intent(in) :: number_symmetries
      !> Symmetry list from the global set of symmetry operations
      integer(i32), intent(in) :: symmetry_list(:)
      !> Lattice point group symmetries in Cartesian coordinates
      real(dp), intent(in)     :: rotations_cartesian(:,:,:)
      !> Reciprocal lattice coordinates
      real(dp), intent(in)     :: bvec(3,3)
      !> Inverse of the reciprocal lattice coordinates
      real(dp), intent(in)     :: binv(3,3)  
      !> Size of the qmesh
      integer(i32), intent(in) :: qmesh_size(3)
      !> Offset of the qmesh
      real(dp), intent(in)     :: qmesh_offset(3)

      ! The symmetry operation index
      integer(i32) :: isym, symop_id

      ! Logical flag
      logical :: found_symmetry 

      ! rotation matrix in reciprocal
      real(dp) :: rot(3,3)
      ! Points
      real(dp) :: q(3), qp(3), q_rot(3)
      
      ! Get the q point 
      q = points_list(1:3,iq)
      ! Get the qp point 
      qp = points_list(1:3,iqp)

      found_symmetry = .false.

      do isym = 1, number_symmetries

          ! Symmetry operation id in full symmetry operations list
          symop_id = symmetry_list(isym)

          ! Get the rotation in reduced coordinates
          ! We are using the spglib way of symmetry here
          ! See Eqs. 15-18 in https://dx.doi.org/10.1088/1361-648X/acd831
          rot = matmul(binv, matmul(transpose(rotations_cartesian(1:3,1:3, symop_id)), bvec))

          ! Rotate q
          q_rot = matmul(rot, q)

          ! If found save and finish the search
          if ( all(get_integer_indexes(q_rot, qmesh_offset, qmesh_size) == &
                   get_integer_indexes(qp, qmesh_offset, qmesh_size)) )  then
              found_symmetry = .true.
              what_maps_q_2_qp = isym
          end if
      end do

      ! Terminate in case no symmetry operation maps one point to the other
      call terminate_if_false( found_symmetry, "Error(what_maps_q_2_qp): there is no operation mapping q to qp")

    end function what_maps_q_2_qp

    !> Find all q-point pairs equivalent to the input
    !> Be aware that it can produce equivalent pairs.
    !> In the case in which the q and k grids are symmetry-breaking
    !> equivalent will be -1. That means that the point is mapped to a point
    !> that is not present on the mesh.
    pure subroutine get_equivalent_qpairs(iq1, iq2, points_list, equivalent, number_symmetries, symmetry_list, &
                                          rotations_cartesian, bvec, binv, qmesh_size, qmesh_offset)

        use math_utils,  only: get_integer_indexes

        !> first q point (FBZ)
        integer(i32), intent(in)  :: iq1
        !> second q point (FBZ)
        integer(i32), intent(in)  :: iq2
        !> reciprocal points in reduced coordinates (FBZ)
        real(dp), intent(in)      :: points_list(:,:)
        !> number of crystal symmetries
        integer(i32), intent(in)  :: number_symmetries
        !> Symmetry list from the global set of symmetry operations
        integer(i32), intent(in)  :: symmetry_list(:)
        !> Lattice point group symmetries in Cartesian coordinates
        real(dp), intent(in)      :: rotations_cartesian(:,:,:)
        !> Reciprocal lattice coordinates
        real(dp), intent(in)      :: bvec(3,3)
        !> Inverse of the reciprocal lattice coordinates
        real(dp), intent(in)      :: binv(3,3)  
        !> Size of the qmesh
        integer(i32), intent(in)  :: qmesh_size(3)
        !> Offset of the qmesh
        real(dp), intent(in)      :: qmesh_offset(3)
        !> Equivalent q-pairs
        integer(i32), intent(out) :: equivalent(2,number_symmetries)

        ! The symmetry operation index
        integer(i32) :: isym, symop_id

        ! rotation matrix in reciprocal
        real(dp) :: rot(3,3)
        ! Points
        real(dp) :: q1(3), q2(3), q1_rot(3), q2_rot(3)

        ! Index
        integer(i32) :: i

        ! Get the q1 point 
        q1 = points_list(1:3,iq1)
        ! Get the q2 point 
        q2 = points_list(1:3,iq2)

        ! In the case in which the q and k grids are symmetry-breaking
        ! equivalent will be -1
        equivalent(:,:) = -1

        do isym = 1, number_symmetries
            ! Symmetry operation id in symmetry list
            symop_id = symmetry_list(isym)
            ! Get the rotation in reduced coordinates
            ! We are using the spglib way of symmetry here
            ! See Eqs. 15-18 in https://dx.doi.org/10.1088/1361-648X/acd831
            rot = matmul(binv, matmul(transpose(rotations_cartesian(1:3,1:3, symop_id)), bvec))

            ! Rotate qs
            q1_rot = matmul(rot, q1)
            q2_rot = matmul(rot, q2)

            ! Get the new indexes
            do i = 1, size(points_list,2)
                if ( all(get_integer_indexes(q1_rot, qmesh_offset, qmesh_size) == &
                         get_integer_indexes(points_list(1:3,i), qmesh_offset, qmesh_size)) ) equivalent(1,isym) = i
                if ( all(get_integer_indexes(q2_rot, qmesh_offset, qmesh_size) == &
                         get_integer_indexes(points_list(1:3,i), qmesh_offset, qmesh_size)) ) equivalent(2,isym) = i
            end do

        end do

    end subroutine get_equivalent_qpairs

end module mod_symmetry

