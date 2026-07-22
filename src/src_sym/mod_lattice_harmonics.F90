! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) exciting Code, SOL group. 2025

! Created February 2025 (Mara Voiculescu)

!> Module for computing lattice harmonics (symmetrized real spherical harmonics).
module mod_lattice_harmonics

  use constants, only: real_zero, sqrt_two, zi, zone, zzero
  use general_matrix_multiplication, only: matrix_multiply
  use grid_utils, only: mesh_1d
  use qr_factorization, only: qr_column_pivot
  use modinput, only: input
  use precision, only: dp, i32
  use vector_multiplication, only: norm

  implicit none

  private
  public :: lattice_harmonics_type, construct_lattice_harmonics_type, &
       set_lattice_harmonics, get_lattice_harmonics, &
       remove_zero_rows, generate_matrix_complex_to_real_spherical_harmonics, &
       construct_lattice_harmonics_coeffs, transform_gaunt_coefficients_lattice_harmonics, &
       transform_mt_potential_lattice_harmonics

  !> Type for storing lattice-harmonics coefficients and number. See [[construct_lattice_harmonics_coeffs]] for a
  !> detailed description of how these arrays are constructed.
  type :: lattice_harmonics_type
     !> Number of lattice harmonics for a specific \( l \)-value and the given index of atoms and species
     integer(i32), allocatable :: number(:,:)
     !> Coefficients for the lattice harmonic expansion, depending on  \( l \), \( m_1 \) and \( m_2 \) values,
     !> as well as the given index of atoms and species
     real(dp), allocatable :: coefficients(:,:,:,:)
     !> Gaunt coefficients in the lattice-harmonic representation, defined as an integral over the product of two
     !> complex spherical harmonics and one lattice harmonic
     complex(dp), allocatable :: gaunt_coefficients(:,:,:,:)
  end type lattice_harmonics_type

  !> Global instance of `[[lattice_harmonics_type]]`.
  type(lattice_harmonics_type) :: global_lattice_harmonics

contains

  !> Setup an instance of `[[lattice_harmonics_type]]`.
  function construct_lattice_harmonics_type(num_lattice_harmonics, lattice_harmonics_coeffs, &
       lattice_harmonics_gaunt_coeffs) result(this)
    !> Number of lattice harmonics for a specific \( l \)-value and the given index of atoms and species
    integer(i32), intent(in) :: num_lattice_harmonics(:,:)
    !> Coefficients for the lattice harmonic expansion, depending on  \( l \), \( m_1 \) and \( m_2 \) values,
    !> as well as the given index of atoms and species
    real(dp), intent(in) :: lattice_harmonics_coeffs(:,:,:,:)
    !> Gaunt coefficients in the lattice-harmonic representation, defined as an integral over the product of two
    !> complex spherical harmonics and one lattice harmonic,
    !> i.e. \(\langle Y_{l_1 m_1} | K_{\nu} | Y_{l_2 m_2} \rangle\).
    complex(dp), intent(in) :: lattice_harmonics_gaunt_coeffs(:,:,:,:)
    type(lattice_harmonics_type) :: this

    this%number = num_lattice_harmonics
    this%coefficients = lattice_harmonics_coeffs
    this%gaunt_coefficients = lattice_harmonics_gaunt_coeffs
  end function construct_lattice_harmonics_type

  !> Set the global `[[lattice_harmonics_type]]` instance
  subroutine set_lattice_harmonics(new_instance)
    type(lattice_harmonics_type), intent(in) :: new_instance
    global_lattice_harmonics = new_instance
  end subroutine set_lattice_harmonics

  !> Retrieve the global `[[lattice_harmonics_type]]` instance
  function get_lattice_harmonics() result(current_instance)
    type(lattice_harmonics_type) :: current_instance
    current_instance = global_lattice_harmonics
  end function get_lattice_harmonics

  !> This subroutine removes zero rows in a given matrix and returns the number of non-zero rows as well as the
  !> reduced matrix.
  subroutine remove_zero_rows(C, num_non_zero_rows, C_reduced)
    !> Input matrix
    real(dp), intent(in) :: C(:, :)
    !> Number of non-zero rows
    integer(i32), intent(out) :: num_non_zero_rows
    !> Reduced matrix
    real(dp), allocatable, intent(out) :: C_reduced(:, :)

    ! Local variables
    real(dp), parameter :: epsilon = 1e-12_dp

    C_reduced = C(pack(mesh_1d(1, size(C, 1)), norm2(C, dim=2) >= epsilon), :)
    num_non_zero_rows = size(C_reduced, 1)

  end subroutine remove_zero_rows

  !> This subroutine computes the matrix which relates complex to real spherical harmonics up to a given maximum angular
  !> momentum quantum number \( l \) :
  !> \[
  !>      S_{lm} = \sum_{m'} A^l_{mm'}Y_{lm} \, .
  !> \]
  subroutine generate_matrix_complex_to_real_spherical_harmonics(lmax, A)
    !> Maximum angular momentum
    integer(i32), intent(in) :: lmax
    !> Matrix which relates complex to real spherical harmonics
    complex(dp), allocatable, intent(out) :: A(:,:,:)

#if !defined(__INTEL_LLVM_COMPILER)
    integer(i32) :: l, m1
#endif
    complex(dp), parameter :: zone_over_sqrt_two = zone / sqrt_two
    complex(dp), parameter :: zi_over_sqrt_two = zi / sqrt_two

    allocate(A(-lmax:lmax, -lmax:lmax, lmax), source = zzero)

    ! N.B.: The code below is valid. However, due to a bug in ifx, an offset is 
    ! needed for correct results. Due to this bug, fill_matrix is used (see issue #240)
    ! TODO: check whenever ifx bug is solved and remove fill_matrix
#if !defined(__INTEL_LLVM_COMPILER)   
    do l = 1, lmax
      A(0, 0, l) = zone
      do m1 = -1, -l, -2
        A(m1, m1, l) = zi_over_sqrt_two
        A(m1, -m1, l) = zi_over_sqrt_two
      end do
      do m1 = -2, -l, -2
        A(m1, m1, l) = zi_over_sqrt_two
        A(m1, -m1, l) = -zi_over_sqrt_two
      end do
      do m1 = 1, l, 2
        A(m1, m1, l) = -zone_over_sqrt_two
        A(m1, -m1, l) = zone_over_sqrt_two
      end do
      do m1 = 2, l, 2
        A(m1, m1, l) = zone_over_sqrt_two
        A(m1, -m1, l) = zone_over_sqrt_two
      end do
    end do
#else
    call fill_matrix( A )
    contains 
      pure subroutine fill_matrix( matrix )
        complex(dp), contiguous, intent(inout) :: matrix(:, :, :)

        integer(i32) :: i, k, k_dim, offset
        
        k_dim = size( A, 3 )
        offset = k_dim + 1
        do k = 1, k_dim
          matrix(offset, offset, k) = zone
          do i = -1, -k, -2
            matrix(i+offset, i+offset, k) = zi_over_sqrt_two
            matrix(i+offset, -i+offset, k) = zi_over_sqrt_two
          end do
          do i = -2, -k, -2
            matrix(i+offset, i+offset, k) = zi_over_sqrt_two
            matrix(i+offset, -i+offset, k) = -zi_over_sqrt_two
          end do
          do i = 1, k, 2
            matrix(i+offset, i+offset, k) = -zone_over_sqrt_two
            matrix(i+offset, -i+offset, k) = zone_over_sqrt_two
          end do
          do i = 2, k, 2
            matrix(i+offset, i+offset, k) = zone_over_sqrt_two
            matrix(i+offset, -i+offset, k) = zone_over_sqrt_two
          end do
        end do
      end subroutine
#endif      
  end subroutine generate_matrix_complex_to_real_spherical_harmonics

  !> This subroutine generates the coefficients used to construct the lattice harmonics as a linear combination of real
  !> spherical harmonics.
  !> \[
  !>      K^{\alpha}_{\nu}(\hat{\textbf{r}}_{\alpha}) = \sum_{m = - l_{\nu}}^{l_{\nu}} C^{\alpha}_{m\nu} \;
  !>      S_{l_{\nu}m}(\hat{\textbf{r}}_{\alpha}) \;,
  !> \]
  !> where \( K^{\alpha}_{\nu} \) are the lattice harmonics, \( C^{\alpha}_{m\nu} \) are the corresponding coefficients
  !> and \( S_{lm} \) are the real spherical harmonics.
  !>
  !> the coefficients  are determined by imposing that the lattice harmonics are invariant under all the operations of
  !> the site symmetry group of a given atom. This is ensured by applying the sum over all the operations in the site
  !> symmetry group to the real spherical harmonics. Therefore, we generate rotation matrices for real spherical
  !> harmonics.
  !> \[
  !>      \Delta^l_{m'm}(\mathbfcal{R}) = \sum_{m'' m'''} A^l_{m'm''} \, D^l_{m''m'''}(\mathbfcal{R}) \,
  !>      \bigl(A^l_{m'''m}\bigr)^{-1} \; ,
  !> \]
  !> where \( \mathbfcal{R} \) is the rotation matrix in cartesian coordinates, \( \textbf{A}^l = [A^l_{mm'}]\) is the
  !> matrix that relates real and complex spherical harmonics
  !> \[
  !>      S_{lm} = \sum_{m'} A^l_{mm'}Y_{lm}
  !> \]
  !> and \( \textbf{D}^l(\mathbfcal{R}) = [D^l_{mm'}(\mathbfcal{R})]\) is the rotation matrix of complex spherical
  !> harmonics,
  !> also known as Wigner <i>D</i> matrix (see [[getdlmm(function)]]).
  !> Then, the coefficients are then constructed as
  !> \[
  !>      C^l_{m'm} = \sum_{\mathbfcal{R}}\Delta^l_{m'm}(\mathbfcal{R}) \; .
  !> \]
  !>
  !> The coefficients are then orthonormalized using the QR factorization method with pivoting
  !> (see [[qr_column_pivot]]). Then linearly dependent and zero norm vectors are discarded.
  subroutine construct_lattice_harmonics_coeffs(natmtot, nsymsite, symlatc, lsplsyms, coeffs, num_non_zero_rows)
    !> Total number of atoms
    integer(i32), intent(in) :: natmtot
    !> Number of site symmetries per atoms and species
    integer(i32), intent(in) :: nsymsite(:)
    !> Bravais lattice point group symmetries in cartesian coordinates
    real(dp), intent(in) :: symlatc(:, :, :)
    !> Site symmetry spatial rotation element in lattice point group
    integer(i32), intent(in) :: lsplsyms(:, :)
    !> Coefficients for the lattice harmonic expansion
    real(dp), allocatable, intent(out) :: coeffs(:, :, :,:)
    !> Number of lattice harmonics
    integer(i32), allocatable, intent(out) :: num_non_zero_rows(:,:)

    !> Fortran function
    complex(dp), external :: getdlmm

    ! Local variables
    integer :: lmax, l, m1, m2, is, ia, ias, isym, lspl, lmaxvr, lmaxinr, lmaxmax, dim_m, i ,j, lmmax
    real(dp) :: sym_op_matrix(3,3), inv_sym_op_matrix(3,3)
    complex(dp), allocatable :: C_real_sph_harm(:,:,:,:), C_complex_sph_harm(:,:,:,:), C_temp(:,:,:,:), A(:,:,:,:), &
         A_l(:,:,:)
    integer(i32), allocatable :: P(:), num_non_zero_rows_before_qr(:,:)
    real(dp), allocatable :: C_reduced_before_qr(:, :, :,:), reduced_matrix(:,:), C_full(:,:,:,:), B(:,:), R(:,:), &
         Q(:,:)
    real(dp), parameter :: epsilon = 1e-12_dp

    lmaxvr = input%groundstate%lmaxvr
    lmaxinr = input%groundstate%lmaxinr
    lmaxmax = max(lmaxvr, lmaxinr)

    allocate(C_real_sph_harm(-lmaxmax:lmaxmax, -lmaxmax:lmaxmax, lmaxmax, natmtot), source = zzero)
    allocate(C_full(-lmaxmax:lmaxmax, -lmaxmax:lmaxmax, lmaxmax, natmtot), source = real_zero)
    allocate(C_complex_sph_harm(-lmaxmax:lmaxmax, -lmaxmax:lmaxmax, lmaxmax, natmtot), source = zzero)
    allocate(C_temp(-lmaxmax:lmaxmax, -lmaxmax:lmaxmax, lmaxmax, natmtot), source = zzero)
    allocate(A(-lmaxmax:lmaxmax, -lmaxmax:lmaxmax, lmaxmax, natmtot), source = zzero)
    allocate(C_reduced_before_qr(2*lmaxmax+1, 2*lmaxmax+1, lmaxmax, natmtot), source = real_zero)
    allocate(B(2*lmaxmax+1, 2*lmaxmax+1))
    allocate(R(2*lmaxmax+1, 2*lmaxmax+1))
    allocate(Q(2*lmaxmax+1, 2*lmaxmax+1))
    allocate(P(2*lmaxmax+1))

    ! Generate coefficients by applying site-symmetry operations to complex spherical harmonics.
    lmax = lmaxmax
    do ias = 1, natmtot
       do isym = 1, nsymsite(ias)
          lspl = lsplsyms (isym, ias)
          sym_op_matrix = symlatc(:,:,lspl)
          call r3minv (sym_op_matrix, inv_sym_op_matrix)
          do l = 1, lmax
             do m2 = -l, l
                do m1 = -l, l
                   C_complex_sph_harm(m1, m2, l, ias) = C_complex_sph_harm(m1, m2, l, ias) + &
                        getdlmm(inv_sym_op_matrix, l, m1, m2)
                end do
             end do
          end do
       end do
    end do

    ! Set values lower than epsilon to zero.
    C_complex_sph_harm = merge(C_complex_sph_harm, zzero, abs(C_complex_sph_harm) >= epsilon)

    ! Build A-matrix which relates real and complex spherical harmonics.
    call generate_matrix_complex_to_real_spherical_harmonics(lmax, A_l)

    ! Multiply complex coefficient with A-matrix in order to get coefficients for real spherical harmonics.
    do ias = 1, natmtot
       do l = 1, lmax
          call matrix_multiply(A_l(:, :, l), C_complex_sph_harm(:, :, l, ias), C_temp(:, :, l, ias))
          call matrix_multiply(C_temp(:, :, l, ias), A_l(:, :, l), C_real_sph_harm(:, :, l, ias), trans_B = "C")
       end do
    end do

    ! Full coefficient matrix for real spherical harmonics before removal of non-zero and linearly dependent rows.
    C_full = real(C_real_sph_harm)

    allocate(num_non_zero_rows_before_qr(lmax, natmtot))

    ! Remove zero rows from coefficient matrix.
    do ias = 1, natmtot
       do l = 1, lmax
          call remove_zero_rows(C_full(-l:l, -l:l, l, ias), num_non_zero_rows_before_qr(l, ias), reduced_matrix)
          C_reduced_before_qr(1:num_non_zero_rows_before_qr(l, ias), 1:2*l+1, l, ias) = reduced_matrix
       end do
    end do

    ! Orthonormalize coefficient matrix and set linearly dependent vectors to zero. The coefficient matrix is transposed
    ! because the QR routine orthonormalizes matrices column-wise and for our purpose the rows need to be
    ! orthonormalized.
    do ias = 1, natmtot
       do l = 1, lmax
          dim_m = 2*l + 1
          B = transpose(C_reduced_before_qr(1:num_non_zero_rows_before_qr(l, ias), 1:dim_m, l, ias))
          call qr_column_pivot(B, P(1:num_non_zero_rows_before_qr(l, ias)), Q(1:dim_m, 1:dim_m), &
               R(1:dim_m, 1:num_non_zero_rows_before_qr(l, ias)))
          do i = 1, num_non_zero_rows_before_qr(l, ias)
             if (abs(R(i, i)) < epsilon) then
                do j = 1, dim_m
                   Q(j ,i) = zzero
                end do
             end if
          end do
          C_reduced_before_qr(1:num_non_zero_rows_before_qr(l, ias), 1:dim_m, l, ias) = &
               transpose(Q(1:dim_m, 1:num_non_zero_rows_before_qr(l, ias)))
       end do
    end do

    allocate(coeffs(maxval(num_non_zero_rows_before_qr), 2*lmaxmax+1, lmaxmax, natmtot), source = real_zero)
    allocate(num_non_zero_rows(lmax, natmtot))

    ! Remove zero rows again after orthonormalization. The result is the final coefficents for lattice harmonics.
    do ias = 1, natmtot
       do l = 1, lmax
          call remove_zero_rows(C_reduced_before_qr(1:num_non_zero_rows_before_qr(l,ias), 1:2*l+1,l, ias), &
               num_non_zero_rows(l,ias), reduced_matrix)
          coeffs(1:num_non_zero_rows(l,ias), 1:2*l+1,l, ias) = reduced_matrix
       end do
    end do

  end subroutine construct_lattice_harmonics_coeffs

  !> Generate index map from an \((l, m)\)-pair to the corresponding lattice-harmonic index \(\nu\).
  !> For each angular momentum quantum number \(l\), \(\nu\) corresponds to the sequential numbering starting from
  !> `idxlm(l, -l)` up to the number of lattice harmonics for the specific \( l \)-value and the given index of
  !> atoms and species.
  subroutine generate_index_map_lm_to_nu(num_lattice_harmonics, idxlm, lmax, idx_nu)
    !> Number of lattice harmonics for given index of atoms and species
    integer(i32), intent(in) :: num_lattice_harmonics(:)
    !> Index to (l,m) pairs
    integer(i32), intent(in) :: idxlm (0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw)
    !> Maximum angular momentum
    integer(i32), intent(in) :: lmax
    !> Map from (l,m) pairs to indices of lattice harmonics
    integer(i32), allocatable, intent(out) :: idx_nu(:)

    ! Local variables
    integer :: i_nu, l, m, i

    if (allocated(idx_nu)) deallocate(idx_nu)
    allocate (idx_nu(sum(num_lattice_harmonics)+1))

    idx_nu(1) = 1
    i_nu = 2
    do l = 1, lmax
       if (num_lattice_harmonics(l) /= 0) then
          idx_nu(i_nu:i_nu + num_lattice_harmonics(l) - 1) = [(idxlm(l, -l) + i - 1, i = 1, num_lattice_harmonics(l))]
          i_nu = i_nu + num_lattice_harmonics(l)
       end if
    end do

  end subroutine generate_index_map_lm_to_nu

  !> This routine switches the representation of a function from a real spherical-harmonic expansion
  !> to a lattice-harmonic expansion:
  !> \[
  !>      f(\mathbf{r}) = \sum_{lm} f_{lm}(r) \, S_{lm}(\hat{\mathbf{r}}) =
  !>     \sum_{\nu} f_{\nu}(r) \, K_{\nu}(\hat{\mathbf{r}}) \;.
  !> \]
  !> The transformation of the expansion coefficients is given by:
  !> \[
  !>    f_{\nu}(r) = \sum_{m=-l_\nu}^{l_\nu} \! \Omega_{m\nu} \, f_{l_\nu m}(r) \;,
  !> \]
  !> where \(\Omega_{m\nu}\) are the lattice-harmonic transformation coefficients
  !> generated by [[construct_lattice_harmonics_coeffs]], for a given index of atoms and species. The indices \(\nu\)
  !> are given by the routine [[generate_index_map_lm_to_nu]].
  subroutine transform_real_expansion_real_sph_to_lat_harm(rflm, coeffs, num_lattice_harmonics, idxlm, &
       lmax, idx_nu, kflm_reduced)
    !> Coefficients of real-spherical-harmonic expansion
    real(dp), intent(in) :: rflm(:)
    !> Lattice harmonics coefficients for given index of atoms and species
    real(dp), intent(in) :: coeffs(:,:,:)
    !> Number of lattice harmonics for given index of atoms and species
    integer(i32), intent(in) :: num_lattice_harmonics(:)
    !> Index to (l,m) pairs
    integer(i32), intent(in) :: idxlm (0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw)
    !> Maximum angular momentum
    integer(i32), intent(in) :: lmax
    !> Map from (l,m) pairs to indices of lattice harmonics
    integer(i32), intent(in) :: idx_nu(:)
    !> Reduced array of lattice-harmonic expansion coefficients
    real(dp), intent(out) :: kflm_reduced(:)

    ! Local variables
    integer :: l, i_nu
    real(dp) :: kflm(size(rflm))

    kflm = real_zero

    kflm(1) = rflm(1)
    do l = 1, lmax
       if (num_lattice_harmonics(l) /= 0) then
          call matrix_multiply(coeffs(1:num_lattice_harmonics(l), 1:2*l+1, l), rflm(idxlm(l, -l):idxlm(l, l)), &
               kflm(idxlm(l, -l):idxlm(l, -l) + num_lattice_harmonics(l) - 1))
       end if
    end do

    do i_nu = 1, size(idx_nu)
       kflm_reduced(i_nu) = kflm(idx_nu(i_nu))
    end do

  end subroutine transform_real_expansion_real_sph_to_lat_harm

  !> This routine is analogous to [[transform_real_expansion_real_sph_to_lat_harm]], but for complex-valued functions
  !> based on real spherical harmonics (e.g., Gaunt coefficients `gaunt_yry` in the `wigner3j_symbol` module).
  subroutine transform_complex_expansion_real_sph_to_lat_harm(zflm, coeffs, num_lattice_harmonics, idxlm, &
       lmax, idx_nu, kflm_reduced)
    !> Coefficients of spherical-harmonic expansion
    complex(dp), intent(in) :: zflm(:)
    !> Lattice harmonics coefficients for given index of atoms and species
    real(dp), intent(in) :: coeffs(:,:,:)
    !> Number of lattice harmonics for given index of atoms and species
    integer(i32), intent(in) :: num_lattice_harmonics(:)
    !> Index to (l,m) pairs
    integer(i32), intent(in) :: idxlm (0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw)
    !> Maximum angular momentum
    integer(i32), intent(in) :: lmax
    !> Map from (l,m) pairs to indices of lattice harmonics
    integer(i32), intent(in) :: idx_nu(:)
    !> Reduced array of lattice-harmonic expansion coefficients
    complex(dp), intent(out) :: kflm_reduced(:)

    ! Local variables
    integer :: l, i_nu
    complex(dp) :: kflm(size(zflm))

    kflm = zzero

    kflm(1) = zflm(1)
    do l = 1, lmax
       if (num_lattice_harmonics(l) /= 0) then
          call matrix_multiply(coeffs(1:num_lattice_harmonics(l), 1:2*l+1, l), zflm(idxlm(l, -l):idxlm(l, l)), &
               kflm(idxlm(l, -l):idxlm(l, -l) + num_lattice_harmonics(l) - 1))
       end if
    end do

    do i_nu = 1, size(idx_nu)
       kflm_reduced(i_nu) = kflm(idx_nu(i_nu))
    end do

  end subroutine transform_complex_expansion_real_sph_to_lat_harm

  !> This routine transforms the MT potential from the real spherical-harmonic representation to the lattice-harmonic
  !> representation (for more details see [[transform_real_expansion_real_sph_to_lat_harm]]).
  subroutine transform_mt_potential_lattice_harmonics(mt_pot, nspecies, natoms, idxas, nrmt, nrmtmax, natmtot, coeffs, &
       num_lattice_harmonics, idxlm, mt_pot_lh)
    !> MT potential coefficients in the standard representation
    real(dp), intent(in) :: mt_pot(:,:,:)
    !> Number of species
    integer(i32), intent(in) :: nspecies
    !> Number of atoms for each species
    integer(i32), intent(in) :: natoms(:)
    !> Map for atoms per species to an atomic index over all atoms in the system
    integer(i32), intent(in) :: idxas(:, :)
    !> Number of muffin-tin radial points for each species
    integer(i32), intent(in) :: nrmt(:)
    !> Maximum nrmt over all the species
    integer(i32), intent(in) :: nrmtmax
    !> Total number of atoms
    integer(i32), intent(in) :: natmtot
    !> Lattice harmonics coefficients
    real(dp), intent(in) :: coeffs(:,:,:,:)
    !> Number of lattice harmonics
    integer(i32), intent(in) :: num_lattice_harmonics(:,:)
    !> Index to (l,m) pairs
    integer(i32), intent(in) :: idxlm (0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw)
    !> MT potential coefficients in the lattice-harmonic representation
    real(dp), intent(out) :: mt_pot_lh(:,:,:)

    ! Local variables
    integer :: is, ia, nr, ir, ias, lmaxvr
    integer(i32), allocatable :: idx_nu(:)

    lmaxvr = input%groundstate%lmaxvr

    do is = 1, nspecies
       nr = nrmt (is)
       do ia = 1, natoms (is)
          ias = idxas (ia, is)
          call generate_index_map_lm_to_nu(num_lattice_harmonics(:,ias), idxlm, lmaxvr, idx_nu)
          do ir = 1, nr
             call transform_real_expansion_real_sph_to_lat_harm(mt_pot(:, ir, ias), &
                  coeffs(:,:,:,ias), num_lattice_harmonics(:, ias), idxlm, lmaxvr, idx_nu, mt_pot_lh(:, ir, ias))
          end do
       end do
    end do

  end subroutine transform_mt_potential_lattice_harmonics

  !> This routine computes Gaunt coefficients of a lattice harmonic and two complex spherical harmonics.
  !> They are determined as
  !> \[
  !>      G^{\nu}_{l_{\lambda}, m_{\lambda}, l_{\lambda'}, m_{\lambda'}} =
  !>     \sum_{m=-l_\nu}^{l_\nu} \! \Omega_{m\nu} \; G^{l_\nu m}_{l_{\lambda}, m_{\lambda}, l_{\lambda'}, m_{\lambda'}},
  !> \]
  !> where \(\Omega_{m\nu}\) are the lattice-harmonic transformation coefficients generated by
  !> [[construct_lattice_harmonics_coeffs]] and \( G^{l_\nu m}_{l_{\lambda}, m_{\lambda}, l_{\lambda'}, m_{\lambda'}} \)
  !> are Gaunt coefficients of a real spherical harmonic and two complex spherical harmonics (see [[wigner3j_symbol]]).
  subroutine transform_gaunt_coefficients_lattice_harmonics(gnt_ryy, natmtot, coeffs, num_lattice_harmonics, idxlm, &
       gnt_kyy)
    !> Gaunt coefficient of a real spherical harmonic and two complex spherical harmonics
    complex(dp), intent(in) :: gnt_ryy(:,:,:)
    !> Total number of atoms
    integer(i32), intent(in) :: natmtot
    !> Lattice harmonics coefficients
    real(dp), intent(in) :: coeffs(:,:,:,:)
    !> Number of lattice harmonics
    integer(i32), intent(in) :: num_lattice_harmonics(:,:)
    !> Index to (l,m) pairs
    integer(i32), intent(in) :: idxlm (0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw)
    !> Gaunt coefficient of a lattice harmonic and two complex spherical harmonics
    complex(dp), allocatable, intent(out) :: gnt_kyy(:,:,:,:)

    ! Local variables
    integer :: l1, m1, lm1, l3, m3, lm3, ias, lmaxapw, lmaxvr, lmmaxapw, max_num_lattice_harmonics
    integer(i32), allocatable :: idx_nu(:)

    lmaxapw = input%groundstate%lmaxapw
    lmaxvr = input%groundstate%lmaxvr
    lmmaxapw = (lmaxapw + 1) ** 2

    max_num_lattice_harmonics = maxval(sum(num_lattice_harmonics, dim=1))
    allocate (gnt_kyy(natmtot, max_num_lattice_harmonics + 1, lmmaxapw, lmmaxapw), source=zzero)

    do l1 = 0, lmaxapw
       do m1 = - l1, l1
          lm1 = idxlm (l1, m1)
          do l3 = 0, lmaxapw
             do m3 = - l3, l3
                lm3 = idxlm (l3, m3)
                do ias = 1, natmtot
                   call generate_index_map_lm_to_nu(num_lattice_harmonics(:,ias), idxlm, lmaxvr, idx_nu)
                   call transform_complex_expansion_real_sph_to_lat_harm(gnt_ryy(:, lm3, lm1), coeffs(:,:,:,ias), &
                        num_lattice_harmonics(:,ias), idxlm, lmaxvr, idx_nu, gnt_kyy(ias, :, lm3, lm1))
                end do
             end do
          end do
       end do
    end do

  end subroutine transform_gaunt_coefficients_lattice_harmonics

end module mod_lattice_harmonics
