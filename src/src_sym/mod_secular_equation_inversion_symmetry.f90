! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) exciting Code, SOL group. 2025

! Created April 2025 (Mara Voiculescu)

!> Module for taking advantage of inversion symmetry when solving the secular equation.
!>
!> In systems with inversion symmetry, the APWs transform in the following way:
!> \[
!>      \phi_{\textbf{G}+\textbf{k}}(-\textbf{r})=\phi^\ast_{\textbf{G}+\textbf{k}}(\textbf{r}) \: .
!> \]
!> This makes the generally complex Hermitian Hamiltonian and overlap matrices become real symmetric.
!> Since LOs generally do not respect this property, they are attached to so-called fictitious planewaves by
!> constructing a coefficient matrix to transform the Hamiltonian and overlap matrices into the real basis.
!> Diagonalization is then performed in this real representation and the resulting eigenvectors are back-transformed
!> to maintain compatibility with the rest of the code.
!>
!> This module gathers all the necessary code.
module mod_secular_equation_inversion_symmetry
  use precision, only: dp, i32
  use constants, only: zzero, zil
  use modinput, only: input
  use matrix_rank, only: matrix_rank_SVD
  use general_matrix_multiplication, only: matrix_multiply
  use qr_factorization, only: qr_column_pivot
#include "asserts.fpp"
  use modfvsystem, only: evsystem
  use generalized_hermitian_eigenproblem, only: solve_generalized_hermitian_eigenproblem
  use modfvsystem, only: deletesystem
  use modmpi, only: terminate_if_false
  use xstring, only: newline

  implicit none
  private
  public :: get_lo_transformation_matrix_inv_sym, set_lo_transformation_matrix_inv_sym, build_coefficient_matrix_lo, &
       solve_secular_equation_inversion_symmetry, transform_eigenvectors_inversion_symmetry, &
       check_usage_of_inversion_symmetry_solver


  !> Coefficients for attaching the LOs to fictitious planewaves
  complex(dp), allocatable :: lo_transformation_matrix_inv_sym (:,:,:,:)

contains

  function get_lo_transformation_matrix_inv_sym(ispn,ik) result(coeff_matrix)
    integer(i32), intent(in) :: ispn
    integer(i32), intent(in) :: ik
    complex(dp), allocatable :: coeff_matrix(:,:)

    CALL_ASSERT(allocated(lo_transformation_matrix_inv_sym), 'lo_transformation_matrix_inv_sym not allocated')
    coeff_matrix = lo_transformation_matrix_inv_sym(:,:,ispn,ik)
  end function get_lo_transformation_matrix_inv_sym

  subroutine set_lo_transformation_matrix_inv_sym(coeff_matrix)
    complex(dp), intent(in) :: coeff_matrix(:,:,:,:)

    if (allocated(lo_transformation_matrix_inv_sym)) deallocate (lo_transformation_matrix_inv_sym)
    lo_transformation_matrix_inv_sym = coeff_matrix

  end subroutine set_lo_transformation_matrix_inv_sym

  !> This subroutine builds the coefficients used to attach the LOs to so-called fictitious planewaves in order to
  !> transform the Hamiltonian and overlap matrices into the real basis. This is done by forming linear combinations
  !> over the LOs,
  !> \[
  !>     \chi_{\textbf{G}+\textbf{k}}^{\mu} (\textbf{r}) = \sum_{\mu'}\Lambda^{\alpha lm}_{\mathbf{G}+\mathbf{k}}
  !>     \: \phi_{\mu'} (\textbf{r}) \;\delta_{ll_{{\mu}}}\;,
  !> \]
  !> where the sum runs over equivalent LOs of atoms of the same species.
  !> The coefficients are built in the following way:
  !> \[
  !>     \Lambda^{\alpha lm}_{\mathbf{G}+\mathbf{k}} = \mathrm{e}^{\mathrm{i}(\mathbf{G}+\mathbf{k})\mathbf{R}_\alpha}
  !>     \; \mathrm{i}^l \; Y^*_{lm} \left( \widehat{\mathbf{G}+\mathbf{k}} \right) \; ,
  !> \]
  !> the reciprocal lattice vectors \( \textbf{G} \) are selected to ensure that the resulting LO combinations,
  !> and thus the rows of the coefficient matrix, are linearly independent. This is done by choosing one
  !> \( \textbf{G} \)-vector at a time, constructing the coefficients, and checking for linear independence with
  !> respect to the previously selected rows. For this, a singular-value-decomposition (SVD) algorithm is used
  !> (see [[matrix_rank_SVD]]).
  !> The construction of the coefficient matrix is optimized by computing coefficients for the first LO of each
  !> angular momentum \( 𝑙 \) of a given atomic species and copying the values for subsequent LOs of the same \( 𝑙 \).
  !> Afterwards, we ensure that the coefficient matrix is unitary, by using a Householder-based QR-factorization
  !> algorithm (see [[qr_column_pivot]]).
  !>
  !> Note: The construction of the coefficients can be further simplified by summing only over atom pairs related by
  !> inversion instead of all atoms of a certain species.
  !> For more details, see the transformation of the MPB functions in:
  !> M. Betzinger et. al., Phys. Rev. B 81 (19 May 2010), p. 195117. doi: 10.1103/PhysRevB.81.195117
  subroutine build_coefficient_matrix_lo(nspecies, natoms, idxas, lmmaxapw, nlorb, lorbl, nlotot, ngp, idxlm, idxlo, &
       sfacgp, tpgpc, Lambda)
    !> Number of species
    integer(i32), intent(in) :: nspecies
    !> Number of atoms for each species
    integer(i32), intent(in) :: natoms (:)
    !> Map for atoms per species to an atomic index over all atoms in the system
    integer(i32), intent(in) :: idxas (:, :)
    ! (lmaxapw+1)^2
    integer(i32), intent(in) :: lmmaxapw
    !> Number of local-orbitals
    integer(i32), intent(in) :: nlorb (:)
    !> Local-orbital angular momentum
    integer(i32), intent(in) :: lorbl (:, :)
    !> Total number of local-orbitals
    integer(i32), intent(in) :: nlotot
    !> Number of G+p-vectors
    integer(i32), intent(in) :: ngp
    !> Index to (l,m) pairs
    integer(i32), intent(in) :: idxlm (0:input%groundstate%lmaxapw, -input%groundstate%lmaxapw:input%groundstate%lmaxapw)
    !> Index to the position of the local-orbitals in the Hamiltonian and overlap matrices
    integer(i32), intent(in) :: idxlo(:, :, :)
    !> Structure factors of G+p-vectors
    complex(dp), intent(in) :: sfacgp (:, :)
    !> (theta, phi) coordinates of G+p-vectors
    real(dp), intent(in) :: tpgpc (:, :)
    !> Coefficient matrix for transforming LOs
    complex(dp), allocatable, intent(out) :: Lambda(:,:)

    ! Local variables
    complex(dp), allocatable :: R(:,:), Q2(:,:), coeffs(:,:)
    integer(i32), allocatable :: P(:)
    integer(i32) :: l, is, ia, ias, igp, ilo, i1, i0, lm1, lm0, i_g, n_lm, i0_current, i1_current, idxg
    complex(dp), allocatable :: ylmgp(:, :)

    ! Index of first row corresponding to the first local orbital per species with a certain l-value
    integer(i32), allocatable :: i_g_first_lo(:, :)
    ! "ilo" value for first local orbital per species with a certain l-value
    integer(i32), allocatable :: index_first_lo(:, :)

    allocate(Lambda(nlotot, nlotot), source=zzero)
    allocate(coeffs(nlotot, nlotot), source=zzero)
    allocate(ylmgp(lmmaxapw, ngp))

    ! Generate complex spherical harmonics for the G+p-vectors
    do igp = 1, ngp
       call genylm (input%groundstate%lmaxapw, tpgpc(:, igp), ylmgp(:, igp))
    end do

    allocate(i_g_first_lo(0:input%groundstate%lmaxapw, nspecies), source = 0)
    allocate(index_first_lo(0:input%groundstate%lmaxapw, nspecies), source = -1)
    i_g = 1

    do is = 1, nspecies
       do ilo = 1, nlorb(is)
          l = lorbl(ilo, is)
          lm0 = idxlm(l, -l)
          lm1 = idxlm(l, l)
          n_lm = 2 * l + 1

          if (i_g_first_lo(l, is) == 0) then
             ! First LO with a certain l for a given species
             i_g_first_lo(l, is) = i_g
             index_first_lo(l, is) = ilo

             ! Compute coefficients for this l and species
             do igp = 1, ngp
                do ia = 1, natoms(is)
                   ias = idxas(ia, is)
                   i0 = idxlo(lm0, ilo, ias)
                   i1 = idxlo(lm1, ilo, ias)
                   coeffs(i0:i1, i_g) = sfacgp(igp, ias) * zil(l) * conjg(ylmgp(lm0:lm1, igp))
                end do

                ! Check linear independence using SVD
                if (matrix_rank_SVD(coeffs(:, i_g_first_lo(l, is):i_g)) == i_g - i_g_first_lo(l, is) + 1) then
                   i_g = i_g + 1
                   if ((i_g - i_g_first_lo(l, is)) == natoms(is) * n_lm) exit
                end if
             end do

          else
             ! Same l-value
             do ia = 1, natoms(is)
                ias = idxas(ia, is)

                ! LO-indices corresponding to first LO
                i0 = idxlo(lm0, index_first_lo(l, is), ias)
                i1 = idxlo(lm1, index_first_lo(l, is), ias)

                ! Current LO-indices
                i0_current = idxlo(lm0, ilo, ias)
                i1_current = idxlo(lm1, ilo, ias)

                ! Copy coefficients from the first section
                do idxg = 0, natoms(is) * n_lm - 1
                   coeffs(i0_current:i1_current, i_g + idxg) = coeffs(i0:i1, i_g_first_lo(l, is) + idxg)
                end do
             end do
             i_g = i_g + natoms(is) * n_lm
          end if
       end do
    end do

    ! Orthonormalize matrix
    allocate(P(nlotot))
    allocate(Q2(nlotot,nlotot))
    allocate(R(nlotot,nlotot))

    call qr_column_pivot(coeffs, P, Q2, R)
    Lambda = Q2

    deallocate(P)
    deallocate(Q2)
    deallocate(R)

  end subroutine build_coefficient_matrix_lo

  !> This subroutine transforms the Hamiltonian and overlap matrices into real ones, using the LO transformation
  !> matrix.
  subroutine transform_matrix_inversion_symmetry(nmatp, ngp, coeff_matrix, matrix_in, matrix_out)
    !> Order of overlap and Hamiltonian matrices
    integer(i32), intent(in) :: nmatp
    !> Number of G+p-vectors
    integer(i32), intent(in) :: ngp
    !> Coefficient matrix for transforming LOs
    complex(dp), intent(in) :: coeff_matrix(:, :)
    !> Input matrix
    complex(dp), intent(in) :: matrix_in(:, :)
    !> Output matrix
    complex(dp), intent(out) :: matrix_out(:, :)

    ! Local variables
    complex (dp), allocatable :: C(:,:)
    integer(i32) :: nlotot

    nlotot = size(coeff_matrix, dim=1)
    allocate(C(nlotot,nlotot))

    matrix_out(1:ngp,1:ngp) = matrix_in(1:ngp,1:ngp)
    call matrix_multiply(coeff_matrix, matrix_in(ngp+1:nmatp,1:ngp), matrix_out(ngp+1:nmatp,1:ngp), trans_A="C")
    call matrix_multiply(matrix_in(1:ngp,ngp+1:nmatp), coeff_matrix, matrix_out(1:ngp,ngp+1:nmatp))
    call matrix_multiply(coeff_matrix, matrix_in(ngp+1:nmatp,ngp+1:nmatp), C, trans_A="C")
    call matrix_multiply(C, coeff_matrix, matrix_out(ngp+1:nmatp,ngp+1:nmatp))

    deallocate(C)

  end subroutine transform_matrix_inversion_symmetry

  !> This subroutine is used to solve the secular equation in the new real representation, by first transforming the
  !> Hamiltonian and overlap matrices and then employing the LAPACK solver for real symmetric matrices
  !> (see [[solve_generalized_hermitian_eigenproblem]]).
  subroutine solve_secular_equation_inversion_symmetry(system, nmatp, ngp, nstfv,nmatmax, coeff_matrix, eigenvalues, &
       eigenvectors_real)
    !> Eigensystem datastructure containing S and H matrices
    type(evsystem), intent(inout) :: system
    !> Order of overlap and Hamiltonian matrices
    integer(i32), intent(in) :: nmatp
    !> Number of G+p-vectors
    integer(i32), intent(in) :: ngp
    !> Number of first-variational states
    integer(i32), intent(in) :: nstfv
    !> Maximum order of overlap and Hamiltonian matrices
    integer(i32), intent(in) :: nmatmax
    !> Coefficient matrix for transforming LOs
    complex(dp), intent(in) :: coeff_matrix(:, :)
    !> Eigenvalues
    real(dp), intent(out) :: eigenvalues(:)
    !> Eigenvectors of the real secular equation
    real(dp), allocatable, intent(out) :: eigenvectors_real(:, :)

    ! Local variables
    complex(dp), allocatable :: hamilton_transformed(:,:), overlap_transformed(:,:)
    real(dp), allocatable :: hamilton_real(:,:), overlap_real(:,:)
    integer(i32) :: nlotot

    allocate(hamilton_transformed(nmatp,nmatp), source=zzero)
    allocate(overlap_transformed(nmatp,nmatp), source=zzero)

    call transform_matrix_inversion_symmetry(nmatp, ngp, coeff_matrix, system%hamilton%za, hamilton_transformed)
    call transform_matrix_inversion_symmetry(nmatp, ngp, coeff_matrix, system%overlap%za, overlap_transformed)

    allocate(hamilton_real(nmatp,nmatp))
    allocate(overlap_real(nmatp,nmatp))
    allocate(eigenvectors_real(nmatmax, nstfv))

    nlotot = size(coeff_matrix,dim=1)

    overlap_real = real(overlap_transformed, kind=dp)
    hamilton_real = real(hamilton_transformed, kind=dp)

    call solve_generalized_hermitian_eigenproblem(hamilton_real, overlap_real, input%groundstate%solver%evaltol, &
         eigenvalues, eigenvectors_real)

    call deletesystem (system)
    deallocate(hamilton_transformed)
    deallocate(hamilton_real)
    deallocate(overlap_transformed)
    deallocate(overlap_real)

  end subroutine solve_secular_equation_inversion_symmetry

  !> This subroutine back-transforms the eigenvectors from the real basis, in order to make them compatible with the
  !> rest of the code.
  subroutine transform_eigenvectors_inversion_symmetry(ngp, coeff_matrix, eigenvectors_in, eigenvectors_out)
    !> Number of G+p-vectors
    integer(i32), intent(in) :: ngp
    !> Coefficient matrix for transforming LOs
    complex(dp), intent(in) :: coeff_matrix(:, :)
    !> Input eigenvectors in the real basis
    real(dp), intent(in) :: eigenvectors_in(:, :)
    !> Output eigenvectors in the complex basis
    complex(dp), intent(out) :: eigenvectors_out(:, :)

    ! Local variables
    complex(dp), allocatable :: eigenvectors_lo(:,:)
    integer(i32) :: nlotot, nmatp, nstfv

    nlotot = size(coeff_matrix, dim=1)
    nmatp = ngp + nlotot
    nstfv = size(eigenvectors_in, dim=2)

    allocate(eigenvectors_lo(nlotot,nstfv))

    call matrix_multiply(coeff_matrix, eigenvectors_in(ngp+1:nmatp,:), eigenvectors_lo)

    eigenvectors_out(1:ngp,:)=cmplx(eigenvectors_in(1:ngp,:), 0._dp, kind=dp)
    eigenvectors_out(ngp+1:nmatp,:)=eigenvectors_lo(:,:)

  end subroutine transform_eigenvectors_inversion_symmetry

  !> Checks the use of the inversion-symmetry solver as specified in the input file.
  !>
  !> If the inversion-symmetry solver is selected, the code checks that:
  !>   - The structure exhibits inversion symmetry.
  !>   - The inversion center is located at the origin of the unit cell.
  !>
  !> If either condition is not met, the code terminates with a relevant error message.
  !>
  !> Conversely, if the standard solver is used but the structure does possess inversion symmetry with the correct
  !> placement of the inversion center, a warning is written in the WARNINGS.OUT file recommending the use of the
  !> inversion-symmetry solver.
  subroutine check_usage_of_inversion_symmetry_solver(spainvsym, inv_sym_no_translation)
    !> spainvsym is .true. if symmetry group contains spatial inversion symmetry
    logical, intent(in) :: spainvsym
    !> inv_sym_no_translation is .true. if corresponding translation vector of inversion symmetry operation is zero
    logical, intent(in) :: inv_sym_no_translation

    ! Local variables
    character(:), allocatable :: error_message

    if (input%groundstate%solver%type == 'Lapack' .and. spainvsym .and. &
         inv_sym_no_translation) then
       call warning('Inversion symmetry is present. Consider using input%groundstate%solver%type = "inversionsymmetry"')
    elseif (input%groundstate%solver%type == 'inversionsymmetry' .and. .not. spainvsym) then
       error_message = 'System does not possess inversion symmetry.' // newline // &
            'Use input%groundstate%solver%type = "Lapack".'
       call terminate_if_false(.false., error_message)
    elseif (input%groundstate%solver%type == 'inversionsymmetry' .and. spainvsym .and. &
         .not. inv_sym_no_translation) then
       error_message = 'Inversion center is not in the origin of the unit cell.' // newline // &
            'Update the structure manually or use input%groundstate%solver%type = "Lapack".'
       call terminate_if_false(.false., error_message)
    endif

  end subroutine check_usage_of_inversion_symmetry_solver

end module mod_secular_equation_inversion_symmetry
