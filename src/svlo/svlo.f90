! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2021

! See also:
! Vona, C., Lubeck, S., Kleine, H., Gulans, A. and Draxl, C., 2023.
! Accurate and efficient treatment of spin-orbit coupling via second variation employing local orbitals.
! arXiv preprint arXiv:2306.02965.

!> svlo - second variation that uses local orbitals explicitly as basis functions
module svlo
  use precision, only: dp
  
  implicit none

  private

  !> number of basis functions used in second-variation
  integer :: num_of_basis_funs_sv

  public :: set_num_of_basis_funs_sv
  public :: get_num_of_basis_funs_sv
  public :: is_input_compatible_with_svlo
  public :: construct_H_and_S_in_evecfv_plus_lo_basis

  

contains

  !> Get the number of basis functions used in second-variation
  pure function get_num_of_basis_funs_sv() result(n)
    integer :: n
    n = num_of_basis_funs_sv
  end function get_num_of_basis_funs_sv

  !> Set the number of basis functions used in second-variation
  subroutine set_num_of_basis_funs_sv(n)
    integer, intent(in) :: n
    num_of_basis_funs_sv = n
  end subroutine set_num_of_basis_funs_sv

  !> The purpose of this function is twofold.
  !> On the one hand, it checks if the input parameters that are provided
  !> by the user via the input.xml file are expected to work in combination
  !> with the current svlo implementation.
  !> On the other hand, it provides a To-Do-List for exciting developers.
  !> Developers are expected to work on the compatibility of a feature that is
  !> listed below and then remove the corresponding check from this function.
  ! TODO(Sven) Establish compatibility between svlo and the feautures listed below 
  subroutine is_input_compatible_with_svlo(input, is_compatible_with_svlo, incompatibility_message)
    use modinput, only: input_type
    !> Information from the input file
    type(input_type), intent(in) :: input
    !> True, if the input is compatible with the svlo feature (otherwise false)
    logical, intent(out) :: is_compatible_with_svlo
    !> Incompatibility message for svlo feature
    character(:), allocatable, intent(out) :: incompatibility_message

    character(:), allocatable :: incompatibility_info
    integer :: is, ia

    is_compatible_with_svlo = .true.
    incompatibility_message = 'No incompatibility found.'

    if (associated(input%relax)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/relax.'
    else if (associated(input%properties)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/properties.'
    else if (associated(input%phonons)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/phonons.'
    else if (associated(input%xs)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/xs.'
    else if (associated(input%gw)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/gw.'
    ! not tested but should work with svlo
    else if (associated(input%groundstate%DFTD2parameters)) then 
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/DFTD2parameters.'
    ! not tested but should work with svlo
    else if (associated(input%groundstate%TSvdWparameters)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/TSvdWparameters.'
    else if (associated(input%groundstate%dfthalf)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/dfthalf.'
    else if (associated(input%groundstate%Hybrid)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/Hybrid.'
    else if (associated(input%groundstate%solver)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/solver.'
    else if (associated(input%groundstate%OEP)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/OEP.'
    ! not tested but could work with svlo
    else if (associated(input%groundstate%output)) then 
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/output.'
    else if (associated(input%groundstate%libxc)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/libxc.'
    else if (input%groundstate%ldapu.ne.'none') then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/@ldapu, unless it is equal to "none".'
    else if (input%groundstate%modifiedsv) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/@modifiedsv="true".'
    else if (input%groundstate%tforce) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/@tforce="true".'
    else if (input%groundstate%useDensityMatrix) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/@useDensityMatrix="true".'
    else if (input%groundstate%spin%spinsprl) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/spin/@spinsprl="true".'
    else if (.not.input%groundstate%spin%spinorb) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/spin/@spinorb="false".'
    else if (any(input%groundstate%spin%bfieldc(:) /= 0._dp)) then
       is_compatible_with_svlo = .false.
       incompatibility_info = '/input/groundstate/spin/@bfieldc, unless /input/groundstate/spin/@bfieldc is zero.'
    end if

    if (is_compatible_with_svlo) then
       do is = 1, size(input%structure%speciesarray(:))
          do ia = 1, size(input%structure%speciesarray(is)%species%atomarray(:))
             if (any(input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) /= 0._dp)) then
                is_compatible_with_svlo = .false.
                incompatibility_info = 'structure/species/atom/@bfcmt, unless structure/species/atom/@bfcmt is zero & 
                                       for all atoms.'
             end if
          end do
       end do
    end if

   if (.not. is_compatible_with_svlo) then 
      incompatibility_message = '/input/groundstate/spin/@svlo="true" is currently not compatible with ' // incompatibility_info
   end if 
  end subroutine

  !> First, we construct Hamiltonian and overlap matrices using the APW (augmented plane wave) + lo (local orbitals) basis.
  !> Afterward, these matrices are transformed into a different basis.
  !> This new basis is comprised of two components.
  !> On one hand, we have the first-variational eigenvectors with the local orbital coefficients set to zero,
  !> hereafter referred to as 'evecfv-on-apw'.
  !> On the other hand, we include the original set of local orbitals.
  Subroutine construct_H_and_S_in_evecfv_plus_lo_basis(ik, apwalm, evecfv, igkig, vgkc, ngk, nmat, system_in_evecfv_on_apw_plus_lo)
    Use modfvsystem, only: evsystem, newsystem, deletesystem

    !> k-point index
    Integer, Intent (In) :: ik
    !> APW matching coefficients
    Complex (dp), Intent (In) :: apwalm (:,:,:,:)
    !> first-variational eigenvectors
    Complex (dp), Intent (In) :: evecfv (:,:)
    !> index from G+k-vectors to G-vectors
    Integer, Intent (In) :: igkig(:,:,:)
    !> G+k-vectors in Cartesian coordinates
    Real (dp), Intent (In) :: vgkc(:,:,:,:)
    !> number of G+k-vectors for augmented plane waves
    Integer, Intent (In) :: ngk (:, :)
    !> order of overlap and Hamiltonian matrices for each k-point
    Integer, Intent (In) :: nmat(:, :)
    !> eigensystem that includes hamiltonian and overlap matrix in the evecfv-on-apw + lo basis
    Type (evsystem), Intent (Out) :: system_in_evecfv_on_apw_plus_lo
    
    ! local parameters
    Integer, Parameter :: ispn = 1
    Logical, Parameter :: packed = .false.

    ! local variables
    Type (evsystem) :: system_in_apw_plus_lo
    
    ! initialize new system 
    Call newsystem (system_in_evecfv_on_apw_plus_lo, packed, num_of_basis_funs_sv)
    Call newsystem (system_in_apw_plus_lo, packed, nmat(ispn,ik))

    ! construct hamilton and overlap in old apw + lo basis
    Call hamiltonsetup (system_in_apw_plus_lo, ngk(ispn,ik), apwalm, igkig(1,ispn,ik), vgkc(1,1,ispn,ik))
    Call overlapsetup (system_in_apw_plus_lo, ngk(ispn,ik), apwalm, igkig(1,ispn,ik), vgkc(1,1,ispn,ik))

    ! transform Hamilton to new evecfv-on-apw + lo basis
    Call transform_matrix(nmat(ispn, ik), ngk(ispn,ik), system_in_apw_plus_lo%hamilton%za, system_in_evecfv_on_apw_plus_lo%hamilton%za, evecfv)

    ! transform overlap to new evecfv-on-apw + lo basis
    Call transform_matrix(nmat(ispn, ik), ngk(ispn,ik), system_in_apw_plus_lo%overlap%za, system_in_evecfv_on_apw_plus_lo%overlap%za, evecfv)

    Call deletesystem(system_in_apw_plus_lo)

  End Subroutine construct_H_and_S_in_evecfv_plus_lo_basis

  !> transforms a matrix from apw + lo basis to evecfv-on-apw + lo basis
  Subroutine transform_matrix(nmat, ngk, matrix_in, matrix_out, evecfv)
    Use constants, only: zone, zzero
    Use general_matrix_multiplication, only: matrix_multiply

    !> order of overlap and Hamiltonian matrices
    Integer, Intent (In) :: nmat
    !> number of G+k-vectors for augmented plane waves
    Integer, Intent (In) :: ngk
    !> input matrix in apw + lo basis
    Complex (dp), Intent(In) :: matrix_in(nmat, nmat)
    !> output matrix in evecfv-on-apw + lo basis
    Complex (dp), Intent(Out) :: matrix_out(num_of_basis_funs_sv, num_of_basis_funs_sv)
    !> first-variational eigenvectors
    Complex (dp), Intent(In) :: evecfv(:,:)
    
    ! local variables
    Complex (dp), Allocatable :: matrix__times__evecfv_on_apw(:, :)
    Integer :: icount
    ! array sizes
    Integer :: nmatmax, nstfv
    nmatmax = size( evecfv, dim=1 )
    nstfv = size( evecfv, dim=2 )
    
    Allocate( matrix__times__evecfv_on_apw(nmat, nstfv) )
    
    ! evecfv_on_apw multiplication from right hand side
    Call matrix_multiply( matrix_in(1:nmat,1:ngk), evecfv(1:ngk,1:nstfv), matrix__times__evecfv_on_apw(1:nmat, 1:nstfv) )
    
    ! evecfv@APW - evecfv@APW part
    call matrix_multiply( evecfv(1:ngk,1:nstfv), matrix__times__evecfv_on_apw(1:ngk, 1:nstfv), matrix_out(1:nstfv, 1:nstfv), trans_A = "C" )
    
    ! LO - evecfv@APW part
    matrix_out(nstfv+1:num_of_basis_funs_sv, 1:nstfv) = matrix__times__evecfv_on_apw(ngk+1:nmat, 1:nstfv)

    ! evecfv@APW - LO part
    matrix_out(1:nstfv, nstfv+1:num_of_basis_funs_sv) = &
         Transpose( Conjg( matrix_out(nstfv+1:num_of_basis_funs_sv, 1:nstfv) ) )

    ! LO-LO part
    matrix_out(nstfv+1:num_of_basis_funs_sv, nstfv+1:num_of_basis_funs_sv) = matrix_in(ngk+1:nmat, ngk+1:nmat)

    Deallocate( matrix__times__evecfv_on_apw )

  End Subroutine transform_matrix
end module svlo
