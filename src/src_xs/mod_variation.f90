! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! !MODULE: modxas
! !DESCRIPTION:
!   Contains additional global variables required for XAS calculations within the BSE EXCITING code.
!
! !REVISION HISTORY:
!
!   Created JUNE 2015 by Christian Vorwerk
module mod_variation
  implicit none
    
contains
  subroutine variation_multiplication(A,B,C,X,dimA,dimB,startA,startB)
!> This subroutine is a wrapper for [[matrix_contraction/contract_A_and_C_with_B_complex_dp]].
!> It configures the slices of the matrices A and C based on user-provided 
!> offsets (startA, startB) and dimensions (dimA, dimB), and then calls
!> contract_A_and_C_with_B_complex_dp with either normal or conjugated C,
!> depending on the global `emat_ccket`.
!>
!> - A is (nstsv x nstsv), potentially containing second-variational states.
!> - B is a matrix with dimension (:,:) containing first-variational momentum matrix elements.
!> - C is (nstsv x nstsv), also second-variational states.
!> - X is the output sub-block of size (dimA x dimB).
!>
!> If `emat_ccket` is true, we conjugate C and apply factor_a=zi, factor_b=-zi.
!> Otherwise, we do a normal call with default factors = 1.0.

    use m_ematqk,                only: emat_ccket
    use constants,               only: zi
    use precision, only: dp
    use mod_eigenvalue_occupancy, only: nstsv
    use svlo, only: get_num_of_basis_functions_sv
    use matrix_contraction,           only: contract_A_and_C_with_B_complex_dp
    implicit none

    integer, intent(in) :: dimA, dimB, startA, startB
    complex(dp), intent(in)  :: A(nstsv, nstsv)
    complex(dp), intent(in)  :: B(:,:)
    complex(dp), intent(in)  :: C(nstsv, nstsv)
    complex(dp), intent(out) :: X(dimA, dimB)

    integer :: endA, endB
    integer :: n_bf_sv
    integer :: i, j 
    
    n_bf_sv = get_num_of_basis_functions_sv()

    ! Determine slice boundaries
    endA = startA + dimA - 1
    endB = startB + dimB - 1
    
    ! emat_ccket-based case separation
    if (.not. emat_ccket) then
       call contract_A_and_C_with_B_complex_dp( &
            A(:, startA:endA),   &
            B,                   &
            C(:, startB:endB),   &
            X                    )
    else
       call contract_A_and_C_with_B_complex_dp( &
            A(:, startA:endA),             &
            B,                             &
            conjg(C(:, startB:endB)),      &
            X,                             &
            factor_a = zi,                &
            factor_b = -zi                )
     end if
   end subroutine variation_multiplication 
!--------------------------------------------------------------------------------
  subroutine getdiffocc(iq, ik, ikq, l1, u1, docc1, docc2)
    ! xssave0 has to be called in advance.
      Use modinput
      Use modmain
      Use modxs
      Use m_genfilname
      use mod_getoccsv, only: getoccsv
      Implicit None
  ! arguments
      Integer, intent (in) :: iq, ik, ikq, l1, u1
      Real (8), intent (out) :: docc1 (u1-l1+1), docc2 (u1-l1+1)
  ! local variables
      Integer :: ist, iqt
      Real (8), Allocatable :: o0 (:), o (:)
      iqt = iq
      Allocate (o0(nstsv), o(nstsv))
  ! eigenvalues and occupancies for k+q-point
      Call getoccsv (vkl(:, ikq), o)
  ! eigenvalues and occupancies for k-point
      Call getoccsv0 (vkl0(1, ik), o0)
  ! loop over band ranges    
    if (input%groundstate%tevecsv) then
      Do ist = l1, u1
        docc1 (ist-l1+1) = o0 (ist) - 1.0d0
        docc2 (ist-l1+1) = o (ist) - 1.d0       
     End Do
    else
      Do ist = l1, u1
        !docc1 (ist-l1+1) = o0 (ist)/2.0d0 - 1.0d0
        !docc2 (ist-l1+1) = o (ist)/2.0d0 - 1.0d0
        docc1 (ist-l1+1) = 1.0d0
        docc2 (ist-l1+1) = 1.0d0       
       
     End Do
    end if
      Deallocate (o0, o)
  end subroutine
 !--------------------------------------------------------------------------------
 !> This subroutine:
 !>  - Finds the ik+q point via \c ikmapikq_ptr
 !>  - Sets up the band ranges for first-variational states
 !>  - Calls \c ematqk to get first-variational matrix elements
 !>  - Retrieves second-variational eigenstates using \c getevecsv0, \c getevecsv1
 !>  - Multiplies them using \c contract_A_and_C_with_B_complex_dp to produce 
 !>    second-variational matrix elements in \a emat
subroutine ematqk_sv(iq, ik, emat, bc)
  use mod_eigenvalue_occupancy, only: nstfv, nstsv
  use modxs, only: bcbs, ngq
  use precision, only: dp
  use mod_ematptr, only: ikmapikq_ptr
  use m_ematqk, only: ematqk,emat_ccket
  use constants, only :zi, zone, zzero
  use m_ematqk, only: ematqk,emat_ccket
  use constants, only :zi, zone, zzero
  use m_getgrst, only: getevecsv0, getevecsv1
  use svlo, only: get_num_of_basis_functions_sv
  use matrix_contraction, only: contract_A_and_C_with_B_complex_dp
  use modinput, only: issvlo
  implicit none

  integer, intent(in) :: iq, ik
  type(bcbs), intent(in) :: bc
  complex(dp), intent(inout) :: emat(:,:,:)  ! final 2nd var results

  type(bcbs) :: bc_
  complex(dp), allocatable :: emat_(:,:,:)
  complex(dp) :: evecsvt0(nstsv,nstsv), evecsvt1(nstsv,nstsv)
  integer :: igq, ikq, num_of_basis_functions_sv
  integer :: endA, endB

  ! find ik+q
  ikq = ikmapikq_ptr(ik, iq)

  num_of_basis_functions_sv = get_num_of_basis_functions_sv()

  ! first var band range
  bc_%n1 = nstfv
  bc_%il1=1
  bc_%iu1=nstfv
  bc_%n2 = nstfv
  bc_%il2=1
  bc_%iu2=nstfv

  ! allocate
  allocate(emat_(num_of_basis_functions_sv,num_of_basis_functions_sv,ngq(iq)),source=zzero)
  call ematqk(iq, ik, emat_, bc_, issvlo())

  ! second variational eigenstates
  call getevecsv0(ik, evecsvt0)
  call getevecsv1(ikq, evecsvt1)

  ! loop over igq
  Do igq=1,ngq(iq)
        call variation_multiplication(evecsvt0,emat_(:,:,igq),evecsvt1,emat(:,:,igq),&
        & bc%n1, bc%n2, bc%il1, bc%il2)
  end Do
      deallocate(emat_)
end subroutine ematqk_sv
 !--------------------------------------------------------------------------------
end module
