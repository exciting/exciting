! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !MODULE: modxas
! !DESCRIPTION:
!   Contains additional global variables required for XAS calculations within the BSE EXCITING code.
!
! !REVISION HISTORY:
!
!   Created JUNE 2015 by Christian Vorwerk
!EOP   
!BOC
#include "maxdefinitions.inc"
module mod_variation
  implicit none
    
contains
!--------------------------------------------------------------------------------
  subroutine variation_multiplication(A,B,C,X,dimA,dimB,startA,startB)
  ! performs multiplication 
  !$$ X_{ij}=\sum_{k=1}^{nstfv}\sum_{l=1}^{nstfv} A_{ki}^{\circ}B_{kl}C_{l,j}+\sum_{k=nstfv+1}^{nstsv}\sum_{l=nstfv+1}^{nstsv} A_{ki}^{\circ}B_{k-nstfv,l-nstfv}C_{l,j}$$
  use m_ematqk, only: emat_ccket
  use mod_constants, only: zi, zone, zzero
  use mod_eigenvalue_occupancy, only: nstfv, nstsv 
  implicit none
  integer, intent (in):: dimA, dimB, startA, startB
  Complex (8), intent (in) :: A(nstsv, nstsv), B(nstfv, nstfv), C(nstsv, nstsv)
  Complex (8), intent (out) :: X(dimA, dimB)
    ! local variables
    Complex (8), dimension(nstfv,dimA) :: inter1, inter3
    Complex (8), dimension(nstfv,dimB) :: inter2, inter4
    Complex (8), dimension(nstfv,dimB) :: xinter
    Integer :: i, j
    Integer :: k, l, ispn, ist, jst
    Complex(8) :: zt, zv
    ! Fill temporary matrices
    if (.not. emat_ccket) then
      do i=1, nstfv
        do j=1, dimA
          inter1(i,j)=A(i,j+startA-1)
          inter3(i,j)=A(i+nstfv,j+startA-1)
        end do
        do j=1, dimB
          inter2(i,j)=C(i,j+startB-1)
          inter4(i,j)=C(i+nstfv,j+startB-1)
        end do
      end do
    else
      do i=1, nstfv
        do j=1, dimA
          inter1(i,j)=A(i,j+startA-1)
          inter3(i,j)=A(i+nstfv,j+startA-1)
        end do
        do j=1, dimB
          inter2(i,j)=zi*conjg(C(i+nstfv,j+startB-1))
          inter4(i,j)=-zi*conjg(C(i,j+startB-1))
        end do
      end do
    end if
    
    call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
               nstfv, &        ! M ... rows of op( A ) = rows of C
               dimB, &         ! N ... cols of op( B ) = cols of C
               nstfv, &        ! K ... cols of op( A ) = rows of op( B )
               zone, &          ! alpha
               B, &         ! A
               nstfv, &          ! LDA ... leading dimension of A
               inter2, &             ! B
               nstfv,&           ! LDB ... leading dimension of B
               zzero, &          ! beta
               xinter, &  ! C
               nstfv & ! LDC ... leading dimension of C
               )
    call zgemm('C', &           ! TRANSA = 'H'  op( A ) = H.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
               dimA, &        ! M ... rows of op( A ) = rows of C
               dimB, &         ! N ... cols of op( B ) = cols of C
               nstfv, &        ! K ... cols of op( A ) = rows of op( B )
               zone, &          ! alpha
               inter1, &         ! A
               nstfv, &          ! LDA ... leading dimension of A
               xinter, &             ! B
               nstfv,&           ! LDB ... leading dimension of B
               zzero, &          ! beta
               X, &  ! C
               dimA & ! LDC ... leading dimension of C
               )
    call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
               nstfv, &        ! M ... rows of op( A ) = rows of C
               dimB, &         ! N ... cols of op( B ) = cols of C
               nstfv, &        ! K ... cols of op( A ) = rows of op( B )
               zone, &          ! alpha
               B, &         ! A
               nstfv, &          ! LDA ... leading dimension of A
               inter4, &             ! B
               nstfv,&           ! LDB ... leading dimension of B
               zzero, &          ! beta
               xinter, &  ! C
               nstfv & ! LDC ... leading dimension of C
               )
    call zgemm('C', &           ! TRANSA = 'H'  op( A ) = H.
               'N', &           ! TRANSB = 'N'  op( B ) = B.
               dimA, &        ! M ... rows of op( A ) = rows of C
               dimB, &         ! N ... cols of op( B ) = cols of C
               nstfv, &        ! K ... cols of op( A ) = rows of op( B )
               zone, &          ! alpha
               inter3, &         ! A
               nstfv, &          ! LDA ... leading dimension of A
               xinter, &             ! B
               nstfv,&           ! LDB ... leading dimension of B
               zone, &          ! beta
               X, &  ! C
               dimA & ! LDC ... leading dimension of C
               )
    return
  end subroutine  

!--------------------------------------------------------------------------------
  subroutine getdiffocc(iq, ik, ikq, l1, u1, docc1, docc2)
    ! xssave0 has to be called in advance.
      Use modinput
      Use modmain
      Use modxs
      Use m_genfilname
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
      Call getoccsv (vkl(1, ikq), o)
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

    subroutine ematqk_sv(iq, ik, emat, bc)
      use mod_eigenvalue_occupancy, only: nstfv, nstsv
      use modxs, only: bcbs, ngq
      use mod_ematptr, only: ikmapikq_ptr
      use m_ematqk, only: ematqk
      use m_getgrst, only: getevecsv0, getevecsv1
      implicit none
          
      ! Arguments
      integer, intent(in) :: iq, ik
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:,:)
      ! local variables
      type(bcbs) :: bc_
      complex(8), allocatable :: emat_(:,:,:)
      complex(8) :: evecsvt0(nstsv,nstsv), evecsvt1(nstsv,nstsv)
      integer :: igq, ikq
      
      ! Find ikq point
      ikq = ikmapikq_ptr(ik, iq)
      
      ! Set intermediate bands to calculate all combinations
      bc_%n1=nstfv
      bc_%il1=1
      bc_%iu1=nstfv
      bc_%n2=nstfv
      bc_%il2=1
      bc_%iu2=nstfv
      !Calculate all first-variational matrix elements
      allocate(emat_(nstfv,nstfv,ngq(iq)))
      call ematqk(iq,ik, emat_, bc_)
      ! Get 2nd variational eigenstates
      Call getevecsv0 (ik, evecsvt0)
      Call getevecsv1 (ikq, evecsvt1)
      ! Calculate 2nd variational matrix elements
      Do igq=1,ngq(iq)
        call variation_multiplication(evecsvt0,emat_(:,:,igq),evecsvt1,emat(:,:,igq),&
          & bc%n1, bc%n2, bc%il1, bc%il2)
      End Do
      deallocate(emat_)
    end subroutine ematqk_sv
 !--------------------------------------------------------------------------------
    
end module
