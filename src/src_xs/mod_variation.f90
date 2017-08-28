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
  use m_b_ematqk, only: emat_ccket
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
  subroutine b_wavefmtsv1(lrstp ,lmax ,is ,ia , ngp, ist1, apwalm, evecfv, evecsv, wfmtsv)

    use modinput
    use m_b_getgrst, only: b_wavefmt1
    use mod_muffin_tin, only: nrmt, lmmaxapw, nrmtmax
    use mod_atoms, only: idxas, natmtot
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use mod_ematptr, only: ngkmax1_ptr
    use mod_APW_LO, only: apwordmax
    use mod_eigensystem, only: nmatmax
    use mod_spin, only: nspinor
    !use modmain
    implicit none

    ! arguments
    integer, intent (in) :: lrstp
    integer, intent (in) :: lmax
    integer, intent (in) :: is
    integer, intent (in) :: ia
    integer, intent (in) :: ngp
    integer,    intent(in)  :: ist1
    complex(8), intent(in)  :: apwalm(ngkmax1_ptr,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: wfmtsv(lmmaxapw,nrmtmax,nspinor)
    ! local variables
    integer    :: ispn, ias
    integer    :: i, j, n, ist, igp
    complex(8) :: zt1
    ! allocatable arrays
    logical,    allocatable :: done (:)
    complex(8), allocatable :: wfmt1(:,:,:)

    wfmtsv(:,:,:) = 0.d0

    allocate(done(nstfv))
    allocate (wfmt1(lmmaxapw,nrmtmax,nstsv))

    wfmt1(:,:,:)=0.0d0
    n = lmmaxapw*nrmt(is)
    ias = idxas(ia,is)
    done(:) = .false.
    ! generate spinor wavefunction from second-variational eigenvectors
    i = 0
    do ispn = 1, nspinor
      do ist = 1, nstfv
        i = i + 1
        zt1 = evecsv(i,ist1)
        if (abs(zt1)>1.d-8) then
          if (.not.done(ist)) then
            call b_wavefmt1(lrstp, lmax, is, ia, &
              & ngp, apwalm, evecfv(:,ist), lmmaxapw, wfmt1(:,:,ist))
            done(ist) = .true.
          end if
          ! add to spinor wavefunction
          call zaxpy(n, zt1, wfmt1(:,:,ist), 1, &
          &          wfmtsv(:,:,ispn), 1)
        end if
      end do ! ist
    end do ! ispn
    deallocate(done, wfmt1)
    return
 end subroutine b_wavefmtsv1
!--------------------------------------------------------------------------------
  subroutine b_wavefmtsv0(lrstp ,lmax ,is ,ia , ngp, ist1, apwalm, evecfv, evecsv, wfmtsv)

    use modinput
    use m_b_getgrst, only: b_wavefmt0
    use mod_muffin_tin, only: nrmt, lmmaxapw, nrmtmax
    use mod_atoms, only: idxas, natmtot
    use mod_eigenvalue_occupancy, only: nstfv, nstsv
    use mod_ematptr, only: ngkmax0_ptr
    use mod_APW_LO, only: apwordmax
    use mod_eigensystem, only: nmatmax
    use mod_spin, only: nspinor
    !use modmain
    implicit none

    ! arguments
    integer, intent (in) :: lrstp
    integer, intent (in) :: lmax
    integer, intent (in) :: is
    integer, intent (in) :: ia
    integer, intent (in) :: ngp
    integer,    intent(in)  :: ist1
    complex(8), intent(in)  :: apwalm(ngkmax0_ptr,apwordmax,lmmaxapw,natmtot)
    complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
    complex(8), intent(in)  :: evecsv(nstsv,nstsv)
    complex(8), intent(out) :: wfmtsv(lmmaxapw,nrmtmax,nspinor)
    ! local variables
    integer    :: ispn, ias
    integer    :: i, j, n, ist, igp
    complex(8) :: zt1
    ! allocatable arrays
    logical,    allocatable :: done (:)
    complex(8), allocatable :: wfmt1(:,:,:)

    wfmtsv(:,:,:) = 0.d0

    allocate(done(nstfv))
    allocate (wfmt1(lmmaxapw,nrmtmax,nstsv))

    wfmt1(:,:,:)=0.0d0
    n = lmmaxapw*nrmt(is)
    ias = idxas(ia,is)
    done(:) = .false.
    ! generate spinor wavefunction from second-variational eigenvectors
    i = 0
    do ispn = 1, nspinor
      do ist = 1, nstfv
        i = i + 1
        zt1 = evecsv(i,ist1)
        if (abs(zt1)>1.d-8) then
          if (.not.done(ist)) then
            call b_wavefmt0(lrstp, lmax, is, ia, &
              & ngp, apwalm, evecfv(:,ist), lmmaxapw, wfmt1(:,:,ist))
            done(ist) = .true.
          end if
          ! add to spinor wavefunction
          call zaxpy(n, zt1, wfmt1(:,:,ist), 1, &
          &          wfmtsv(:,:,ispn), 1)
        end if
      end do ! ist
    end do ! ispn
    deallocate(done, wfmt1)
    return
 end subroutine b_wavefmtsv0
 !--------------------------------------------------------------------------------

    subroutine getmoo_sv(iq, ik, emat, bc)
      use mod_eigenvalue_occupancy, only: nstfv, nstsv
      use modxs, only: bcbs, ngq
      use mod_ematptr, only: vkl0_ptr
      use m_b_ematqk, only: b_ematqk
      implicit none
          
      ! Arguments
      integer, intent(in) :: iq, ik
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:,:)
      ! local variables
      type(bcbs) :: emat_
      complex(8), allocatable :: moo_(:,:,:)
      complex(8) :: evecsvt0(nstsv,nstsv), evecsvt1(nstsv,nstsv)
      integer :: igq
      ! Set intermediate bands to calculate all combinations
      emat_%n1=nstfv
      emat_%il1=1
      emat_%iu1=nstfv
      emat_%n2=nstfv
      emat_%il2=1
      emat_%iu2=nstfv
      !Calculate all first-variational matrix elements
      allocate(moo_(nstfv,nstfv,ngq(iq)))
      call b_ematqk(iq,ik, moo_, emat_)
      ! Get 2nd variational eigenstates
      Call getevecsv (vkl0_ptr(1, ik), evecsvt1)
      Call getevecsv (vkl0_ptr(1, ik), evecsvt0)
      ! Calculate 2nd variational matrix elements
      Do igq=1,ngq(iq)
        call variation_multiplication(evecsvt1,moo_(:,:,igq),evecsvt0,emat(:,:,igq),&
          & bc%n1, bc%n2, bc%il1, bc%il2)
      End Do
      deallocate(moo_)
    end subroutine getmoo_sv
 !--------------------------------------------------------------------------------
    
    subroutine getmuu_sv(iq, ik, emat, bc)
      use mod_eigenvalue_occupancy, only: nstfv, nstsv
      use modxs, only: bcbs
      use mod_ematptr, only: vkl1_ptr, ikmapikq_ptr
      use modxs, only: ngq, bcbs
      use m_b_ematqk, only: b_ematqk

      implicit none
          
      ! Arguments
      integer, intent(in) :: iq, ik
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:,:)
      ! local variables
      type(bcbs) :: emat_
      complex(8), allocatable :: muu_(:,:,:)
      complex(8) :: evecsvt0(nstsv,nstsv), evecsvt1(nstsv,nstsv)
      integer :: igq, ikq

      ! Find ikq point
      ikq = ikmapikq_ptr(ik, iq)
      ! Set intermediate bands to calculate all combinations
      emat_%n1=nstfv
      emat_%il1=1
      emat_%iu1=nstfv
      emat_%n2=nstfv
      emat_%il2=1
      emat_%iu2=nstfv
      !Calculate all first-variational matrix elements
      allocate(muu_(nstfv,nstfv,ngq(iq)))
      call b_ematqk(iq,ik, muu_, emat_)
      ! Get 2nd variational eigenstates
      Call getevecsv (vkl1_ptr(1, ikq), evecsvt1)
      Call getevecsv (vkl1_ptr(1, ikq), evecsvt0)
      ! Calculate 2nd variational matrix elements
      Do igq=1,ngq(iq)
        call variation_multiplication(evecsvt1,muu_(:,:,igq),evecsvt0,emat(:,:,igq),&
          & bc%n1, bc%n2, bc%il1, bc%il2)
      End Do
      deallocate(muu_)
    end subroutine getmuu_sv
 !--------------------------------------------------------------------------------
    
    subroutine getmou_sv(iq, ik, emat, bc)
      use mod_eigenvalue_occupancy, only: nstfv, nstsv
      use modxs, only: bcbs
      use mod_ematptr, only: vkl0_ptr, vkl1_ptr, ikmapikq_ptr
      use modxs, only: ngq, bcbs
      use m_b_ematqk, only: b_ematqk



      implicit none
          
      ! Arguments
      integer, intent(in) :: iq, ik
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:,:)
      ! local variables
      type(bcbs) :: emat_
      complex(8), allocatable :: mou_(:,:,:)
      complex(8) :: evecsvt0(nstsv,nstsv), evecsvt1(nstsv,nstsv)
      integer :: igq, ikq

      ! Find ikq point
      ikq = ikmapikq_ptr(ik, iq)
      ! Set intermediate bands to calculate all combinations
      emat_%n1=nstfv
      emat_%il1=1
      emat_%iu1=nstfv
      emat_%n2=nstfv
      emat_%il2=1
      emat_%iu2=nstfv
      !Calculate all first-variational matrix elements
      allocate(mou_(nstfv,nstfv,ngq(iq)))
      call b_ematqk(iq,ik, mou_, emat_)
      ! Get 2nd variational eigenstates
      Call getevecsv (vkl0_ptr(1, ik), evecsvt1)
      Call getevecsv (vkl1_ptr(1, ikq), evecsvt0)
      ! Calculate 2nd variational matrix elements
      Do igq=1,ngq(iq)
        call variation_multiplication(evecsvt1,mou_(:,:,igq),evecsvt0,emat(:,:,igq),&
          & bc%n1, bc%n2, bc%il1, bc%il2)
      End Do
      deallocate(mou_)
    end subroutine getmou_sv
    
 !--------------------------------------------------------------------------------
    
    subroutine getmuo_sv(iq, ik, emat, bc)
      use mod_eigenvalue_occupancy, only: nstfv, nstsv
      use modxs, only: bcbs
      use mod_ematptr, only: vkl0_ptr, vkl1_ptr, ikmapikq_ptr
      use modxs, only: ngq, bcbs
      use m_b_ematqk, only: b_ematqk



      implicit none
          
      ! Arguments
      integer, intent(in) :: iq, ik
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:,:)
      ! local variables
      type(bcbs) :: emat_
      complex(8), allocatable :: muo_(:,:,:)
      complex(8) :: evecsvt0(nstsv,nstsv), evecsvt1(nstsv,nstsv)
      integer :: igq, ikq

      ! Find ikq point
      ikq = ikmapikq_ptr(ik, iq)
      ! Set intermediate bands to calculate all combinations
      emat_%n1=nstfv
      emat_%il1=1
      emat_%iu1=nstfv
      emat_%n2=nstfv
      emat_%il2=1
      emat_%iu2=nstfv
      !Calculate all first-variational matrix elements
      allocate(muo_(nstfv,nstfv,ngq(iq)))
      call b_ematqk(iq,ik, muo_, emat_)
      ! Get 2nd variational eigenstates
      Call getevecsv (vkl1_ptr(1, ik), evecsvt1)
      Call getevecsv (vkl0_ptr(1, ikq), evecsvt0)
      ! Calculate 2nd variational matrix elements
      Do igq=1,ngq(iq)
        call variation_multiplication(evecsvt1,muo_(:,:,igq),evecsvt0,emat(:,:,igq),&
          & bc%n1, bc%n2, bc%il1, bc%il2)
      End Do
      deallocate(muo_)
    end subroutine getmuo_sv
    
end module
