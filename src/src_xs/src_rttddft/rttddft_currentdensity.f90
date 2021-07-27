! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created: May 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module that deals with the Current Density in RT-TDDFT calculations
module rttddft_CurrentDensity
  implicit none

  private

  public :: UpdateCurrentDensity

contains

  !> Here, we calculate the paramagnetic part of the current density at time 
  !> \( t \). It is calculated as:
  !> \[
  !>    \mathbf{J}(t) = \frac{\mathrm{i}}{\Omega} \sum_{j\mathbf{k}}
  !>      w_{\mathbf{k}}f_{j\mathbf{k}} \left\langle \psi_{j\mathbf{k}}(t) \big|
  !>      \nabla \big|\psi_{j\mathbf{k}}(t)\right\rangle
  !>  \]
  !>  where \( N \) is the number of valence electrons in the unit cell with
  !>  volume \( \Omega \), \( w_{\mathbf{k}} \) is the weight of the considered
  !>  k-point, and \( f_{j\mathbf{k}} \) is the occupation number of the
  !>  corresponding KS state.
  subroutine UpdateCurrentDensity( first_kpt, last_kpt, evec, jpara )
    use precision, only: dp
    use constants, only: zzero, zone
    use modinput, only: input
    use modmpi
    use mod_lattice, only: omega
    use mod_eigenvalue_occupancy, only: occsv, nstfv
    use mod_eigensystem, only: nmat, nmatmax
    use rttddft_GlobalVariables, only: pmat

    implicit none

    !> index of the first `k-point` to be considered in the sum
    integer,intent(in)        :: first_kpt
    !> index of the last `k-point` considered
    integer,intent(in)        :: last_kpt
    !> Basis-expansion coefficients of the KS-wavefunctions at time \( t \).
    !> Dimensions: `nmatmax`, `nstfv`, `first_kpt:last_kpt`
    complex(dp), intent(in)   :: evec(:, :, first_kpt:)
    !> `x`, `y` and `z` components of the parametic current density
    real(8), intent(out)      :: jpara(3)

    integer                   :: ik, ist, j

    real(dp)                  :: weight
    real(dp)                  :: acc(nstfv)
    real(dp)                  :: aux2(3)
    real(dp), allocatable     :: aux(:,:)
    complex(dp), allocatable  :: scratch(:,:)
    ! Blas subroutines
    real(dp)                  :: ddot
    complex(dp)               :: zdotc

    allocate( scratch(nmatmax,nstfv) )
    allocate( aux(3,first_kpt:last_kpt) )

    aux(:,:) = 0._dp

    ! Summation weight of each kpoint
    ! TODO(Ronaldo): consider symmetries in the k-grid
    weight = 1._dp/dble(product(input%groundstate%ngridk))

    ! For the x, y, and z components ...
    do j = 1,3
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP& PRIVATE(ik,ist,scratch,acc), &
!$OMP& SHARED(j,first_kpt,last_kpt,aux,nmatmax,nstfv,pmat), &
!$OMP& SHARED(evec,occsv,nmat,input)
!$OMP DO
#endif
      do ik = first_kpt, last_kpt
        ! C := alpha*A*B + beta*C, A hermitian
        ! ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
        call ZHEMM('L','U',nmat(1,ik),nstfv, zone,pmat(:,:,j,ik),nmatmax, &
          & evec(:,:,ik),nmatmax, zzero,scratch(:,:),nmatmax)
        do ist = 1, nstfv
          ! If the occupation of a certain state is small, we can consider
          ! that all the others above it will have occsv(ist,ik) zero
          if ( occsv(ist, ik)  <= input%groundstate%epsocc ) exit
          acc(ist) = dble(zdotc(nmat(1,ik),evec(:,ist,ik),1,scratch(:,ist),1))
        end do
        aux(j,ik) = -ddot(ist-1,occsv(:,ik),1,acc(:),1)
      end do ! do ik = first_kpt, last_kpt
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
    end do !do j = 1,3

    do j = 1, 3
      aux2(j) = sum( aux(j, first_kpt:last_kpt) )
    end do
#ifdef MPI
    call MPI_ALLREDUCE(aux2, jpara, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    jpara = jpara*weight/omega
#else
    jpara = aux2*weight/omega
#endif
    deallocate( scratch, aux )

  end subroutine UpdateCurrentDensity

end module rttddft_CurrentDensity
