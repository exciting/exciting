! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Reference: https://arxiv.org/abs/2102.02630

!> Module to obtain the number of excitations
!> The number of excited electrons after the interaction with a laser pulse
!> In RT-TDDFT, the occupation number \( f_{j\mathbf{k}} \) of a KS state is
!> kept fixed to its initial value. As the wavefunctions evolve, they are not
!> any longer eigenstates of \( \hat{H}(t) \). It is possible to describe
!> the number of excitations by projecting \( | \psi_{i\mathbf{k}}(t)\rangle \)
!> onto the reference ground state at \( t=0 \).
!> For a given k-point, we define the number of electrons that have
!> been excited to an unoccupied KS state, labeled  \( j \), as
!> \[
!> 	m_{j\mathbf{k}}(t)= \sum_{i} f_{i\mathbf{k}}| \langle \psi_{j\mathbf{k}}(0)
!>	         | \psi_{i\mathbf{k}}(t)\rangle |^2.
!> Similarly, the number of holes created in an occupied KS \( j \) state can
!> specified as
!> 	\[
!> 	m_{j'\mathbf{k}}(t)= f_{j'\mathbf{k}} - \sum_{i}
!> 	f_{i\mathbf{k}}	| \langle \psi_{j'\mathbf{k}}(0)| \psi_{i\mathbf{k}}(t)\rangle |^2.
!> 	\]
!> Thus, the total number of excited electrons in a unit cell can be
!> obtained by considering all the unoccupied states
!> \[
!> 	N_{exc}(t)=
!> 	\sum_{j\mathbf{k}}^{j\, unocc}
!> 	w_\mathbf{k} m_{j\mathbf{k}}(t) = \sum_{j'\mathbf{k}}^{j'\, occ}
!> 	w_\mathbf{k} m_{j'\mathbf{k}}(t) .
!> 	\]
! HISTORY
! Created July 2019 (Ronaldo)
module rttddft_NumberExcitations
  use precision, only: dp

  implicit none

  private

  public :: Obtain_number_excitations

contains

  !> With this subroutine, we obtain the number of excitations
  !> @param[in]   first_kpt   first k-point to be considered in the average
  !> @param[in]   last_kpt    last k-point
  !> @param[in]   evec_init   coefficients of the KS-wavefunctions at t=0
  !>                          Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  !> @param[in]   evec_time   coefficients of the KS-wavefunctions at current time
  !>                          Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  !> @param[in]   overlap     overlap matrix
  !>                          Dimensions: nmatmax, nmatmax, first_kpt:last_kpt
  !> @param[out]  nex         number of excited electrons
  !> @param[out]  ngs         number of electrons on the groundstate state
  !> @param[out]  nt          total number of electrons, i.e., ngs + nex
subroutine Obtain_number_excitations( first_kpt, last_kpt, evec_init, evec_time, overlap, &
    & nex, ngs, nt )
  use mod_kpoint, only: wkpt
  use modinput, only: input
  use mod_eigenvalue_occupancy, only: occsv, nstfv
  use mod_eigensystem, only: nmatmax, nmat
  use modmpi
  use constants, only: zzero, zone
  implicit none

  !> first_kpt   first k-point to be considered in the average
  !> last_kpt    last k-point
  integer,intent(in)        :: first_kpt, last_kpt
  !> evec_init   coefficients of the KS-wavefunctions at t=0
  !>             Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  complex(dp), intent(in)   :: evec_init(:, :,first_kpt:)
  !> evec_time   coefficients of the KS-wavefunctions at current time
  !>             Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  complex(dp), intent(in)   :: evec_time(:, :, first_kpt:)
  !> overlap     overlap matrix, Dimensions: nmatmax, nmatmax, first_kpt:last_kpt
  complex(dp), intent(in)   :: overlap(:, :, first_kpt:)
  !> nex         number of excited electrons
  real(dp), intent(out)     :: nex
  !> ngs         number of electrons on the groundstate state
  real(dp), intent(out)     :: ngs
  !> nt          total number of electrons, i.e., ngs + nex
  real(dp), intent(out)     :: nt

  integer                   :: ik, ist, jst, nmatp
  real(dp)                  :: aux
  complex(dp), allocatable  :: scratch(:,:),proj(:,:)


  allocate( scratch(nmatmax, nstfv) )
  allocate( proj(nstfv, nstfv) )

  nex = 0._dp
  ngs = 0._dp
  nt  = 0._dp
  do ik = first_kpt, last_kpt
    nmatp = nmat(1,ik)
    ! C := alpha*A*B + beta*C, A hermitian
    ! ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    ! scratch = overlap*evect_time
    call ZHEMM( 'L', 'U', nmatp, nstfv, zone, overlap(:,:,ik), nmatmax, &
      & evec_time(:,:,ik), nmatmax, zzero, scratch, nmatmax )
    ! Matrix multiplication: proj = (evec_init**H)*scratch
    call ZGEMM( 'C', 'N', nstfv, nstfv, nmatp, zone, evec_init(:,:,ik), nmatmax, &
      & scratch(:,:), nmatmax, zzero, proj(:,:), nstfv )
    do ist = 1, nstfv
      do jst = 1, nstfv
        ! If the occupation is small, we assume that the current and
        ! all other states with higher "jst" will be unoccupied
        if ( occsv(jst, ik)  <= input%groundstate%epsocc ) exit
        ! If the state "jst" is occupied, then we follow
        aux = wkpt(ik) * occsv(jst, ik) * ( abs( proj(ist, jst) )**2 )
        nt = nt + aux
        if ( occsv(ist, ik)  <= input%groundstate%epsocc ) then
          nex = nex + aux
        else
          ngs = ngs + aux
        endif
      end do
    end do
  end do
#ifdef MPI
  call MPI_ALLREDUCE(nex, aux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  nex = aux
  call MPI_ALLREDUCE(ngs, aux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  ngs = aux
  call MPI_ALLREDUCE(nt, aux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  nt = aux
#endif

  deallocate(scratch)
  deallocate(proj)
end subroutine Obtain_number_excitations
end module rttddft_NumberExcitations
