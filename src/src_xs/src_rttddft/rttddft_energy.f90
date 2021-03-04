! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020


! Created July 2019 (Ronaldo)
! Modified Jan 2021 (Ronaldo): recoded as a module
! Reference: https://arxiv.org/abs/2102.02630

!> Module to compute the total energy \( E \) for the case of RT-TDDFT calculations
!> \[
!>      E = E_{XC} + E_{Madelung} + E_{eig,core} + E_{ham} - \frac{1}{2}E_{vcl} - E_{vxc}
!> \]
!> where, these components are obtained as follows.
!> 1. XC energy
!> \[
!>    E_{XC} = \int n(\mathbf{r})(e_{XC}) \mathrm{d}\mathbf{r}
!> \]
!> \( e_{XC} \) means the XC energy per particle
!> 2. Madelung energy
!> \[
!>      E_{Madelung}=\frac{1}{2}\sum_{\alpha}z_{\alpha}R_{\alpha},
!> \]
!   where for each atom \( \alpha \) with nuclear charge \( z_{\alpha} \)
!> \[
!>    R_{\alpha} = \lim_{r\rightarrow 0} \left(v^{\rm C}_{\alpha,00}(r)Y_{00}
!>      + \frac{z_{\alpha}}{r} \right)
!> \]
!> with \( v^{\rm C}_{\alpha,00} \) being the \( l=0 \) component of the
!> spherical harmonic expansion of \( v_{\rm C} \) in the muffin-tin region.
!> 3. Contribution of core eigenvalues, \( E_{eig,core} \): sum over all atoms of
!> the eigenvalues obtained for core states
!> 4. Contribution of the hamiltonian: this corresponds to what in groundstate
!> calculations would be the contribution from the valence eigenvalues
!> (for RT-TDDFT, eigenvalues of the hamiltonian do not have the same meaning
!> as in the groundstate)
!> \[
!>    E_{ham} = \sum_{n\mathbf{k}} w_{\mathbf{k}} f_{n\mathbf{k}}
!>    \langle \psi_{n\mathbf{k}} | \hat{h} | \psi_{n\mathbf{k}} \rangle
!> \]
!> where \( w_{\mathbf{k}} \) is the weight of the k-point \( \mathbf{k} \),
!> \( f_{n\mathbf{k}} \) is the occupation of the KS state \( n \) with k-point
!> \mathbf{k}, \( \psi_{n\mathbf{k}} \) is the KS wavefunction.
!> 5. Contribution of the Coulomb potentials:
!> \[
!>    E_{vcl} = \int n(\mathbf{r})v_H(\mathbf{r})\mathrm{d}\mathbf{r}
!> \]
!> where \( v_H \) is the Hartree potential
!> 6. Contribution of the XC correlation potential
!> \[
!>    E_{vxc} = \int n(\mathbf{r}) v_{XC}(\mathbf{r}) \mathrm{d}\mathbf{r}
!> \]

module rttddft_Energy

  use precision, only: dp

  implicit none

  private

  public :: obtain_energy_rttddft, TotalEnergy

  !> This type encapsulates all contributions to the total energy
  type :: TotalEnergy
    !> The following terms account for each different contribution
    !> Exchange \( E_X \)
    real(dp)  :: exchange
    !> Correlation \( E_C \)
    real(dp)  :: correlation
    !> Hartree \( E_{vcl} = \)
    real(dp)  :: Coulomb
    !> XC potential \( \int n(\mathbf{r})v_{XC}(\mathbf{r}) d\mathbf{r} \)
    real(dp)  :: integral_vxc_times_density
    !> Eigenvalues of core states
    real(dp)  :: eigenvalues_core
    !> Madelung
    real(dp)  :: madelung
    !> Hamiltonian. This corresponds to what in groundstate calculations would
    !> be the contribution from the valence eigenvalues (for RT-TDDFT, eigenvalues of the
    !> hamiltonian do not have the same meaning as in the groundstate)
    real(dp)  :: hamiltonian

    !> Finally the total energy itself
    real(dp)  :: total_energy

    contains
      !> total_energy is evaluated from its components
      procedure :: sum_contributions => sum_contr
  end type TotalEnergy

contains
  !> Here, we take into account all the contributions to the total energy to evaluate it
  subroutine sum_contr(this)
    class(TotalEnergy), intent (inout) :: this
    this%total_energy = this%exchange + this%correlation + this%madelung &
      & + this%eigenvalues_core + this%hamiltonian - (0.5_dp)*this%Coulomb &
      & - this%integral_vxc_times_density
  end subroutine sum_contr


  !> Subroutine that calculates the total energy for RT-TDDFT calculations
  !> Adapted from src/energy.f90
  !> @param[in]   first_kpt   first k-point to be considered in the average
  !>                          of the \( E_{ham} \)
  !> @param[in]   last_kpt    last k-point to be considered
  !> @param[in]   ham         Hamiltonian matrix at time \( t \). Dimensions
  !>                          assumed for it: nmatmax, nmatmax, first_kpt:last_kpt
  !> @param[in]   evec        array of coefficients of each KS state in terms
  !>                          of the basis.
  !>                          Dimensions: nmatmax, nstfv, first_kpt:last_kpt
  !> @param[out] rt_tddft_energy
  !>                          Type with the total energy and its components
  subroutine obtain_energy_rttddft(first_kpt,last_kpt,ham,evec,rt_tddft_energy)
    use modinput, only: input
    use mod_kpoint, only: wkpt
    use mod_eigenvalue_occupancy, only: occsv, nstfv
    use mod_eigensystem, only: nmatmax, nmat
    use mod_atoms, only: idxas, natoms, spzn, spnst, nspecies, spcore, spocc
    use mod_potential_and_density, &
      only: rhomt,rhoir,vclmt,vclir,vxcmt,vxcir,exmt,exir,ecmt,ecir, vmad
    use modmpi
    use constants, only: zzero, zone
    use mod_corestate, only: evalcr

    implicit none

    !> first_kpt  first k-point to be considered in the average of the \( E_{ham} \)
    !> last_kpt   last k-point to be considered
    integer,intent(in)              :: first_kpt, last_kpt
    complex(dp), intent(in)         :: ham(:, :, first_kpt:)
    !> evec: array of coefficients of each KS state in terms of the basis.
    !>       Dimensions: nmatmax, nstfv, first_kpt:last_kpt
    complex(dp), intent(in)         :: evec(:, :, first_kpt:)
    !> rt_tddft_energy: Type with the total energy and its components
    type(TotalEnergy), intent(out)  :: rt_tddft_energy


    integer                         :: ik, ist, is, ia, ias, nmatp
    real(dp), allocatable           :: aux(:)
    real(dp)                        :: rfinp
    complex(dp)                     :: acc(nstfv)
    complex(dp),allocatable         :: scratch(:,:),occcmplx(:)
    complex(dp)                     :: zdotc


    allocate(scratch(nmatmax,nmatmax))


    ! contribution of XC and Coulomb potentials, \( v_H \) and \( v_{XC} \), respectively
    rt_tddft_energy%Coulomb = rfinp (1, rhomt, vclmt, rhoir, vclir)
    rt_tddft_energy%integral_vxc_times_density = rfinp (1, rhomt, vxcmt, rhoir, vxcir)

    ! XC energy
    rt_tddft_energy%exchange = rfinp (1, rhomt, exmt, rhoir, exir)
    rt_tddft_energy%correlation = rfinp (1, rhomt, ecmt, rhoir, ecir)

    ! contribution from core eigenvalues
    rt_tddft_energy%eigenvalues_core = 0._dp
    do is = 1, nspecies
      do ia = 1, natoms (is)
        ias = idxas (ia, is)
          do ist = 1, spnst (is)
            if (spcore(ist, is)) rt_tddft_energy%eigenvalues_core = rt_tddft_energy%eigenvalues_core + spocc (ist, is) &
            & * evalcr (ist, ias)
          end do
       end do
    end do

    ! Contribution from the eigenvalues (valence)
    ! They are obtained as the average value of the hamiltonian matrix
    rt_tddft_energy%hamiltonian = 0._dp
    allocate(aux(first_kpt:last_kpt))
    allocate(occcmplx(nstfv))
    aux(:) = 0._dp
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP& PRIVATE(ik, ist, occcmplx, scratch, acc, nmatp), &
!$OMP& SHARED(first_kpt,last_kpt,aux,nmatmax,nstfv,nmat,ham,evec,occsv,wkpt,input)
!$OMP DO
#endif
    do ik = first_kpt, last_kpt
      nmatp = nmat(1,ik)
      ! C := alpha*A*B + beta*C, A hermitian
      ! ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      call ZHEMM( 'L', 'U', nmatp, nstfv, zone, ham(:,:,ik), nmatmax, evec(:,:,ik), &
        & nmatmax, zzero, scratch, nmatmax )

      do ist = 1, nstfv
        ! If the occupation is small, we assume that the current and
        ! all other states with higher "ist" will be unoccupied
        if ( occsv(ist,ik) <= input%groundstate%epsocc ) exit
        acc(ist) = ZDOTC( nmatp, evec(:, ist, ik), 1, scratch(:,ist), 1)
      end do
      occcmplx(:) = occsv(:,ik)
      aux(ik) = dble(wkpt(ik)*ZDOTC( ist, occcmplx(:), 1, acc(:), 1))
    end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
#ifdef MPI
    call MPI_ALLREDUCE(sum(aux), rt_tddft_energy%hamiltonian, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    rt_tddft_energy%hamiltonian = sum(aux)
#endif

    ! Madelung energy
    rt_tddft_energy%madelung = 0._dp
    do is = 1, nspecies
      ! compute the bare nucleus potential at the origin
      do ia = 1, natoms (is)
        ias = idxas (ia, is)
        rt_tddft_energy%madelung = rt_tddft_energy%madelung + 0.5d0 * spzn (is) * vmad(ias)
      end do
    end do

    ! Total energy
    call rt_tddft_energy%sum_contributions()

  end subroutine obtain_energy_rttddft

end module rttddft_Energy
