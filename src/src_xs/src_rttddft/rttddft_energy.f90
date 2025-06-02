! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020


! Created July 2019 (Ronaldo)
! Modified Jan 2021 (Ronaldo): recoded as a module
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module to compute the total energy \( E \) for the case of RT-TDDFT calculations
!> \[
!>      E = E_{XC} + E_{Madelung} + E_{eig,core} + E_{ham} - \frac{1}{2}E_{vcl} 
!>          - E_{vxc}
!> \]
!> where, these components are obtained as follows.
!> <ol>
!> <li> XC energy
!> \[
!>    E_{XC} = \int n(\mathbf{r})e_{XC}(\mathbf{r}) \mathrm{d}\mathbf{r}
!> \]
!> \( e_{XC}(\mathbf{r}) \) means the XC energy per particle
!> </li>
!> <li> Madelung energy
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
!> </li>
!> <li> Contribution of core eigenvalues, \( E_{eig,core} \): sum over all atoms of
!> the eigenvalues obtained for core states
!> </li>
!> <li> Contribution of the hamiltonian: this corresponds to what in groundstate
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
!> </li>
!> <li> Contribution of the Coulomb potentials:
!> \[
!>    E_{vcl} = \int n(\mathbf{r})v_H(\mathbf{r})\mathrm{d}\mathbf{r}
!> \]
!> where \( v_H \) is the Hartree potential
!> </li>
!> <li> Contribution of the XC correlation potential
!> \[
!>    E_{vxc} = \int n(\mathbf{r}) v_{XC}(\mathbf{r}) \mathrm{d}\mathbf{r}
!> \]
!> </li>
!> </ol>
module rttddft_Energy
  use asserts, only: assert
  use constants, only: real_zero
  use exciting_mpi, only: mpiinfo, xmpi_allreduce
  use hermitian_matrix_multiplication, only: hermitian_matrix_multiply
  use modinput, only: input
  use mod_atoms, only: idxas, natoms, spzn, spnst, nspecies, spcore, spocc
  use mod_corestate, only: evalcr
  use mod_eigensystem, only: nmatmax, nmat
  use mod_potential_and_density, only: rhomt, rhoir, vclmt, vclir, vxcmt, &
    vxcir, exmt, exir, ecmt, ecir, vmad  
  use precision, only: dp, i32
  use rttddft_Wavefunction, only: wavefunction_set
  use vector_multiplication, only: dot_multiply

  implicit none

  private

  public :: obtain_energy_rttddft, TotalEnergy

  !> This type encapsulates all contributions to the total energy.  
  !> The following terms account for each different contribution
  type :: TotalEnergy
    !> Exchange \( E_X \)
    real(dp) :: exchange
    !> Correlation \( E_C \)
    real(dp) :: correlation
    !> Hartree \( E_{vcl} \)
    real(dp) :: Coulomb
    !> XC potential \( \int n(\mathbf{r})v_{XC}(\mathbf{r}) d\mathbf{r} \)
    real(dp) :: integral_vxc_times_density
    !> Eigenvalues of core states
    real(dp) :: eigenvalues_core
    !> Madelung
    real(dp) :: madelung
    !> Hamiltonian. This corresponds to what in groundstate calculations would
    !> be the contribution from the valence eigenvalues (for RT-TDDFT, 
    !> eigenvalues of the hamiltonian do not have the same meaning as in the 
    !> groundstate)
    real(dp) :: hamiltonian
    !> The total energy itself
    real(dp) :: total_energy

    contains
      !> total_energy is evaluated from its components
      procedure :: sum_contributions => sum_contr
  end type TotalEnergy

contains
  !> Here, we take into account all the contributions to the total energy to 
  !> evaluate it as
  !> \[
  !>      E = E_{XC} + E_{Madelung} + E_{eig,core} + E_{ham} - \frac{1}{2}E_{vcl} 
  !>          - E_{vxc}
  !> \]
  subroutine sum_contr(this)
    class(TotalEnergy), intent (inout) :: this
    this%total_energy = this%exchange + this%correlation + this%madelung &
      & + this%eigenvalues_core + this%hamiltonian - (0.5_dp)*this%Coulomb &
      & - this%integral_vxc_times_density
  end subroutine sum_contr


  !> Subroutine that calculates the total energy for RT-TDDFT calculations
  !> Adapted from `src/energy.f90`
  subroutine obtain_energy_rttddft(first_kpt, ham, psi, occupations, initial_ks_energies, &
      mpi_env, kpt_weights, rt_tddft_energy )

    implicit none

    !> index of the first `k-point` to be considered in the sum appearing in 
    !> \( E_{ham} \)
    integer,intent(in) :: first_kpt
    !> Hamiltonian matrix at time \( t \). 
    !> Dimensions: `nmatmax`, `nmatmax`, `first_kpt:last_kpt`
    complex(dp), intent(in) :: ham(:, :, :)
    !> Basis-expansion coefficients of the KS-wavefunctions at time \( t \)
    class(wavefunction_set), intent(in) :: psi
    !> Initial occupations array
    real(dp), intent(in) :: occupations(:, :)
    !> Initial KS energies array
    real(dp), intent(in) :: initial_ks_energies(:, :)
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env
    !> k points weights array
    real(dp), intent(in) :: kpt_weights(:)
    !> Type with the total energy and its components
    type(TotalEnergy), intent(out) :: rt_tddft_energy


    integer(i32) :: ik, ist, is, ia, ias, nmatp, real_kpt, n_kpt, first_active, n_states, n_basis, n_frozen
    real(dp), allocatable :: aux(:)
    real(dp) :: rfinp
    complex(dp), allocatable :: acc(:), scratch(:, :), occcmplx(:)

    n_kpt = psi%n_kpts()
    n_states = psi%n_active()
    n_basis = psi%n_basis()
    first_active = psi%first_active()
    n_frozen = psi%n_frozen()

    call assert( n_kpt == size( occupations, 2 ), 'psi and occupations have different nkpts' )

    allocate( scratch(n_basis, n_states) )
    allocate( acc(n_states) )

    ! contribution of XC and Coulomb potentials, \( v_H \) and \( v_{XC} \), respectively
    rt_tddft_energy%Coulomb = rfinp (1, rhomt, vclmt, rhoir, vclir)
    rt_tddft_energy%integral_vxc_times_density = rfinp (1, rhomt, vxcmt, rhoir, vxcir)

    ! XC energy
    rt_tddft_energy%exchange = rfinp (1, rhomt, exmt, rhoir, exir)
    rt_tddft_energy%correlation = rfinp (1, rhomt, ecmt, rhoir, ecir)

    ! contribution from core eigenvalues
    rt_tddft_energy%eigenvalues_core = 0._dp
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
          do ist = 1, spnst(is)
            if ( spcore(ist, is) ) rt_tddft_energy%eigenvalues_core = rt_tddft_energy%eigenvalues_core + spocc(ist, is) &
            & * evalcr (ist, ias)
          end do
       end do
    end do

    ! Contribution from the eigenvalues (valence)
    ! They are obtained as the average value of the hamiltonian matrix
    rt_tddft_energy%hamiltonian = 0._dp
    allocate( aux(n_kpt), source = real_zero )
    allocate( occcmplx(n_states) )
    !$omp parallel default(none), &
    !$omp private(ik, ist, occcmplx, scratch, acc, nmatp, real_kpt), &
    !$omp shared(first_kpt, aux, first_active, nmat, ham, psi, occupations, &
    !$omp kpt_weights, input, n_kpt, n_states, initial_ks_energies, n_frozen)
    !$omp do
    do ik = 1, n_kpt
      real_kpt = ik + first_kpt - 1
      if ( psi%expanded_in_lapwlo() ) then
        nmatp = nmat(1, real_kpt)
      else
        nmatp = psi%n_basis()
      end if
      call hermitian_matrix_multiply( ham(:, :, ik), psi%active(:, :, ik), scratch, side='L', uplo='U' )

      do ist = 1, n_states
        ! If the occupation is small, we assume that the current and
        ! all other states with higher "ist" will be unoccupied
        if ( occupations(n_frozen + ist, ik) <= input%groundstate%epsocc ) exit
        acc(ist) = dot_multiply( psi%active(1:nmatp, ist, ik), scratch(1:nmatp, ist), conjg_a=.true. )
      end do
      occcmplx = occupations(first_active : n_frozen + n_states, ik)
      aux(ik) = kpt_weights(ik) * real( dot_multiply( occcmplx(1:ist - 1), acc(1:ist - 1) ), dp )
      if ( psi%has_frozen() ) aux(ik) = aux(ik) + kpt_weights(ik) * &
        dot_multiply( occupations(1 : n_frozen, ik), initial_ks_energies(1 : n_frozen, ik) )
    end do
    !$omp end parallel
    rt_tddft_energy%hamiltonian = sum( aux )

    call xmpi_allreduce( rt_tddft_energy%hamiltonian, mpi_env )

    ! Madelung energy
    rt_tddft_energy%madelung = 0._dp
    do is = 1, nspecies
      ! compute the bare nucleus potential at the origin
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        rt_tddft_energy%madelung = rt_tddft_energy%madelung + 0.5_dp * spzn(is) * vmad(ias)
      end do
    end do

    ! Total energy
    call rt_tddft_energy%sum_contributions()

  end subroutine obtain_energy_rttddft

end module rttddft_Energy
