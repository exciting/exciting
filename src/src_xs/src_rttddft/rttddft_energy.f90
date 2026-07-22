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
  use rttddft_Hamiltonian, only: hamiltonian_set
  use rttddft_Wavefunction, only: wavefunction_set
  use vector_multiplication, only: dot_multiply

  implicit none

  private

  public :: Total_Energy

  !> This type encapsulates all contributions to the total energy.  
  !> The following terms account for each different contribution
  type :: Total_Energy
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

    contains
      !> total_energy is evaluated from its components
      procedure, public :: total_energy => sum_contr
      !> calculate each contribution
      procedure, public :: calculate => obtain_energy_rttddft
  end type Total_Energy

contains
  !> Here, we take into account all the contributions to the total energy to 
  !> evaluate it as
  !> \[
  !>      E = E_{XC} + E_{Madelung} + E_{eig,core} + E_{ham} - \frac{1}{2}E_{vcl} 
  !>          - E_{vxc}
  !> \]
  real(dp) pure function sum_contr(this)
    class(Total_Energy), intent (in) :: this
    sum_contr = this%exchange + this%correlation + this%madelung &
      & + this%eigenvalues_core + this%hamiltonian - (0.5_dp)*this%Coulomb &
      & - this%integral_vxc_times_density
end function sum_contr


  !> Subroutine that calculates all contributions to the total energy for RT-TDDFT calculations. 
  !> Adapted from `src/energy.f90`
  subroutine obtain_energy_rttddft(this, H, psi, mpi_env )
    class(Total_Energy), intent(inout) :: this
    !> Hamiltonian matrix at time \( t \). 
    !> Object that packs information about the Hamiltonian
    class(hamiltonian_set), intent(in) :: H
    !> Basis-expansion coefficients of the KS-wavefunctions at time \( t \)
    class(wavefunction_set), intent(in) :: psi
    !> MPI environment
    type(mpiinfo), intent(in) :: mpi_env

    integer(i32) :: ik, ist, is, ia, ias, nmatp, first_active, n_states, n_basis, n_frozen
    real(dp), allocatable :: aux(:)
    real(dp) :: rfinp
    complex(dp), allocatable :: acc(:), scratch(:, :), occcmplx(:)

    n_states = psi%n_active()
    n_basis = psi%n_basis()
    first_active = psi%first_active()
    n_frozen = psi%n_frozen()

    allocate( scratch(n_basis, n_states) )
    allocate( acc(n_states) )

    ! contribution of XC and Coulomb potentials, \( v_H \) and \( v_{XC} \), respectively
    this%Coulomb = rfinp (1, rhomt, vclmt, rhoir, vclir)
    this%integral_vxc_times_density = rfinp (1, rhomt, vxcmt, rhoir, vxcir)

    ! XC energy
    this%exchange = rfinp (1, rhomt, exmt, rhoir, exir)
    this%correlation = rfinp (1, rhomt, ecmt, rhoir, ecir)

    ! contribution from core eigenvalues
    this%eigenvalues_core = 0._dp
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
          do ist = 1, spnst(is)
            if ( spcore(ist, is) ) this%eigenvalues_core = &
              this%eigenvalues_core + spocc(ist, is)*evalcr(ist, ias)
          end do
       end do
    end do

    ! Contribution from the eigenvalues (valence): obtained as the expected value 
    ! of the hamiltonian matrix
    this%hamiltonian = 0._dp
    allocate( aux(psi%first_kpt():psi%last_kpt()), source = real_zero )
    allocate( occcmplx(n_states) )
    !$omp parallel default(none), &
    !$omp private(ik, ist, occcmplx, scratch, acc, nmatp), &
    !$omp shared(aux, first_active, nmat, H, psi,input, n_states, n_frozen)
    !$omp do
    do ik = psi%first_kpt(), psi%last_kpt()
      if ( psi%expanded_in_lapwlo() ) then
        nmatp = nmat(1, ik)
      else
        nmatp = psi%n_basis()
      end if
      call hermitian_matrix_multiply( H%H_t%array(:, :, ik), psi%active(:, :, ik), scratch, side='L', uplo='U' )

      do ist = 1, n_states
        ! If the occupation is small, we assume that the current and
        ! all other states with higher "ist" will be unoccupied
        if ( psi%occupations(n_frozen + ist, ik) <= psi%eps_occ ) exit
        acc(ist) = dot_multiply( psi%active(1:nmatp, ist, ik), scratch(1:nmatp, ist), conjg_a=.true. )
      end do
      occcmplx = psi%occupations(first_active : n_frozen + n_states, ik)
      aux(ik) = psi%kset%wkpt(ik) * real( dot_multiply( occcmplx(1:ist - 1), acc(1:ist - 1) ), dp )
      if ( psi%has_frozen() ) aux(ik) = aux(ik) + psi%kset%wkpt(ik) * &
        dot_multiply( psi%occupations(1 : n_frozen, ik), H%initial_eigenvalues(1 : n_frozen, ik) )
    end do
    !$omp end parallel
    this%hamiltonian = sum( aux )

    call xmpi_allreduce( this%hamiltonian, mpi_env )

    ! Madelung energy
    this%madelung = 0._dp
    do is = 1, nspecies
      ! compute the bare nucleus potential at the origin
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        this%madelung = this%madelung + 0.5_dp * spzn(is) * vmad(ias)
      end do
    end do
  end subroutine obtain_energy_rttddft

end module rttddft_Energy
