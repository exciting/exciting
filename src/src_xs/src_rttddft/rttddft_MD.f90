! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! Created Apr 2021 (Ronaldo)
!> Module for Ehrenfest Dynamics in RT-TDDFT
module rttddft_MD
  use asserts, only: assert
  use constants, only: zone, zzero
  use exciting_mpi, only: mpiinfo, xmpi_allreduce
  use MD, only: MD_input_keys, MD_timing, force, obtain_core_corrections, force_ext, &
    obtain_Hellmann_Feynman_force, obtain_valence_corrections_part1, &
    val_corr_pt2_given_atom_and_kpt => obtain_valence_corrections_part2
  use mod_atoms, only: atposc, idxas, natoms, natmtot, nspecies, &
    spcore, spmass, spocc, spr, spzn
  use mod_corestate, only: rhocr
  use mod_eigensystem, only: nmat, nmatmax
  use mod_eigenvalue_occupancy, only: occsv
  use mod_gkvector, only: ngk, ngkmax, gkc, tpgkc, sfacgk, vgkc
  use mod_gvector, only: ngvec, vgc, sfacg
  use mod_kpoint, only: nkpt, wkpt
  use mod_lattice, only: ainv
  use mod_muffin_tin, only: nrmt
  use mod_potential_and_density, only: vclmt, veffmt, rhomt
  use mod_spin, only: nspnfv
  use modinput, only: input
  use modmpi, only: rank, mpi_env_k, distribute_loop
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_electric_field, only: Electric_Field
  use rttddft_GlobalVariables, only: mathcalH, mathcalB
  use rttddft_timings, only: Print_Timings, timesec_RTTDDFT
  use rttddft_VectorPotential, only: Vector_Potential_Field
  use vector_multiplication, only: dot_multiply

  implicit none 

  private

  public :: force_rttdft, move_ions, allocate_global_arrays, deallocate_global_arrays, &
    evaluate_charge_val

  !> valence charge of each species
  real(dp), allocatable :: charge_val(:)

contains
  !> Allocate all global arrays from this module
  subroutine allocate_global_arrays( n_species )
    !> number of species
    integer(i32), intent(in)           :: n_species

    allocate( charge_val(n_species), source = 0._dp )
  end subroutine

  !> Deallocate all global arrays from this module
  subroutine deallocate_global_arrays( )
    deallocate( charge_val )
  end subroutine

  !> Calculate the valence charge (of each species) 
  !> and store in the global array `charge_val`
  subroutine evaluate_charge_val
    integer(i32) :: is, n_species

    n_species = size( charge_val, 1 )
    ! Sum over occupied states listed as core=false in the species file
    forall( is = 1:n_species ) charge_val(is) = sum( spocc(:, is), mask=(.not.spcore(:, is)) )
  end subroutine

  subroutine force_rttdft( forces, a_tot, e_field, MD_input, evecfv_time, overlap, &
    ham_time, printTimings, t_MD )
    !> Object that packs information about the total forces
    type(force), intent(inout)      :: forces
    !> `x`, `y`, and `z` components of the (total) vector potential
    type(Vector_Potential_Field), intent(in)  :: a_tot
    !> Electric field
    type(Electric_Field), intent(in) :: e_field
    !> Object that contains the inputs keys given in the MD element
    type(MD_input_keys), intent(in) :: MD_input
    !> Basis-expansion coefficients of the KS-WFs at time \(t\)
    complex(dp), intent(in) :: evecfv_time(:, :, :)
    !> Overlap matrix (of basis functions)
    complex(dp), intent(in) :: overlap(:, :, :)
    !> Hamiltonian matrix at current time \(t\)
    complex(dp), intent(in) :: ham_time(:, :, :)
    !> Object that packs information about printing of timings [[Print_Timings]]
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings spent in MD
    class(MD_timing), optional, intent(inout)     :: t_MD

    integer(i32)     :: is, ia, ias, nr, first_kpt, last_kpt
    real(dp)         :: fact, ti
    logical          :: tDetail
    

    ! Check optional (timing) arguments
    tDetail = .false.
    if( present(printTimings) ) tDetail = printTimings%detailed()
    if( tDetail ) then
      call assert( present(t_MD), 't_MD must be present when tDetail is true' )
      call timesec( ti )
    end if

    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)
    fact = dot_multiply(a_tot%components, a_tot%components)/2_dp/c**2
    do is = 1, nspecies
      nr = nrmt(is)
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        forces%EXT(:,ias) = force_ext( charge_val(is), e_field )
        ! Z = -spzn(is): Z is negative in species file
        call obtain_Hellmann_Feynman_force( -spzn(is), spr(1:nr,is), &
          vclmt(:,1:nr,ias), forces%HF(:,ias) )

        if ( MD_input%core_corrections ) &
          call obtain_core_corrections( spr(1:nr,is), rhocr(1:nr,ias), &
            veffmt(:,1:nr,ias), forces%core(:,ias) )

        ! Valence corrections 1: integral of nabla rho_v times (v_KS+A**2/2c**2)
        if ( MD_input%valence_corrections ) &
          call obtain_valence_corrections_part1( spr(1:nr, is), rhocr(1:nr, ias), &
            rhomt(:, 1:nr, ias), veffmt(:,1:nr,ias), fact, forces%val(:,ias) )
      end do ! do ia = 1, natoms (is)
    end do ! do is = 1, nspecies

    if( tDetail ) call timesec_RTTDDFT( ti, t_MD%t_MD_1st )

    ! Valence corrections: second part
    if( MD_input%valence_corrections ) &
      call obtain_valence_corrections_part2( first_kpt, last_kpt, mpi_env_k, &
        evecfv_time, overlap, ham_time, forces%val )
    if( tDetail ) call timesec_RTTDDFT( ti, t_MD%t_MD_2nd )
    ! sum all contributions to total force and store it
    call forces%evaluate_total_force()

    if( tDetail ) call timesec_RTTDDFT( ti, t_MD%t_MD_sumforces )

  end subroutine

  !> Wrapper for calling val_corr_pt2_given_atom_and_kpt
  subroutine obtain_valence_corrections_part2( first_kpt, last_kpt, mpi_env, &
    evecfv_time, overlap, ham_time, forces_val )
    !> index of the first `k-point` to be considered in the sum
    integer(i32),intent(in)        :: first_kpt
    !> index of the last `k-point` considered
    integer(i32),intent(in)        :: last_kpt
    !> MPI environment
    type(mpiinfo), intent(in)      :: mpi_env
    !> Basis-expansion coefficients of the KS-WFs at time \(t\)
    !> (nmatmax, nstfv, first_kpt : last_kpt)
    complex(dp), intent(in) :: evecfv_time(:, :, first_kpt :)
    !> Overlap matrix (of basis functions)
    !> (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), intent(in) :: overlap(:, :, first_kpt :)
    !> Hamiltonian matrix at current time \(t\)
    !> (nmatmax, nmatmax, first_kpt : last_kpt)
    complex(dp), intent(in) :: ham_time(:, :, first_kpt :)
    !> valence corrections to the total force
    real(dp), intent(inout)        :: forces_val(:, :)
    
    integer   :: ik, nmatp, last_occupied, ias
    real(dp)  :: aux(3, natmtot, first_kpt:last_kpt )
    real(dp)  :: sumaux(3, natmtot)
    complex(dp), allocatable :: mathcalS(:, :, :, :)

    call assert( size(forces_val, 1) == 3, 'forces_val must have size = 3 along 1st dim' )
    call assert( size(forces_val, 2) == natmtot, 'forces_val must have size = natmtot along 2nd dim' )
    
    allocate( mathcalS(nmatmax, nmatmax, 3, natmtot) )
    aux = 0._dp
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), &
!$OMP& PRIVATE(ik,ias,nmatp,last_occupied,mathcalS), SHARED(nmat,mathcalH,mathcalB), &
!$OMP& SHARED(natmtot,wkpt,evecfv_time,occsv,aux,first_kpt,last_kpt,overlap,ham_time)
!$OMP DO
#endif
      do ik = first_kpt, last_kpt
        nmatp = nmat(1, ik)
        call obtain_mathcalS( mathcalS, mathcalB(:,:,:,:,ik), overlap(:,:,ik), ham_time(:,:,ik), nmatp )
        last_occupied = first_match( occsv(:, ik) < 1e-4_dp, .true. ) -1
        do ias = 1, natmtot
          call val_corr_pt2_given_atom_and_kpt( mathcalH(1:nmatp,1:nmatp,:,ias,ik), &
            mathcalS(1:nmatp,1:nmatp,:,ias), evecfv_time(1:nmatp,1:last_occupied,ik), &
            occsv(1:last_occupied,ik), aux(:, ias, ik) )
        end do ! do ias = 1, natmtot
        aux(:, :, ik) = -wkpt(ik)*aux(:, :, ik)
      end do ! do ik = 1,nkpt
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      sumaux = sum( aux, dim=3 ) ! sum over kpt
      call xmpi_allreduce( sumaux, mpi_env )
      forces_val = forces_val + sumaux
  end subroutine

  subroutine move_ions( first_kpt, forces, forces_old, dt, atoms_velocities, apwalm, &
    printTimings, t_MD )
    !> The first k point
    integer(i32), intent(in) :: first_kpt
    !> Forces acting on each atom at time \( t \)
    real(dp), intent(in)            :: forces(:, :)
    !> Forces acting on each atom at time \( t - \Delta t \)
    real(dp), intent(in)            :: forces_old(:, :)
    !> Time step for the molecular dynamics
    real(dp), intent(in)            :: dt
    !> Velocities of the nuclei at time \( t \)
    real(dp), intent(inout)         :: atoms_velocities(:,:)
    !> Matching coefficients of the (L)APWs
    !> (ngkmax, apwordmax, lmmaxapw, natmtot, first_kpt : last_kpt)
    complex(dp), intent(inout) :: apwalm(:, :, :, :, first_kpt :)
    !> Object that packs information about printing of timings [[Print_Timings]]
    type(Print_Timings), optional, intent(in) :: printTimings
    !> Object that packs information about timings spent in MD
    class(MD_timing), optional, intent(inout)     :: t_MD
  
    logical                         :: tDetail

    integer                         :: ia, ias, is, ik, ispn, last_kpt
    real(dp)                        :: ti

    call assert( size(forces, 1) == 3, 'forces must have size = 3 along dim = 1' )
    call assert( size(forces_old, 1) == 3, 'forces_old must have size = 3 along dim = 1' )
    call assert( size(atoms_velocities, 1) == 3, 'atoms_velocities must have size = 3 along dim = 1' )
    call assert( size(forces, 2) == natmtot, 'forces must have size = natmtot along dim = 2' )
    call assert( size(forces_old, 2) == natmtot, 'forces_old must have size = natmtot along dim = 2' )
    call assert( size(atoms_velocities, 2) == natmtot, 'atoms_velocities must have size = natmtot along dim = 2' )
  
    ! Check optional (timing) arguments
    tDetail = .False.
    if ( present(printTimings) ) tDetail = printTimings%detailed()
    if( tDetail ) then 
      call assert( present(t_MD), 't_MD must be present when tDetail is true' )
      call timesec( ti )
    end if

    last_kpt = ubound( apwalm, 5 )
  
    do is = 1, nspecies
      do ia = 1, natoms (is)
        ias = idxas (ia, is)
        call update_position_velocity( dt, forces(:,ias)/spmass(is), &
          forces_old(:,ias)/spmass(is), atoms_velocities(:, ias), atposc(:, ia, is) )
        ! obtain the lattice coordinates with the new positions of the ions
        ! call r3mv(ainv,atposc(:,ia,is),input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
      end do
    end do
    if ( tDetail ) call timesec_RTTDDFT( ti, t_MD%t_MD_moveions )
    ! lattice and symmetry set up
    !call findsymcrys ! find the crystal symmetries and shift atomic positions if required
    !call findsymsite ! find the site symmetries
    call checkmt     ! check for overlapping muffin-tins
    call gencfun     ! generate the characteristic function
    call energynn    ! determine the nuclear-nuclear energy
    ! generate structure factors for G and G+k-vectors
    call gensfacgp (ngvec, vgc, ngvec, sfacg)
    do ik = first_kpt, last_kpt
      do ispn = 1, nspnfv
        call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), ngkmax, sfacgk(:, :, ispn, ik))
      end do
    end do
    call gencore          ! generate the core wavefunctions and densities
    call linengy          ! find the new linearization energies
    if (rank == 0) call writelinen
    call genapwfr         ! generate the APW radial functions
    call genlofr          ! generate the local-orbital radial functions
    call olprad
    ! Matching coefficients (apwalm)
    do ik = first_kpt, last_kpt
      call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm(:,:,:,:,ik))
    end do
  
    if ( tDetail ) call timesec_RTTDDFT( ti, t_MD%t_MD_updateBasis )

  end subroutine

  !TODO(Ronaldo): this function can be replaced by findloc( array, condition )
  !However old versions of gfortran (<=8) do not support it
  pure function first_match( array, condition ) 
    logical, intent(in) :: array(:)
    logical, intent(in) :: condition
    integer             :: first_match
    integer             :: i, n
    n = size( array )
    do i = 1, n
      if( array(i) .eqv. condition ) exit
    end do
    first_match = i
  end function

  subroutine update_position_velocity( dt, a_past, a, v, position )
    real(dp), intent(in)    :: dt
    real(dp), intent(in)    :: a(3)
    real(dp), intent(in)    :: a_past(3)
    real(dp), intent(inout) :: v(3)
    real(dp), intent(inout) :: position(3)

    real(dp) :: v_save(3)
    v_save = v
    v = v + 0.5*dt*( a + a_past )
    position = position + 0.5*dt*(v + v_save)
  end subroutine
  
  !> obtain `mathcalS`: an auxiliary matrix used to obtain the valence 
  !> corrections to the total force
  !> \[ \mathcal{S}_{J\mathbf{k}} = 
  !> H_{\mathbf{k}}(S_{\mathbf{k}})^{-1} \mathcal{B}_{J\mathbf{k}} +
  !> \mathcal{B}_{J\mathbf{k}})^{\dagger}(S_{\mathbf{k}})^{-1}H_{\mathbf{k}}
  !> \]
  !> where J labels the atoms; \(H\) and \(S\) are the hamiltonian and overlap 
  !> matrices, respectively
  subroutine obtain_mathcalS( mathcalS, B, S, H, nmatp )
    !> matrix \(\mathcal{S}_{J\mathbf{k}}\) as described before
    complex(dp), intent(out)  :: mathcalS(:, :, :, :)
    !> matrix \(\mathcal{B}_{J\mathbf{k}}\) as described before
    complex(dp), intent(in)   :: B(:, :, :, :)
    !> matrix \(S_{\mathbf{k}}\) as described before
    complex(dp), intent(in)   :: S(:, :)
    !> matrix \(H_{\mathbf{k}}\) as described before
    complex(dp), intent(in)   :: H(:, :)
    !> number of non-zero elements along dim=1 and 2 for all matrices
    integer(i32), intent(in)  :: nmatp

    integer(i32) :: ias, i, j, info, n_atoms, ld
    complex(dp), allocatable  :: aux(:,:), prod(:,:)
  
    call assert( size( mathcalS, 3 ) == 3, 'mathcalS must have size = 3 along dim = 3' )
    do i = 1, 4
      call assert( size( B, i ) == size( mathcalS, i ), 'B and mathcalS must have same size along all dimensions' )
    end do
    do i = 1, 2
      call assert( size( H, i ) == size( S, i ), 'H and S must have same size along all dimensions' )
      call assert( size( H, i ) == size( B, i ), 'H and B must have same size along all dimensions' )
      call assert( size( H, i) >= nmatp, 'H must have size >= nmatp along all dimenstions' )
    end do

    n_atoms = size( mathcalS, 4 )
    ld = size( mathcalS, 1 )
    allocate(aux(nmatp,nmatp), source=S(1:nmatp,1:nmatp))
    allocate(prod(nmatp,nmatp), source=H(1:nmatp,1:nmatp))
    mathcalS = zzero
    ! prod = ((S)**-1)*(H)
    call ZPOSV( 'U', nmatp, nmatp, aux, nmatp, prod, nmatp, info )
    do ias = 1, n_atoms
      ! Loop over x,y,z components
      do j = 1, 3
        ! ZHER2K: C = alpha*B^H*A + cnjg(alpha)*A^H*B + beta*C
        call ZHER2K( 'U', 'C', nmatp, nmatp, zone, prod, nmatp, B(:,:,j,ias), ld, &
          zzero, mathcalS(:,:,j,ias), ld )
        ! hermitize
        do i = 1, nmatp-1
          mathcalS(i+1:nmatp,i,j,ias) = conjg( mathcalS(i,i+1:nmatp,j,ias) )
        end do
      end do
    end do ! do ias = 1, natmtot
  
  end subroutine

end module rttddft_MD