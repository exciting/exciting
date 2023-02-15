!> This module carries global DFPT phonon variables that are used across multiple modules.
module phonons_variables
  use dfpt_variables
  use phonons_symmetry, only: irrep_basis
  use phonons_parallelization, only: ph_part

  use precision, only: dp
  use modmpi
  use asserts, only: assert
  use mod_kpointset, only: k_set, Gk_set, G_set

  implicit none
  private

  ! ** general
  !> \({\bf k}\)- and \({\bf q}\)-point set are commensurate
  logical, public :: ph_commensurate = .false.
  !> use canonical displacement patterns
  logical, public :: ph_canonical = .false.
  !> number of phonon modes
  integer, public :: ph_nmodes

  ! ** reciprocal space
  !> the phonon \({\bf q}\)-point set
  type(k_set), public :: ph_qset
  !> the \({\bf q}\)-dependent electronic \({\bf k}\) point set \(\text{IBZ}({\bf q})\)
  type(k_set), public :: ph_kset
  !> the set of \({\bf k+q}\)-points
  type(k_set), public :: ph_kqset

  ! ** interstitial wave function and operator expansion
  !> set of \({\bf G+q}\) vectors corresponding to `dfpt_Gset`
  type(G_set), public :: ph_Gqset
  !> set of \({\bf G+q}\) vectors corresponding to `dfpt_2Gset`
  type(G_set), public :: ph_2Gqset
  !> set of \({\bf G+k}\) vectors
  type(Gk_set), public :: ph_Gkset
  !> set of \({\bf G+k+q}\) vectors
  type(Gk_set), public :: ph_Gkqset

  ! ** symmetry adapted basis of displacement patterns
  !> symmetry adapted phonon expansion basis for each \({\bf q}\) point
  type(irrep_basis), allocatable, public :: ph_irrep_basis(:)

  ! ** auxiliary quantities
  !> spherical Bessel functions \(j_l(|{\bf G+q}|R_\alpha)\)
  real(dp), allocatable, public :: ph_jlgqr(:,:,:)
  !> spherical harmonics \(Y_{lm}(\widehat{{\bf G+q}})\)
  complex(dp), allocatable, public :: ph_ylmgq(:,:)
  !> structure factors \({\rm e}^{{\rm i} ({\bf G+q})\cdot{\bf \tau}_\alpha}\)
  complex(dp), allocatable, public :: ph_sfacgq(:,:)

  ! ** parallelization variables
  !> number of processes
  integer, public :: ph_numprocs
  !> independent parts of phonon calculation (to be calculated)
  type(ph_part), allocatable, public :: ph_parts(:)
  !> independent parts of phonon calculation (all)
  type(ph_part), allocatable, public :: ph_parts_all(:)
  !> flag if \(({\bf q},I)\) pair has already been precomputed
  logical, allocatable, public :: ph_parts_done(:,:)
  !> parts this (global) process works on in their respective order
  integer, allocatable, public :: ph_parts_per_rank(:)
  !> schedule telling which process works on which part in which order
  integer, allocatable, public :: ph_schedule(:,:)

  public :: ph_var_init, ph_var_free
  public :: ph_var_init_q, ph_var_free_q

  contains

    !> This subroutine initializes global variables that remain constant
    !> during the entire phonon calculation.
    subroutine ph_var_init
      use phonons_symmetry, only: ph_sym_find_irreps
      use mod_kpointset, only: generate_k_vectors
      use mod_atoms, only: natmtot
      use modinput

      integer :: iq
      real(dp) :: v1(3), v2(3)

      ! check if canonical displacement patterns should be used
      ph_canonical = input%phonons%canonical

      ! generate q-point set
      call generate_k_vectors( ph_qset, dfpt_kset%bvec, &
             input%phonons%ngridq, &
             [0.0_dp, 0.0_dp, 0.0_dp], &
             input%phonons%reduceq, &
             uselibzint=.false. )

      ! check if k- and q-point set are commensurate
      v1 = dble( dfpt_kset%ngridk ) / dble( ph_qset%ngridk )
      v2 = ph_qset%vkloff * v1
      ph_commensurate = (norm2( v1(:) - nint( v1(:) ) ) < input%structure%epslat ) .and. &
                        (norm2( v2(:) - nint( v2(:) ) ) < input%structure%epslat )
      ! TODO for the moment we require the grids to be commensurate
      ! in the future, we want to scrap that necessity
      call terminate_if_false( ph_commensurate, '(ph_var_init) &
        k- and q-grid are not commensuarate.' )

     ! find basis of irreducible representations (symmetry adapted displacement patterns)
     ph_nmodes = 3 * natmtot
     allocate( ph_irrep_basis(0:ph_qset%nkpt) )
     do iq = 1, ph_qset%nkpt
       call ph_sym_find_irreps( ph_qset%vkl(:, iq), ph_irrep_basis(iq), canonical=ph_canonical )
     end do
     ! the 0th basis contains the symmetries of the unpertubed system
     call ph_sym_find_irreps( [0.0_dp, 0.0_dp, 0.0_dp], ph_irrep_basis(0), canonical=.false., unperturbed=.true. )

     ! set number of processes
     ph_numprocs = mpiglobal%procs

     ! allocate list of done parts
     allocate( ph_parts_done(ph_qset%nkpt, 3*natmtot), source=.false. )
    end subroutine ph_var_init

    !> This subroutine frees memory from global variables.
    subroutine ph_var_free
      use mod_kpointset, only: delete_k_vectors
      integer :: iq
      if( allocated( ph_parts_done ) ) deallocate( ph_parts_done )
      if( allocated( ph_parts ) ) deallocate( ph_parts )
      if( allocated( ph_parts_all ) ) deallocate( ph_parts_all )
      if( allocated( ph_irrep_basis ) ) then
        do iq = 0, ph_qset%nkpt
          call ph_irrep_basis(iq)%free
        end do
        deallocate( ph_irrep_basis )
      end if
      call delete_k_vectors( ph_qset )
    end subroutine ph_var_free

    !> This subroutine initializes global variables that are \({\bf q}\)-dependent.
    subroutine ph_var_init_q( iq )
      use mod_kpointset, only: generate_k_vectors, generate_Gk_vectors, generate_G_vectors
      use mod_symmetry, only: maxsymcrys, lsplsymc, nsymcrys
      use mod_muffin_tin, only: rmt
      use mod_atoms, only: natmtot, nspecies, natoms, idxas, atposc
      use modinput
      !> index of point in \({\bf q}\) point set
      !> (if `iq=0` the Gamma point is assumed)
      integer, intent(in) :: iq

      integer :: is, ia, ias, ikq, ig, lmmax
      integer :: lsplsymc_(maxsymcrys), nsymcrys_
      real(dp) :: vql(3), vqc(3), v(3), gc, t1

      ! copy q-vector
      call assert( iq >= 0 .and. iq <= ph_qset%nkpt, '(ph_var_init_q): &
        Invalid q-point index.' )
      if( iq == 0 ) then
        vql = 0.0_dp
        vqc = 0.0_dp
      else
        vql = ph_qset%vkl(:, iq)
        vqc = ph_qset%vkc(:, iq)
      end if

      ! store global variables and overwrite them
      lsplsymc_ = lsplsymc
      nsymcrys_ = nsymcrys
      do is = 1, ph_irrep_basis(iq)%nsym
        lsplsymc(is) = lsplsymc_(ph_irrep_basis(iq)%isym(is))
      end do
      nsymcrys = ph_irrep_basis(iq)%nsym

      ! build k-vectors (IBZ of q)
      call generate_k_vectors( ph_kset, dfpt_kset%bvec, &
             dfpt_kset%ngridk, &
             dfpt_kset%vkloff, &
             dfpt_kset%isreduced, &
             uselibzint=.false. )

      ! generate G+k vectors
      call generate_Gk_vectors( ph_Gkset, ph_kset, dfpt_Gset, dfpt_Gkset%gkmax )

      ! restore global variables
      lsplsymc = lsplsymc_
      nsymcrys = nsymcrys_

      ! build k+q-vectors
      ph_kqset = ph_kset
      do ikq = 1, ph_kqset%nkpt
        ph_kqset%vkl(:, ikq) = ph_kqset%vkl(:, ikq) + vql
        ph_kqset%vkc(:, ikq) = ph_kqset%vkc(:, ikq) + vqc
      end do

      ! build G+q and G+k+q-vectors
      call generate_G_vectors( ph_Gqset, dfpt_Gset%bvec, dfpt_Gset%intgv, dfpt_Gset%gmaxvr, vpl=vql )
      call generate_G_vectors( ph_2Gqset, dfpt_2Gset%bvec, dfpt_2Gset%intgv, dfpt_2Gset%gmaxvr, vpl=vql )
      call generate_Gk_vectors( ph_Gkqset, ph_kqset, dfpt_Gset, dfpt_Gkset%gkmax )
      
      ! build auxiliary variables
      lmmax = (dfpt_lmaxmp + 1)**2
      allocate( ph_jlgqr(0:(dfpt_lmaxmp+input%groundstate%npsden+1), ph_Gqset%ngvec, nspecies) )
      allocate( ph_ylmgq(lmmax, ph_Gqset%ngvec) )
      allocate( ph_sfacgq(ph_Gqset%ngvec, natmtot) )
      
!$omp parallel default( shared ) private( ig, v, gc, is, ia, ias, t1 )
!$omp do
      do ig = 1, ph_Gqset%ngvec
        v = ph_Gqset%vgc(:, ig)
        gc = ph_Gqset%gc(ig)
        call ylm( v, dfpt_lmaxmp, ph_ylmgq(:, ig) )
        do is = 1, nspecies
          call sbessel( dfpt_lmaxmp+input%groundstate%npsden+1, gc*rmt(is), ph_jlgqr(:, ig, is) )
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            t1 = dot_product( v, atposc(:, ia, is) )
            ph_sfacgq(ig, ias) = cmplx( cos( t1 ), sin( t1 ), dp )
          end do
        end do
      end do
!$omp end do
!$omp end parallel
    end subroutine ph_var_init_q

    !> This subroutine frees memory from global \({\bf q}\)-dependent variables.
    subroutine ph_var_free_q
      use mod_kpointset, only: delete_k_vectors, delete_G_vectors, delete_Gk_vectors
      call delete_k_vectors( ph_kset )
      call delete_k_vectors( ph_kqset )
      call delete_Gk_vectors( ph_Gkset )
      call delete_Gk_vectors( ph_Gkqset )
      call delete_G_vectors( ph_Gqset )
      call delete_G_vectors( ph_2Gqset )
      if( allocated( ph_jlgqr ) ) deallocate( ph_jlgqr )
      if( allocated( ph_ylmgq ) ) deallocate( ph_ylmgq )
      if( allocated( ph_sfacgq ) ) deallocate( ph_sfacgq )
    end subroutine ph_var_free_q

end module phonons_variables
