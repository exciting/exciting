!> This module carries global DFPT variables that are used across multiple modules.
module dfpt_variables
  use precision, only: dp
  use modinput
  use mod_kpointset, only: k_set, G_set, Gk_set
  use muffin_tin_basis, only: mt_basis_type
  use block_data_file, only: block_data_file_type

  implicit none
  private

  ! ** reciprocal space
  !> set of electronic \({\bf k}\) vectors
  type(k_set), public :: dfpt_kset

  ! ** interstitial wave function and operator expansion
  !> set of \({\bf G}\) vectors for the expansion of Fourier series
  !> with cutoff `gmaxvr`
  type(G_set), public :: dfpt_Gset
  !> set of \({\bf G}\) vectors with twice the cutoff `2*gmaxvr`
  type(G_set), public :: dfpt_2Gset
  !> set of \({\bf G+k}\) vectors for the (L)APW expansion
  type(Gk_set), public :: dfpt_Gkset

  ! ** muffin tin wave function and operator expansion
  !> maximum angular momentum \(l\) in APW expansion
  integer, public :: dfpt_lmaxapw
  !> maximum angular momentum \(l\) in potential and density expansion
  integer, public :: dfpt_lmaxvr
  !> maximum number of \((l,m)\) pairs in potential and density expansion
  integer, public :: dfpt_lmmaxvr
  !> maximum number of \((l,m)\) pairs for multipoles in Weinert's method
  integer, public :: dfpt_lmaxmp
  !> muffin-tin basis functions
  type(mt_basis_type), public :: mt_basis

  ! ** i/o files
  !> parallel i/o file objects for eigenvalues and eigenvectors on global DFPT \({\bf k}\)-grid
  type(block_data_file_type), public :: fevalk0, feveck0

  public :: dfpt_var_init, dfpt_var_free

  contains

    !> This subroutine initializes global variables that remain constant
    !> during the entire calculation.
    subroutine dfpt_var_init
      use mod_atoms, only: spr, nspecies
      use mod_muffin_tin, only: nrmt
      use mod_APW_LO, only: apwfr, lofr, apword, nlorb, lorbl
      use mod_lattice, only: bvec
      use mod_Gkvector, only: gkmax
      use mod_kpointset, only: generate_k_vectors, generate_G_vectors, generate_Gk_vectors

      ! initialize global variables
      call init0
      call init1

      ! set limits for muffin tin expansion
      dfpt_lmaxapw = input%groundstate%lmaxapw
      dfpt_lmaxvr = input%groundstate%lmaxvr
      dfpt_lmmaxvr = (dfpt_lmaxvr + 1)**2
      dfpt_lmaxmp = max( dfpt_lmaxvr, 30 )

      ! read potential and density from file
      call readstate

      ! set up muffin-tin basis
      call linengy   ! find linearizations energies
      call genapwfr  ! generate APW radial functions
      call genlofr   ! generate LO radial functions
      mt_basis = mt_basis_type( spr(:, 1:nspecies), nrmt(1:nspecies), apwfr, lofr, dfpt_lmaxapw, &
        apword(:, 1:nspecies), nlorb(1:nspecies), lorbl(:, 1:nspecies) )

      ! generate G-vectors for plane wave expansion
      call generate_G_vectors( dfpt_Gset,  bvec, [[0,0,0],[0,0,0]],   input%groundstate%gmaxvr, auto_intgv=.true. )
      call generate_G_vectors( dfpt_2Gset, bvec, [[0,0,0],[0,0,0]], 2*input%groundstate%gmaxvr, auto_intgv=.true. )

      ! generate k-point set
      call generate_k_vectors( dfpt_kset, bvec, &
             input%groundstate%ngridk, &
             input%groundstate%vkloff, &
             input%groundstate%reducek, &
             uselibzint=.false. )
           
      ! generate G+k vectors
      call generate_Gk_vectors( dfpt_Gkset, dfpt_kset, dfpt_Gset, gkmax )
    end subroutine dfpt_var_init

    !> This subroutine frees memory from global variables.
    subroutine dfpt_var_free
      use mod_kpointset, only: delete_k_vectors, delete_G_vectors, delete_Gk_vectors
      call delete_k_vectors( dfpt_kset )
      call delete_G_vectors( dfpt_Gset )
      call delete_G_vectors( dfpt_2Gset )
      call delete_Gk_vectors( dfpt_Gkset )
      call mt_basis%destroy
    end subroutine dfpt_var_free

end module dfpt_variables
