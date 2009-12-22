!
!
#include "maxdefinitions.inc"
Module mod_Gkvector
!----------------------------------!
!     G+k-vector set variables     !
!----------------------------------!
! smallest muffin-tin radius times gkmax
!replaced by inputstructurereal(8)::rgkmax
! maximum |G+k| cut-off for APW functions
      Real (8) :: gkmax
! number of G+k-vectors for augmented plane waves
      Integer, Allocatable :: ngk (:, :)
! maximum number of G+k-vectors over all k-points
      Integer :: ngkmax
! index from G+k-vectors to G-vectors
      Integer, Allocatable :: igkig (:, :, :)
! G+k-vectors in lattice coordinates
      Real (8), Allocatable :: vgkl (:, :, :, :)
! G+k-vectors in Cartesian coordinates
      Real (8), Allocatable :: vgkc (:, :, :, :)
! length of G+k-vectors
      Real (8), Allocatable :: gkc (:, :, :)
! (theta, phi) coordinates of G+k-vectors
      Real (8), Allocatable :: tpgkc (:, :, :, :)
! structure factor for the G+k-vectors
      Complex (8), Allocatable :: sfacgk (:, :, :, :)
End Module
!
