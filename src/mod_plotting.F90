!
!
#include "maxdefinitions.inc"
Module mod_plotting
!-------------------------------------!
!     1D/2D/3D plotting variables     !
!-------------------------------------!
! number of vertices in 1D plot
      Integer :: nvp1d
! total number of points in 1D plot
      Integer :: npp1d
! vertices in lattice coordinates for 1D plot
      Real (8), Allocatable :: vvlp1d (:, :)
! distance to vertices in 1D plot
      Real (8), Allocatable :: dvp1d (:)
! plot vectors in lattice coordinates for 1D plot
      Real (8), Allocatable :: vplp1d (:, :)
! distance to points in 1D plot
      Real (8), Allocatable :: dpp1d (:)
! corner vectors of 2D plot in lattice coordinates
      Real (8) :: vclp2d (3, 3)
! grid sizes of 2D plot
      Integer :: np2d (2)
! corner vectors of 3D plot in lattice coordinates
      Real (8) :: vclp3d (3, 4)
! grid sizes of 3D plot
      Integer :: np3d (3)
! number of states for plotting Fermi surface
!replaced by inputstructureinteger::nstfsp
End Module
!
