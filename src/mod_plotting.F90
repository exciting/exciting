

#include "maxdefinitions.inc"
module mod_plotting
!-------------------------------------!
!     1D/2D/3D plotting variables     !
!-------------------------------------!
! number of vertices in 1D plot
integer::nvp1d
! total number of points in 1D plot
integer::npp1d
! vertices in lattice coordinates for 1D plot
real(8), allocatable :: vvlp1d(:, :)
! distance to vertices in 1D plot
real(8), allocatable :: dvp1d(:)
! plot vectors in lattice coordinates for 1D plot
real(8), allocatable :: vplp1d(:, :)
! distance to points in 1D plot
real(8), allocatable :: dpp1d(:)
! corner vectors of 2D plot in lattice coordinates
real(8)::vclp2d(3, 3)
! grid sizes of 2D plot
integer::np2d(2)
! corner vectors of 3D plot in lattice coordinates
real(8)::vclp3d(3, 4)
! grid sizes of 3D plot
integer::np3d(3)
! number of states for plotting Fermi surface
!replaced by inputstructureinteger::nstfsp
end module

