

#include "maxdefinitions.inc"
module mod_Gvector
!--------------------------------!
!     G-vector set variables     !
!--------------------------------!
! G-vector cut-off for interstitial potential and density
real(8)::gmaxvr
! G-vector grid sizes
integer::ngrid(3)
! total number of G-vectors
integer::ngrtot
! integer grid intervals for each direction
integer::intgv(3, 2)
! number of G-vectors with G < gmaxvr
integer::ngvec
! G-vector integer coordinates
integer, allocatable :: ivg(:, :)
! map from integer grid to G-vector array
integer, allocatable :: ivgig(:, :, :)
! map from G-vector array to FFT array
integer, allocatable :: igfft(:)
! G-vectors in Cartesian coordinates
real(8), allocatable :: vgc(:, :)
! length of G-vectors
real(8), allocatable :: gc(:)
! spherical harmonics of the G-vectors
complex(8), allocatable :: ylmg(:, :)
! structure factor for the G-vectors
complex(8), allocatable :: sfacg(:, :)
! G-space characteristic function: 0 inside the muffin-tins and 1 outside
complex(8), allocatable :: cfunig(:)
! real-space characteristic function: 0 inside the muffin-tins and 1 outside
real(8), allocatable :: cfunir(:)
! damping coefficient for characteristic function
real(8)::cfdamp
end module
