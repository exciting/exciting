
#include "maxdefinitions.inc"
module mod_SHT
!-----------------------------------------------------!
!     spherical harmonic transform (SHT) matrices     !
!-----------------------------------------------------!
! real backward SHT matrix for lmaxapw
real(8), allocatable :: rbshtapw(:, :)
! real forward SHT matrix for lmmaxapw
real(8), allocatable :: rfshtapw(:, :)
! real backward SHT matrix for lmaxvr
real(8), allocatable :: rbshtvr(:, :)
! real forward SHT matrix for lmaxvr
real(8), allocatable :: rfshtvr(:, :)
! complex backward SHT matrix for lmaxapw
complex(8), allocatable :: zbshtapw(:, :)
! complex forward SHT matrix for lmaxapw
complex(8), allocatable :: zfshtapw(:, :)
! complex backward SHT matrix for lmaxvr
complex(8), allocatable :: zbshtvr(:, :)
! complex forward SHT matrix for lmaxvr
complex(8), allocatable :: zfshtvr(:, :)
end module
