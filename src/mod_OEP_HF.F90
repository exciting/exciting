
#include "maxdefinitions.inc"
module mod_OEP_HF
!----------------------------------------!
!     OEP and Hartree-Fock variables     !
!----------------------------------------!
! maximum number of core states over all species
integer::ncrmax
! maximum number of OEP iterations
integer::maxitoep
! initial value and scaling factors for OEP step size
real(8)::tauoep(3)
! magnitude of the OEP residual
real(8)::resoep
! kinetic matrix elements
complex(8), allocatable :: kinmatc(:, :, :)
! complex versions of the exchange potential and field
complex(8), allocatable :: zvxmt(:, :, :)
complex(8), allocatable :: zvxir(:)
complex(8), allocatable :: zbxmt(:, :, :, :)
complex(8), allocatable :: zbxir(:, :)
end module
